;+
; NAME:
;   fiberflat
;
; PURPOSE:
;   Construct the flat-field vectors from an extracted flat-field image.
;
; CALLING SEQUENCE:
;   fflat = fiberflat( flux, fluxivar, wset, [ fibermask=fibermask, $
;    minval=, ncoeff=, pixspace=, /dospline, nord=, lower=, upper=,
;    /dospline, plottitle= ])
;
; INPUTS:
;   flux       - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   fluxivar   - Inverse variance map for FLUX.
;   wset       - Wavelength solution
;
; OPTIONAL KEYWORDS:
;   fibermask  - Mask of 0 for bad fibers and 1 for good fibers [NFIBER]
;   minval     - Minimum value to use in fits to flat-field vectors;
;                default to 3% of the median of FLUX.
;   ncoeff     - Number of coefficients used in constructing FFLAT;
;                default to 3 (cubic)
;   pixspace   - Approximate spacing in pixels for break points in the
;                spline fits to individual fibers; default to 10 pixels.
;   dospline   - If this keyword is set, then fit the flat-field vectors
;                to splines (using PIXSPACE) rather than to a Legendre
;                polynomial (using NCOEFF).
;                This is now what we use?
;   plottitle  - Title for QA plot; if not set, then do not plot.
;
; PARAMETERS FOR SLATEC_SPLINEFIT:
;   nord
;   lower
;   upper
;
; OUTPUTS:
;   fflat      - Array of flat-field flat-field vectors for each fiber
;                that remove relative flat-field variations as a function
;                of wavelength between fibers [Nrow, Ntrace]
;
; OPTIONAL OUTPUTS:
;   fibermask  - (Modified)
;
; COMMENTS:
;   The user should first "flat-field" the input array to take out
;   pixel-to-pixel variations.
;
;   The parameters for SLATEC_SPLINEFIT are only used when generating the
;   "superflat".
;
;   The 'BADFLAT' bit is set in FIBERMASK if the mean throughput for
;   a fiber is less than 0.7 times the median of all good-fiber throughputs.
;
;   In any given fiber, set FFLAT=0 wherever there are at least 5 contiguous
;   bad pixels.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   fibermask_bits()
;   slatec_bvalu()
;   slatec_splinefit()
;   splog
;   superflat
;   traceset2xy
;   xy2traceset
;
; REVISION HISTORY:
;   14-Oct-1999  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------

function fiberflat, flux, fluxivar, wset, fibermask=fibermask, $
 minval=minval, ncoeff=ncoeff, pixspace=pixspace, nord=nord, $
 lower=lower, upper=upper, dospline=dospline, plottitle=plottitle

   dims = size(flux, /dimens)
   ny = dims[0]
   ntrace = dims[1]
   fflat = fltarr(ny,ntrace)

   if (NOT keyword_set(minval)) then minval = 0.03 * median(flux)
   if (N_elements(pixspace) EQ 0) then pixspace = 10
   if (N_elements(ncoeff) EQ 0) then ncoeff = 3
   if (N_elements(nord) EQ 0) then nord = 4
   if (N_elements(lower) EQ 0) then lower = 10.0
   if (N_elements(upper) EQ 0) then upper = 10.0
   if (N_elements(fibermask) NE ntrace) then fibermask = bytarr(ntrace)

   igood = where(fibermask EQ 0, ngood)
   if (ngood EQ 0) then begin
     splog, 'WARNING: No good fibers according to FIBERMASK'
     return, -1
   endif 

   ;----------
   ; Compute the wavelengths for all flat vectors from the trace set

   traceset2xy, wset, xx, loglam

   ;----------
   ; Construct the "superflat" vector

   superflat, flux, fluxivar, wset, afullbkpt, acoeff, $
    fibermask=fibermask, minval=minval, lower=lower, upper=upper, $
    medval=medval, title=plottitle

   if (n_elements(allfullbkpt) EQ 1) then begin
      splog, 'WARNING: Spline fit failed' 
      return, -1
   endif

   fit2  = slatec_bvalu(loglam, afullbkpt, acoeff)

   ;----------

   if (keyword_set(dospline)) then begin

      ;------------------------------------------------------------------------
      ; SPLINE FIT TO FFLAT VECTORS
      ;------------------------------------------------------------------------

      ; Always select the same break points in log-wavelength for all fibers
      nbkpts = fix(ny / pixspace) + 2
      bkpt = findgen(nbkpts) * (max(loglam) - min(loglam)) / (nbkpts-1) $
       + min(loglam)

      for i=0, ntrace-1 do begin
         print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])

         ; Evaluate "superflat" spline fit at exactly the same wavelengths
         ; Let's divide out superflat first to make fitting smoother
         ; Larger breakpoint separations and less hassles 

         ; Locate only unmasked points
;         indx = where(fluxivar[*,i] GT 0.0 AND flux[*,i] GT minval $
;          AND fit2[*,i] GT minval, ct)
         indx = where(fluxivar[*,i] GT 0.0 AND flux[*,i] GT minval, ct)

         if (ct GE 5) then begin ; Require at least 5 data points

            ; The following should work for either ascending or descending
            ; wavelengths since BKPT is always sorted ascending.

            istart = (where(bkpt GT min(loglam[indx,i])))[0]
            istart = (istart - 1) > 0
            iend = (where(bkpt GT max(loglam[indx,i])))[0]
            if (iend EQ -1) then iend = nbkpts-1

            ratio = flux[indx,i] / fit2[indx,i]
            ratioivar = fluxivar[indx,i] * fit2[indx,i]^2

            ; Dispose of leading or trailing points with zero weight
            fullbkpt = slatec_splinefit(loglam[indx,i], ratio, coeff, $
             maxiter=maxiter, upper=upper, lower=lower, /eachgroup, $
             invvar=ratioivar, nord=nord, bkpt=bkpt[istart:iend], mask=mask)

            ; Evaluate spline fit to this fiber
            if (N_elements(fullbkpt) GT 1) then begin
               fflat[*,i] = slatec_bvalu(loglam[*,i], fullbkpt, coeff)
            endif else begin
               fflat[*,i] = 0
            endelse

            ; Replace leading or trailing masked points with the first or last
            ; unmasked value.
            ; SLATEC_BVALU already does the check below
            ;if (indx[0] NE 0) then $
            ; fit1[0:indx[0]-1] = fit1[indx[0]]
            ;if (indx[ct-1] NE ny-1) then $
            ; fit1[indx[ct-1]+1:ny-1] = fit1[indx[ct-1]]

         endif else begin

            fflat[*,i] = 0

         endelse

      endfor

   endif else begin

      ;------------------------------------------------------------------------
      ; LEGENDRE FIT TO FFLAT VECTORS
      ;------------------------------------------------------------------------

      ratimg = flux / fit2
      rativar = fluxivar * fit2^2

      ;----------
      ; Replace each flat-field vector with a cubic fit to that vector

      xmask = fluxivar GT 0 AND flux GT minval ; Mask bad pixels in the fit
      xy2traceset, loglam, ratimg, fset, func='legendre', ncoeff=ncoeff, $
       ; invvar=rativar, $ ; Weight all points equally instead
        maxiter=10, /singlerej, xmask=xmask, yfit=fflat

      ;----------
      ; For flat vectors that are completely bad, replace with zeros.

      indx = where(total(xmask,1) EQ 0)
      if (indx[0] NE -1) then fflat[*,indx] = 0

   endelse

   ;----------
   ; Set FFLAT=0 only when there are at least 5 bad pixels in a row. 
   ; Smaller gaps should be OK with our spline-fitting across them.

   sz = 5
   for i=0, ntrace-1 do begin
      indx = where(smooth( (smooth((fluxivar[*,i] NE 0)*sz, sz) EQ 0)*sz, sz ))
      if (indx[0] NE -1) then fflat[indx,i] = 0
   endfor

   ;----------
   ; Check to see if there are fewer good fibers

   igood = where(fibermask EQ 0 AND total(fflat,1) GT 0, ngood)
   if (ngood EQ 0) then begin
      splog, 'WARNING: All flat fibers have been rejected!'
      return, -1
   endif
 
   ;----------
   ; Divide FFLAT by a global median of all (good) fibers

   globalmed = median([medval[igood]]) ; Global median for all vectors
   fflat = fflat / globalmed 

   junk = where(fflat LE 0, nz)
   splog, 'Number of fiberflat points LE 0 = ', nz

   ;----------
   ;  Set flatfield bit in FIBERMASK if needed

   indx = where(medval LT 0.7 * globalmed)
   if (indx[0] NE -1) then $
    fibermask[indx] = fibermask[indx] OR fibermask_bits('BADFLAT')

   return, fflat
end
;------------------------------------------------------------------------------
