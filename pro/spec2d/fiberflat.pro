;+
; NAME:
;   fiberflat
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   fflat = fiberflat( flux, fluxivar, wset, [ fibermask=fibermask, $
;    minval=, ncoeff=, pixspace=, /dospline, nord=, lower=, upper=, /dospline ] )
;
; INPUTS:
;   flux       - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   fluxivar   - Inverse variance map for FLUX.
;   wset       - Wavelength solution
;
; OPTIONAL KEYWORDS:
;   fibermask  - Mask of 0 for bad fibers and 1 for good fibers [NFIBER]
;   minval     - Minimum value to use in fits to flat-field vectors;
;                default to 0.03.
;   ncoeff     - Number of coefficients used in constructing FFLAT;
;                default to 3 (cubic)
;   pixspace   - Approximate spacing in pixels for break points in the
;                spline fits to individual fibers; default to 10 pixels.
;   dospline   - If this keyword is set, then fit the flat-field vectors
;                to splines (using PIXSPACE) rather than to a Legendre
;                polynomial (using NCOEFF).  This is **not** recommended.
;
; PARAMTERS FOR SLATEC_SPLINEFIT:
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
;
; COMMENTS:
;   The user should first "flat-field" the input array to take out
;   pixel-to-pixel variations.
;
;   The parameters for SLATEC_SPLINEFIT are only used when generating the
;   "superflat".
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   slatec_bvalu()
;   slatec_splinefit()
;   traceset2xy
;
; REVISION HISTORY:
;   14-Oct-1999  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------

function fiberflat, flux, fluxivar, wset, fibermask=fibermask, $
 minval=minval, ncoeff=ncoeff, pixspace=pixspace, dospline=dospline, $
 nord=nord, lower=lower, upper=upper

   dims = size(flux, /dimens)
   ny = dims[0]
   ntrace = dims[1]
   fflat = fltarr(ny,ntrace)

   if (NOT keyword_set(minval)) then minval = 0.03
   if (N_elements(pixspace) EQ 0) then pixspace = 10
   if (N_elements(ncoeff) EQ 0) then ncoeff = 3
   if (N_elements(nord) EQ 0) then nord = 4
   if (N_elements(lower) EQ 0) then lower = 10.0
   if (N_elements(upper) EQ 0) then upper = 10.0
   if (N_elements(fibermask) NE ntrace) then fibermask = bytarr(ntrace) + 1

   igood = where(fibermask NE 0, ngood)
   if (ngood EQ 0) then $
    message, 'No good fibers according to FIBERMASK'

;   ; If any flux points have zero or negative flux, set the weight to zero
;   ibad = where(flux LE 0)
;   if (ibad[0] NE -1) then fluxivar[ibad] = 0

   ; Compute the wavelengths for all flat vectors from the trace set
   traceset2xy, wset, xx, loglam

   ; Determine the range of wavelengths, [LOGMIN,LOGMAX] in common w/all fibers
   if (loglam[1,0] GT loglam[0,0]) then begin ; Ascending wavelengths
      logmin = max(loglam[0,igood])
      logmax = min(loglam[ny-1,igood])
   endif else begin ; Descending wavelengths
      logmin = max(loglam[ny-1,igood])
      logmax = min(loglam[0,igood])
   endelse

   ; Find the approximate scalings between all fibers
   ; Do this with a straight median value for all wavelengths in common
   qq = loglam GE logmin AND loglam LE logmax
   medval = fltarr(ntrace)
   for i=0, ntrace-1 do $
    medval[i] = median( flux[where(qq[*,i]),i] )
   izero = where(medval LE 0)
   if (izero[0] NE -1) then medval[izero] = 1.0

   ; Create a version of flux (and fluxivar) that has all fibers
   ; approximately scaled to have a median value of 1
   scalef = fltarr(ny,ntrace)
   scalefivar = fltarr(ny,ntrace)
   for i=0, ntrace-1 do $
    scalef[*,i] = flux[*,i] / medval[i]
   for i=0, ntrace-1 do $
    scalefivar[*,i] = fluxivar[*,i] * (medval[i])^2

   ; Create a "superflat" spectrum, analogous to the "supersky"
   splog, 'Creating superflat from ', ngood, ' fibers'
   isort = sort(loglam[*,igood])
   allwave = (loglam[*,igood])[isort]
   allflux = (scalef[*,igood])[isort]
   allivar = (scalefivar[*,igood])[isort]
   indx = where(allflux GT minval)
   if (indx[0] EQ -1) then $
    message, 'No points above MINVAL'
   afullbkpt = slatec_splinefit(allwave[indx], allflux[indx], acoeff, $
    maxiter=maxiter, upper=upper, lower=lower, $
    invvar=allivar[indx], nord=4, nbkpts=ny, mask=mask)

   fit2  = slatec_bvalu(loglam, afullbkpt, acoeff)

   ;------
   ; QA plot of superflat
   ; Plot sampled every 1 Ang

   wmin = fix(10^min(allwave))
   wmax = ceil(10^max(allwave))
   plot_lam = wmin + lindgen(wmax-wmin+1)
   plot_fit  = slatec_bvalu(alog10(plot_lam), afullbkpt, acoeff)
 
   djs_plot, plot_lam, plot_fit, xrange=[wmin,wmax], xstyle=1, $
    xtitle='\lambda [A]', ytitle='Normalized flux', $
    title='Superflat for ???'

   ; Overplot pixels masked from the fit
   ii = where(mask EQ 0)
   djs_oplot, 10^allwave[indx[ii]], allflux[indx[ii]], ps=3, color='red' 

   ;------

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
         indx = where(fluxivar[*,i] GT 0.0 AND flux[*,i] GT minval $
          AND fit2 GT minval, ct)

         if (ct GT 0) then begin

; Does the following work for either ascending or descending wavelengths???
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
            fflat[*,i] = slatec_bvalu(loglam[*,i], fullbkpt, coeff)

            ; Replace leading or trailing masked points with the first or last
            ; unmasked value
            ; bvalu already does the check below
            ;if (indx[0] NE 0) then fit1[0:indx[0]-1] = fit1[indx[0]]
            ;if (indx[ct-1] NE ny-1) then fit1[indx[ct-1]+1:ny-1] = fit1[indx[ct-1]]

         endif else begin

            fflat[*,i] = 1.0
            fibermask[i] = 0

         endelse

      endfor

   endif else begin

      ;------------------------------------------------------------------------
      ; LEGENDRE FIT TO FFLAT VECTORS
      ;------------------------------------------------------------------------

      ratimg = flux / fit2
      rativar = fluxivar * fit2^2

      ; Replace each flat-field vector with a cubic fit to that vector
      xmask = fluxivar GT 0 ; Mask bad pixels in the fit
      xy2traceset, loglam, ratimg, fset, func='legendre', ncoeff=ncoeff, $
       ; invvar=rativar, $ ; Weight all points equally instead
        maxiter=10, /singlerej, xmask=xmask, yfit=fflat

   endelse

   ; Check to see if fibermask has changed
   igood = where(fibermask NE 0, ngood)

   if (igood[0] EQ -1) then message,'All flat fibers have been !!'
 
   ; Divide fflat by a global median of all fibers
   globalmed = median(medval[igood]) ; Global median for all vectors
   fflat = fflat / globalmed 

   junk = where(fflat LE 0, nz)
   splog, 'Number of fiberflat points LE 0 = ', nz

; Take out radial dependence...???
;   fibermed = djs_median(fflat,1)
;   plot,fibermed,ps=1
;   if (keyword_set(plugmap)) then begin
;      ; Fitting out radial dependence
;      r2 = (plugmap.xFocal^2 + plugmap.yFocal^2)*1.0e-6  ; units m^2
;      plot,r2,fibermed,ps=1
;      radialcoeff = polyfitw(r2,fibermed,fibermask,1,radialfit)
;      for i=0,ntrace-1 do fflat[*,i] = fflat[*,i]/radialfit[i]
;      splog, radialcoeff
;   endif

   return, fflat
end
;------------------------------------------------------------------------------
