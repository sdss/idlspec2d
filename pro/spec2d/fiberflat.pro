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
;    /dospline, /nonorm, plottitle=, badflatfracthresh= ])
;
; INPUTS:
;   flux       - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   fluxivar   - Inverse variance map for FLUX.
;   wset       - Wavelength solution
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
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
;   nonorm     - Do not normalize the fluxes in FFLAT by the super-flat.
;   superflatset-Bspline set to reconstruct superflat
;   badflatfracthresh - the fraction of bad flat pixels to decide if
;                       the flat is bad
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
;   bspline_valu()
;   bspline_iterfit()
;   splog
;   superflat()
;   traceset2xy
;   xy2traceset
;
; REVISION HISTORY:
;   14-Oct-1999  Written by D. Schlegel, APO
;    3-Oct-2000  Changed over to IDL bspline routines
;-
;------------------------------------------------------------------------------
function scale_master, flux, master_flux, loglam, minlam, maxlam, outside=outside, npix = npix, red=red
    idx = sort(loglam)
    sll = loglam[idx]
    sf  = flux[idx]
    mf  = master_flux[idx]
    inside = where((sll GE minlam) and (sll LE maxlam), cti)
    if cti gt 10 then nipix = 10 else nipix =cti
    sf = sf[inside]
    mf = mf[inside]
    npix = 10

    if keyword_set(red) then begin
        outside = where(sll GT maxlam, cto)
        if cto eq 0 then return, -1
        
        idx_f = where(finite(sf), ct)
        if ct eq 0 then return, -1
        sf = sf[idx_f]
        idx_f = where(finite(mf), ct)
        if ct eq 0 then return, -1
        mf = mf[idx_f]
        if n_elements(mf) lt 10 then npix = n_elements(mf)
        
        scale = mean(sf[-1*nipix:*])/mean(mf[-1*nipix:*])
        if idx[0] ne 0 then outside = where(loglam GT maxlam,cto)
    endif else begin
        outside = where(loglam LT minlam, cto)
        if cto eq 0 then return, -1

        idx_f = where(finite(sf), ct)
        if ct eq 0 then return, -1
        sf = sf[idx_f]
        idx_f = where(finite(mf), ct)
        if ct eq 0 then return, -1
        mf = mf[idx_f]
        if n_elements(mf) lt 10 then npix =  n_elements(mf)
        scale = mean(sf[0:nipix-1])/mean(mf[0:nipix-1])
        if idx[0] ne 0 then outside = where(loglam LT minlam,cto)
    endelse
    
    if not finite(scale) then scale = -1
    return, scale

end

;------------------------------------------------------------------------------
function fiberflat, flux, fluxivar, wset, fibermask=fibermask, $
 minval=minval, ncoeff=ncoeff, pixspace=pixspace, nord=nord, $
 lower=lower, upper=upper, dospline=dospline, plottitle=plottitle, $
 nonorm=nonorm, superflatset=superflatset, superflat_minval=superflat_minval, $
 badflatfracthresh=badflatfracthresh, configuration =configuration, $
 master_flat=master_flat, pad_blue=pad_blue, pad_red=pad_red

   pad_blue = 0
   pad_red = 0

   dims = size(flux, /dimens)
   ny = dims[0]
   ntrace = dims[1]
   fflat = fltarr(ny,ntrace)
   fratio_pad = fltarr(ny,ntrace)

   if (NOT keyword_set(minval)) then minval = 0.03 * median(flux)
   if (N_elements(pixspace) EQ 0) then pixspace = 10
   if (N_elements(ncoeff) EQ 0) then ncoeff = 3
   if (N_elements(nord) EQ 0) then nord = 4
   if (N_elements(lower) EQ 0) then lower = 10.0
   if (N_elements(upper) EQ 0) then upper = 10.0
   if (N_elements(fibermask) NE ntrace) then fibermask = bytarr(ntrace)
   if (not keyword_set(badflatfracthresh)) then badflatfracthresh=0.7

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

   superflatset = superflat(flux, fluxivar, wset, $
    fibermask=fibermask, minval=minval, lower=3.0, upper=3.0, $
    medval=medval, title=plottitle)

   if (NOT keyword_set(superflatset)) then begin
      splog, 'WARNING: Spline fit failed' 
      return, -1
   endif

   fit2  = bspline_valu(loglam, superflatset)

   if keyword_set(master_flat) then begin
        master_ff = mrdfits(master_flat,0,hdr, status=status, /SILENT)
        if status ne 0 then master_flat = 0
   endif
   ;----------

   if (keyword_set(dospline)) then begin
      if keyword_set(master_flat) then begin
        splog,'Padding with FFlat from '+FILE_BASENAME(master_flat)
        master_fflat = mrdfits(master_flat,0, /SILENT)
        master_flux = mrdfits(master_flat,'FLUX', /SILENT)
        master_ivar = mrdfits(master_flat,'IVAR', /SILENT)
        master_wset = mrdfits(master_flat,'WSET', /SILENT)
        traceset2xy, master_wset, xx, master_loglam
        superflat_minval = fltarr(ntrace)
        
        ; Always select the same break points in log-wavelength for all fibers
        nbkpts = fix(ny / pixspace) + 2
        bkpt = findgen(nbkpts) * (max(loglam) - min(loglam)) / (nbkpts-1) $
                + min(loglam)

        for i=0, ntrace-1 do begin
            print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])

            ; The following should work for either ascending or descending
            ; wavelengths since BKPT is always sorted ascending.
 
            indx = where(fluxivar[*,i] GT 0.0 AND flux[*,i] GT minval $
                      AND fit2[*,i] GT 0.0, ct)
            minlam = min(loglam[indx,i])
            maxlam = max(loglam[indx,i])
            istart = (where(bkpt GT minlam))[0]
            istart = (istart - 1) > 0
            iend = (where(bkpt GT maxlam))[0]
            if (iend EQ -1) then iend = nbkpts-1
            
            ; this interpolate the master flat to this flat, it should be close, but might not be exact
            IF (MIN(loglam[*,i]) GT MAX(loglam[*,i])) THEN begin
                tloglam = REVERSE(loglam[*,i])
            endif else tloglam = loglam[*,i]
            IF (MIN(master_loglam[*,i]) GT MAX(master_loglam[*,i])) THEN $
                master_loglam[*,i] = REVERSE(master_loglam[*,i])

            master_flux[*,i] = INTERPOL(master_flux[*,i], master_loglam[*,i], tloglam)
            master_ivar[*,i] = INTERPOL(master_ivar[*,i], master_loglam[*,i], tloglam)
        
            IF (MIN(loglam[*,i]) GT MAX(loglam[*,i])) then begin
                master_flux[*,i] = REVERSE(master_flux[*,i])
                master_ivar[*,i] = REVERSE(master_ivar[*,i])
            endif
        
            scale_b = scale_master(flux[*,i], master_flux[*,i], loglam[*,i], $
                                    minlam, maxlam, outside=outside_b)
            if scale_b ne -1 then begin
                if (total(master_fflat[*,i]) ne 0) and (total(flux[*,i]) gt 0) then begin
                    flux[outside_b,i] = master_flux[outside_b,i]*scale_b
                    fluxivar[outside_b,i] = master_ivar[outside_b,i]/scale_b^2
                    if keyword_set(pad_blue) then pad_blue = [pad_blue, 10.0^minlam] $
                    else pad_blue = [10.0^minlam]
                endif
            endif
            
            scale_r = scale_master(flux[*,i], master_flux[*,i], loglam[*,i], $
                                    minlam, maxlam, outside=outside_r, /red)
            if scale_r ne -1 then begin
                if (total(master_fflat[*,i]) ne 0) and (total(flux[*,i]) gt 0) then begin
                    flux[outside_r,i] = master_flux[outside_r,i]*scale_r
                    fluxivar[outside_r,i] = master_ivar[outside_r,i]/scale_r^2
                    if keyword_set(pad_red) then pad_red = [pad_red, 10.0^maxlam] $
                    else pad_red = [10.0^maxlam]
                endif
            endif
            fratio_pad[*,i] = flux[*,i]/fit2
            
        endfor
        ffit2 = fit2
        if keyword_set(pad_red) then pad_red = mean(pad_red)
        if keyword_set(pad_blue) then pad_blue = mean(pad_blue)
        
        minval = configuration->spcalib_fiberflat_minval(flux)
        superflatset = superflat(flux, fluxivar, wset, $
                                 fibermask=fibermask, $
                                 minval=configuration->spcalib_fiberflat_minval(flux), $
                                 lower=3.0, upper=3.0, $
                                 medval=medval, title=plottitle+' padded')

        if (NOT keyword_set(superflatset)) then begin
            splog, 'WARNING: Spline fit failed'
            return, -1
        endif

        fit2  = bspline_valu(loglam, superflatset)
      endif
     
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
         if keyword_set(pad_red) or keyword_set(pad_blue) then begin
            indx = where(fluxivar[*,i] GT 0.0 $
                         AND fratio_pad[*,i] GT 0.5*median(fratio_pad[*,i]) $
                         AND fit2[*,i] GT 0.0, ct)
         endif else begin
            indx = where(fluxivar[*,i] GT 0.0 AND flux[*,i] GT minval $
                         AND fit2[*,i] GT 0.0, ct)
         endelse

         if (ct GE 5) then begin ; Require at least 5 data points

            ; The following should work for either ascending or descending
            ; wavelengths since BKPT is always sorted ascending.
 
            minlam = min(loglam[indx,i])
            maxlam = max(loglam[indx,i])
            istart = (where(bkpt GT minlam))[0]
            istart = (istart - 1) > 0
            iend = (where(bkpt GT maxlam))[0]
            if (iend EQ -1) then iend = nbkpts-1

            ratio  = flux[indx,i] / fit2[indx,i]
            ratioivar  = fluxivar[indx,i] * fit2[indx,i]^2

            ; Dispose of leading or trailing points with zero weight

            ratioset = bspline_iterfit(loglam[indx,i],ratio,invvar=ratioivar, $
             maxiter=maxiter, upper=upper, lower=lower, $
             groupsize=n_elements(indx), $
             nord=nord, bkpt=bkpt[istart:iend], requiren=2)

            inside = where(loglam[*,i] GE minlam AND loglam[*,i] LE maxlam)
            if inside[0] NE -1 then $
              fflat[inside,i] = bspline_valu(loglam[inside,i], ratioset)
            if keyword_set(superflat_minval) then superflat_minval[i] = min(fit2[inside,i],/NAN)
         endif else begin

            fflat[*,i] = 0

         endelse

      endfor
      if keyword_set(superflat_minval) then begin
            superflat_minval = mean(superflat_minval, /NAN)
            if not FINITE(superflat_minval) then superflat_minval = 0
            if keyword_set(superflat_minval) then superflat_minval= max([superflat_minval,0.005])
            splog, 'DEBUG', superflat_minval
      endif
   endif else begin

      ;------------------------------------------------------------------------
      ; LEGENDRE FIT TO FFLAT VECTORS
      ;------------------------------------------------------------------------

      ratimg = flux / fit2
      rativar = fluxivar * fit2^2

      ;----------
      ; Replace each flat-field vector with a cubic fit to that vector

      inmask = fluxivar GT 0 AND flux GT minval ; Mask bad pixels in the fit
      xy2traceset, loglam, ratimg, fset, func='legendre', ncoeff=ncoeff, $
       ; invvar=rativar, $ ; Weight all points equally instead
       maxiter=100, maxrej=1, /sticky, $
       inmask=inmask, outmask=xmask, yfit=fflat

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

   indx = where(medval LT 0.7 * globalmed $
    OR total(fflat GT 0,1) LT badflatfracthresh*ny)
   if (indx[0] NE -1) then $
    fibermask[indx] = fibermask[indx] OR fibermask_bits('BADFLAT')

   if (keyword_set(nonorm)) then return, fflat * fit2 $
    else return, fflat
end
;------------------------------------------------------------------------------
