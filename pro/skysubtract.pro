;+
; NAME:
;   skysubtract
;
; PURPOSE:
;   Sky-subtract an image and modify the variance
;
; CALLING SEQUENCE:
;   skysubtract, obj, objivar, plugsort, wset, objsub, objsubivar, $
;    [iskies= , fibermask=, nord=, upper=, lower=, maxiter=, pixelmask= ]
;
; INPUTS:
;   obj        - Image
;   objivar    - Inverse variance for OBJ
;   plugsort   - Plugmap structure trimmed to one element per fiber
;   wset       - Wavelength solution
;
; OPTIONAL KEYWORDS:
;   fibermask  - Mask of 0 for good fibers and non-zero for bad fibers [NFIBER]
;   pixelmask  - Mask of 0 for good pixels [NPIX,NFIBER]
;
; PARAMETERS FOR SLATEC_SPLINEFIT (for supersky fit):
;   nord       -
;   upper      -
;   lower      -
;   maxiter    -
;
; OUTPUTS:
;   objsub     - Image (OBJ) after sky-subtraction
;   objsubivar - Inverse variance (OBJIVAR) after sky-subtraction
;
; OPTIONAL OUTPUTS:
;   iskies=    - array of good sky fibers
;
; COMMENTS:
;   Construct a "supersky" spectrum by spline-fitting the (good) sky fibers,
;   resampling this at every wavelength in the extracted image, then
;   subtracting.  We then measure the variance of the sky-subtracted sky
;   fibers in terms of chi^2.  Wherever chi^2 > 1, we increase the variance
;   of all fibers such that chi^2=1 in the sky fibers.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_iterstat
;   djs_oplot
;   djs_oploterr
;   djs_plot
;   slatec_splinefit()
;   slatec_bvalue()
;
; REVISION HISTORY:
;   16-Sep-1999  Written by S. Burles, APO
;   30-Dec-1999  Modified by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro skysubtract, obj, objivar, plugsort, wset, objsub, objsubivar, $
   iskies=iskies, fibermask=fibermask, nord=nord, upper=upper, $
   lower=lower, maxiter=maxiter, pixelmask=pixelmask

   if (size(obj, /n_dimen) NE 2) then message, 'OBJIVAR is not 2-D'
   if (size(objivar, /n_dimen) NE 2) then message, 'OBJIVAR is not 2-D'

   dims = size(obj, /dimens)
   ncol = dims[0]
   nrow = dims[1]

   if (n_elements(fibermask) NE nrow) then fibermask = bytarr(nrow) 

   if ((size(plugsort, /dimens))[0] NE nrow) then $
    message, 'PLUGMAP does not have same size as nrow'

   if ( (size(wset.coeff, /dimens))[1] NE nrow) then $
    message, 'WSET does not have same size as nrow'

   ;----------
   ; Solve for wavelength of each pixel

   traceset2xy, wset, pixnorm, wave

   ;----------
   ; Find sky fibers
     
   iskies = where(plugsort.objtype EQ 'SKY' AND plugsort.fiberid GT 0 AND $
     (fibermask EQ 0), nskies)
   splog, 'Number of sky fibers = ', nskies
   if (nskies EQ 0) then message, 'No sky fibers in PLUGMAP'

   skywave = wave[*,iskies]
   skyflux = obj[*,iskies]
   skyivar = objivar[*,iskies]

   ;----------
   ; Sort sky points by wavelengths

   isort = sort(skywave)
   skywave = skywave[isort]
   skyflux = skyflux[isort]
   skyivar = skyivar[isort]

   ;----------
   ; Compute "supersky" with a spline fit
   ; Use the EVERYN parameter to space the spline points according to
   ; the density of data points.

   ; Return BKPT below for use later
   fullbkpt = slatec_splinefit(skywave, skyflux, coeff, invvar=skyivar, $
    nord=nord, upper=upper, lower=lower, maxiter=maxiter, $
    /eachgroup, everyn=nskies, bkpt=bkpt)

   ;----------
   ; Sky-subtract the entire image

   fullfit = slatec_bvalu(wave, fullbkpt, coeff) 
;   objsub = obj - fullfit * (objivar GT 0.0) ; No need to do this.


   objsub = obj - fullfit

   if (N_elements(pixelmask) EQ N_elements(obj)) then begin
       badsky = where(fullfit GT 10.0 * abs(objsub))   ;?? Does this make sense
       if (badsky[0] NE -1) then pixelmask[badsky] = $
                             pixelmask[badsky] OR pixelmask_bits('SKYLEVEL')
   endif

   ;----------
   ; Now attempt to model variance with residuals on sky fibers.
   ; This is difficult since variance has noise, so only do this if there
   ; are at least 3 sky fibers.

   if (nskies GE 3) then begin

      skyfit = (fullfit[*,iskies])[isort]
;      skyfit  = slatec_bvalu(skywave, fullbkpt, coeff) ; Same thing, more clear
      skychi2 = (skyflux-skyfit)^2 * skyivar

      ; Bin according to the break points used in the supersky fit.

      nbin = N_elements(bkpt) - 1
      relwave = fltarr(nbin)
      relchi2 = fltarr(nbin)
      for ibin=0, nbin-1 do begin
         ; Locate data points in this bin, exluding masked pixels
         ii = where(skywave GE bkpt[ibin] AND skywave LT bkpt[ibin+1] $
          AND skyivar GT 0, nn)

         if (nn GT 2) then begin
            ; Find the mean wavelength for these points
            relwave[ibin] = total(skywave[ii]) / nn

            ; Find the mean relative chi^2, assuming gaussian statistics.
            ; But this evaluation is wrecked by any outliers.
;            relchi2[ibin] = total(skychi2[ii]) / (nn-1)

            ; The following evaluation looks at the 67th percentile of the
            ; points, which is much more robust.
            pos67 = ceil(2.*nn/3.) - 1
            tmpchi2 = skychi2[ii]
            relchi2[ibin] = tmpchi2[ (sort(tmpchi2))[pos67] ]

            ; Schlegel counter of bin number...
            print, format='("Bin ",i4," of ",i4,a1,$)', $
             ibin, nbin, string(13b)

         endif
      endfor

      ; Trim to only those bins where we set RELCHI2
      ii = where(relwave NE 0, nbin)
      relwave = relwave[ii]
 
      ; Spline fit RELCHI2, only for the benefit of getting a smooth function
      ; Also, force the fit to always be >= 1, such that we never reduce the
      ; formal errors.
      fullbkpt = slatec_splinefit(relwave, relchi2, coeff, $
       maxiter=maxiter, upper=30, lower=30, everyn=2, nord=3)
      relchi2fit = slatec_bvalu(wave, fullbkpt, coeff)
      relchi2fit = relchi2fit > 1 ; Never let drop below 1

      splog, 'Median sky-residual chi2 = ', median(relchi2)
      splog, 'Max sky-residual chi2 = ', max(relchi2)

   endif else begin

      splog, 'WARNING: Too few sky fibers to model sky-sub variance'
      relchi2fit = 1

   endelse

   ;----------
   ; Modify OBJSUBIVAR with the relative variance

   objsubivar = objivar / relchi2fit

   ; Reselect the values of SKYIVAR from OBJSUBIVAR
;   skyivar = (objsubivar[*,iskies])[*]
;   skyivar = skyivar[isort]

   ;----------
   ; If any pixels on the image are outside of the wavelength range
   ; covered by the "supersky", then the formal errors are infinite
   ; for those pixels after skysubtraction.  Set SKYSUBIVAR=0.

   ii = where(skyivar GT 0, ni) ; Note that SKYWAVE is already sorted
   iout = where(wave LT skywave[ii[0]] OR wave GT skywave[ii[ni-1]])
   if (iout[0] NE -1) then objsubivar[iout] = 0.0

   return
end
;------------------------------------------------------------------------------
