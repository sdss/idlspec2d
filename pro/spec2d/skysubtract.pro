;+
; NAME:
;   skysubtract
;
; PURPOSE:
;   Sky-subtract an image and modify the variance
;
; CALLING SEQUENCE:
;   skystruct = skysubtract(obj, objivar, plugsort, wset, objsub, objsubivar, $
;    [ iskies= , fibermask=, nord=, upper=, lower=, maxiter=, pixelmask=, $
;    dispset=, /novariance, relchi2struct=, npoly=, tai= ])
;
; INPUTS:
;   obj        - Image
;   objivar    - Inverse variance for OBJ
;   plugsort   - Plugmap structure trimmed to one element per fiber
;   wset       - Wavelength solution
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   pixelmask  - Mask of 0 for good pixels [NPIX,NFIBER]
;   dispset    - Dispersion trace-set; if present, then solve for the
;                super-sky vector using a variable PSF described by this
;                structure.
;   novariance - Set keyword to prevent variance correction for sky residuals
;   relchi2struct- Structure containing information of chi^2 fitting
;   npoly      - Polynomial order of 2d fit.
;   tai        - TAI of plate exposure, if supplied a linear airmass correction
;                  is applied
;   nbkpt      - Number of bkpts to use for full sky spectrum fit.
;                This gives us the freedom to use less points for the blue side.
;
; PARAMETERS FOR SLATEC_SPLINEFIT (for supersky fit):
;   nord       -
;   upper      -
;   lower      -
;   maxiter    -
;
; OUTPUTS:
;   skystruct  - structure containing sorted sky wavelengths,flux,fluxivar
;                      +  bkpts and coeffs from fitting
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
;   pixelmask_bits()
;   bspline_iterfit()
;   bspline_valu()
;   splog
;   redmonster
;   tai2airmass()
;   traceset2xy
;
; REVISION HISTORY:
;   16-Sep-1999  Written by S. Burles, APO
;   30-Dec-1999  Modified by D. Schlegel, Princeton
;    4-Oct-2000  Changed to bspline_iterfit 
;   26-Jul-2001  Implemented Higher-order fits, if dispset is set it fits
;                  against fiber number (I know this doesn't make sense)
;                Also upped the tolerance on rejection by 5 times!
;                Too many pixels were masked for high order fit!
;-
;------------------------------------------------------------------------------

function skysubtract, obj, objivar, plugsort, wset, objsub, objsubivar, $
   iskies=iskies, fibermask=fibermask, nord=nord, upper=upper, $
   lower=lower, maxiter=maxiter, pixelmask=pixelmask, $
   dispset=dispset, npoly=npoly, relchi2struct=relchi2struct, $
   novariance=novariance, tai=tai, nbkpt=nbkpt

   if (size(obj, /n_dimen) NE 2) then message, 'OBJIVAR is not 2-D'
   if (size(objivar, /n_dimen) NE 2) then message, 'OBJIVAR is not 2-D'

   dims = size(obj, /dimens)
   ncol = dims[0]
   nrow = dims[1]

   if NOT keyword_set(upper) then upper = 10.0
   if NOT keyword_set(lower) then lower = 10.0
   if NOT keyword_set(nbkpt) then nbkpt = ncol
   thresh = 3.0

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
     
   iskies = where(strtrim(plugsort.objtype, 2)  EQ 'SKY' $
    AND plugsort.fiberid GT 0 AND (fibermask EQ 0), nskies)
   splog, 'Number of sky fibers = ', nskies
   if (nskies EQ 0) then begin
      splog, 'ABORT: No unmasked sky fibers in PLUGMAP'
      return, 0
   endif

   if NOT keyword_set(tai) then airmass = replicate(1.0, nrow) $
    else airmass = float(tai2airmass(plugsort.ra, plugsort.dec, tai=tai))

   minairmass = min(airmass, max=maxairmass)
   splog, (maxairmass GT 2.5) ? 'WARNING: ' : '', $
    'Airmass range = ', minairmass, maxairmass

   airmass_correction = replicate(1.0,ncol) # airmass

   skywave = wave[*,iskies]
   skyflux = obj[*,iskies]
   skyivar = objivar[*,iskies]

   divideflat, skyflux, invvar=skyivar, airmass_correction[*,iskies]

   ;----------
   ; Mask any sky pixels where LOWFLAT or NEARBADPIXEL are set.

   if (keyword_set(pixelmask)) then begin
      qbad = ((pixelmask[*,iskies] AND pixelmask_bits('LOWFLAT')) NE 0) $
           OR ((pixelmask[*,iskies] AND pixelmask_bits('NEARBADPIXEL')) NE 0)
      ibad = where(qbad, nbad)
      splog, 'Discarding ', float(nbad)/n_elements(skywave), $
       ' (fractional) of the sky pixels as bad'
      if (nbad GT 0) then skyivar[ibad] = 0
   endif

   ;----------
   ; Sort sky points by wavelengths

   isort = sort(skywave)
   skywave = skywave[isort]
   skyflux = skyflux[isort]
   skyivar = skyivar[isort]

   ;----------
   ; Compute "supersky" with a spline fit.
   ; Use the EVERYN parameter to space the spline points according to
   ; the density of data points.

   bkpt = 0
   firstset = bspline_iterfit(skywave, skyflux, invvar=skyivar, $
    nord=nord, upper=upper, lower=lower, maxiter=maxiter, $
    /eachgroup, everyn=2*nskies/3, bkpt=bkpt, outmask=outmask, yfit=skyfit)

   if (NOT keyword_set(skyfit)) then begin
      splog, 'ABORT: Fit sky is all zeros'
      return, 0
   endif

   ;----------
   ; This code has been added to slightly reduce the number of
   ; breakpoints, from 3100 to NCOL, where we now space them such
   ; that there is an equal amount of S/N between each.
   ; This also gives better behavior near the boundaries.

   igood = where(skyivar GT 0 AND outmask NE 0, ngoodpix)
   minwave  = min(skywave[igood])
   snsqrt = sqrt((skyfit * sqrt(skyivar) > 0))
   ipos = where(snsqrt GT 0, npos)
   gkern = gauss_kernel(2.0*nskies)

   if (npos GT n_elements(gkern) AND nbkpt GE 16) then begin
      snsqrt[ipos] = convol(snsqrt[ipos], gauss_kernel(2.0*nskies))

      snsum = snsqrt
      for i=1L, n_elements(snsqrt)-1 do snsum[i] = snsum[i] + snsum[i-1]
      place = long(snsum/max(snsum)*nbkpt) < (nbkpt-2)
      newbkpt = uniq(place)
      bkpt = skywave[newbkpt]
      bkpt[0] = minwave
      bkpt[7:nbkpt-8] = (smooth(skywave[newbkpt],7))[7:nbkpt-8]

      ; Why the factor of 10 below ???
      newnbk = n_elements(bkpt)
      lowdiff = (bkpt[1] - bkpt[0]) * 10.0
      highdiff = (bkpt[newnbk-1] - bkpt[newnbk-2]) * 10.0
      fullbkpt = [ bkpt[0] - lowdiff, bkpt, max(bkpt)+highdiff ]
      for i=2, nord-1 do $
       fullbkpt= [ fullbkpt[0] - lowdiff, fullbkpt, max(fullbkpt)+ highdiff ]
   endif

   ;----------
   ; Re-do the fit with the new break points.
 
   sset = bspline_iterfit(skywave, skyflux, invvar=skyivar, $
    nord=nord, upper=upper*1.5, lower=lower*1.5, maxiter=maxiter, $
    /eachgroup, fullbkpt=fullbkpt, yfit=skyfit2, outmask=outmask)
 
   fullfit = bspline_valu(wave, sset) 

   if (keyword_set(dispset)) then begin

      if (NOT keyword_set(npoly)) then npoly = 4L

      ; SIGMA is smooth fit to widths of arclines ???
      traceset2xy, dispset, pixnorm, sigma
      sigma = sigma - 1.0
      skysigma = (sigma[*,iskies])[isort]

;      sset = bspline_iterfit(skywave, skyflux, invvar=skyivar*outmask, $
;       nord=nord, upper=upper, lower=lower, maxiter=maxiter, x2=skysigma, $
;       npoly = npoly, /eachgroup, fullbkpt=fullbkpt, $
;       yfit=skyfit4)
;
;      fullfit = bspline_valu(wave, sset, x2=sigma) 
;
      ssetorig = sset
      fullx2 = replicate(1.0,ncol) # findgen(nrow)
      x2     = (fullx2[*,iskies])[isort]

      sset = bspline_iterfit(skywave, skyflux, invvar=skyivar*outmask, $
       nord=nord, upper=upper, lower=lower, maxiter=maxiter, $
       npoly = npoly, /eachgroup, fullbkpt=fullbkpt, $
       yfit=skyfit3, x2=x2, xmin=0., xmax=nrow)

      fullfit = bspline_valu(wave, sset, x2=fullx2) 

;
;      The commented out code below fits the residuals "flux-fit" to
;       higher order polynomial terms.  Not used for now
;
;       yfit=skyfit4, x2=x2, xmin=0., xmax=nrow, funcname='poly1', outmask=o)

;      sset4 = bspline_iterfit(skywave, skyflux-skyfit2, invvar=skyivar*outmask, $
;       nord=nord, upper=upper, lower=lower, maxiter=maxiter, $
;       npoly = npoly-1, /eachgroup, fullbkpt=fullbkpt, $
;       yfit=skyfit4, x2=x2, funcname='poly1', outmask=o)
;
;     fullfit = bspline_valu(wave, sset4, x2=fullx2) + fullfit
;       yfit=skyfit4, x2=skysigma, funcname='poly1',outmask=o)
;
;     fullfit = bspline_valu(wave, sset4, x2=sigma) + fullfit
   endif

   ;----------
   ; Sky-subtract the entire image

;   objsub = obj - fullfit * (objivar GT 0.0) ; No need to do this.
   objsub = obj - fullfit * airmass_correction

   ;----------
   ; Fit to sky variance (not inverse variance)

   posvar = where(skyivar GT 0)
   if (posvar[0] NE -1) then begin
      skyvariance = 1.0 / skyivar[posvar]
      skyvarset = bspline_iterfit(skywave[posvar], skyvariance, $
        invvar=skyivar[posvar], nord=nord, upper=upper, lower=lower, $
        maxiter=maxiter, /eachgroup, bkpt=bkpt)

;      skyvarbkpt = slatec_splinefit(skywave[posvar], skyvariance, skyvarcoeff, $
;       invvar=skyivar[posvar], nord=nord, upper=upper, lower=lower, $
;       maxiter=maxiter, /eachgroup, bkpt=bkpt)
      skyvarfit = bspline_valu(wave, skyvarset) * airmass_correction

   endif

   ;----------
   ; Store "super" sky information in a structure
   ; We can't name it, because it could change size each time

   skystruct = create_struct( $
    'ISKIES', iskies, $
    'WAVE', skywave, $
    'FLUX', skyflux, $
    'INVVAR', skyivar, $
    'FULLBKPT', sset.fullbkpt, $
    'COEFFS', sset.coeff)

   ;----------
   ; Now attempt to model variance with residuals on sky fibers.
   ; This is difficult since variance has noise, so only do this if there
   ; are at least 3 sky fibers.

   if (nskies GE 3 AND ngoodpix GT ncol $
    AND NOT keyword_set(novariance)) then begin

      skyfit = (fullfit[*,iskies])[isort]
;      skyfit  = slatec_bvalu(skywave, fullbkpt, coeff) ; Same thing, more clear
      skychi2 = (skyflux-skyfit)^2 * skyivar

      ; Bin according to the break points used in the supersky fit.

      nbin = N_elements(bkpt) - 1
      relwave = fltarr(nbin)
      relchi2 = fltarr(nbin)
      for ibin=0, nbin-1 do begin
         ; Locate data points in this bin, excluding masked pixels
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

            ; Burles counter of bin number...
            print, format='("Bin ",i4," of ",i4,a1,$)', $
             ibin, nbin, string(13b)

         endif
      endfor

      ; Trim to only those bins where we set RELCHI2
      ; ?? Do we need to trim relchi2 also, to have matching size arrays??

      ii = where(relwave NE 0, nbin)
      relwave = relwave[ii]
      relchi2 = relchi2[ii]
 
      ;----------
      ; Spline fit RELCHI2, only for the benefit of getting a smooth function
      ; Also, force the fit to always be >= 1, such that we never reduce the
      ; formal errors.

      relchi2set = bspline_iterfit(relwave, relchi2, nord=3, $
        upper=30, lower=30, maxiter=maxiter, everyn=2)

      relchi2fit = bspline_valu(wave, relchi2set) > 1

;      fullbkpt = slatec_splinefit(relwave, relchi2, coeff, $
;       maxiter=maxiter, upper=30, lower=30, everyn=2, nord=3)
;      relchi2fit = slatec_bvalu(wave, fullbkpt, coeff)
;      relchi2fit = relchi2fit > 1 ; Never let drop below 1

      splog, 'Median sky-residual chi2 = ', median(relchi2)
      splog, 'Max sky-residual chi2 = ', max(relchi2)

      ;----------
      ; Store Relative Chi2 information in a structure
      ; Add in other information we want to write to disk??
      
      relchi2struct = create_struct( $
          'WAVE', relwave, $
          'CHI2', relchi2, $
          'FULLBKPT', relchi2set.fullbkpt, $
          'COEFF', relchi2set.coeff)

   endif else begin

      splog, 'WARNING: Too few sky fibers or pixels to model sky-sub variance'
      relchi2fit = 1

   endelse

   ;----------
   ; Modify OBJSUBIVAR with the relative variance

   objsubivar = objivar
   if (keyword_set(skyvarfit) AND n_elements(relchi2fit) GT 1) then begin
      posvar = where(objivar GT 0)
      if (posvar[0] NE -1) then begin
        objvar = 1.0 / objivar[posvar]
        objsubivar[posvar] = $
         1.0 / (objvar + ((relchi2fit[posvar] > 1)-1.0)*skyvarfit[posvar])
      endif
   endif 

   ; Reselect the values of SKYIVAR from OBJSUBIVAR
;   skyivar = (objsubivar[*,iskies])[*]
;   skyivar = skyivar[isort]

   ;----------
   ; Set the BADSKYCHI mask bit at any wavelength where the relative chi^2
   ; for the sky fibers is greater than 4.
   ; Also, look for the Red Monster (any bad region contiguous in wavelength)

   if (keyword_set(relchi2)) then begin
      if (keyword_set(pixelmask)) then $
       pixelmask = pixelmask OR pixelmask_bits('BADSKYCHI') $
        * (relchi2fit GT thresh)

      redmonster, relwave, relchi2, wave, pixelmask=pixelmask
   endif

   ;----------
   ; If any pixels on the image are outside of the wavelength range
   ; covered by the "supersky", then the formal errors are infinite
   ; for those pixels after skysubtraction.  Set the mask bit 'NOSKY'
   ; and set SKYSUBIVAR=0.
   ;
   ; Also, set the mask bit 'BRIGHTSKY' for any pixels where the sky
   ; level is more than the (sky-subtracted) object flux + 10 * error,
   ; and where the sky is greater than 1.25 times a median sky.
   ; Grow this mask by 1 neighboring pixel in each direction.

   ii = where(skyivar GT 0, ni) ; Note that SKYWAVE is already sorted
   iout = where(wave LT skywave[ii[0]] OR wave GT skywave[ii[ni-1]])
   if (iout[0] NE -1) then objsubivar[iout] = 0.0

   if (keyword_set(pixelmask)) then begin
      ; Compute a median sky vector for each fiber
      medsky = 0 * fullfit
      for irow=0, nrow-1 do $
       medsky[*,irow] = djs_median(fullfit[*,irow], width=99, $
        boundary='reflect')

      qbright = (objsubivar NE 0) $
       AND (fullfit GT objsub + 10.0 / sqrt(objsubivar + (objsubivar EQ 0))) $
       AND (fullfit GT 1.25 * medsky)
      qbright = convol(qbright, [1,1,1], /center, /edge_truncate)
      ibright = where(qbright)
      if (ibright[0] NE -1) then $
       pixelmask[ibright] = pixelmask[ibright] OR pixelmask_bits('BRIGHTSKY')

      if (iout[0] NE -1) then $
       pixelmask[iout] = pixelmask[iout] OR pixelmask_bits('NOSKY')
   endif

   return, skystruct
end
;------------------------------------------------------------------------------
