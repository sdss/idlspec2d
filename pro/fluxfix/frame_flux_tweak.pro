;+
; NAME:
;   smearcorr
;
; PURPOSE:
;   Create smear correction vectors
;
; CALLING SEQUENCE:
;   smearcorr, bsmearfile, rsmearfile, bcalibfile, rcalibfile, $
;              bscifile, rscifile, corrfile, noplot = noplot
;
; INPUTS:
;   bsmearfile - spFrame FITS file chosen as blue smear image
;   rsmearfile - spFrame FITS file chosen as red smear image
;   bcalibfile - spFrame FITS file chosen as blue calib image
;   rcalibfile - spFrame FITS file chosen as red calib image
;   bscifile   - spFrame FITS file(s) chosen as blue science image
;   rscifile   - spFrame FITS file(s) chosen as red science image
;   corrfile   - FITS file to output flux correction vectors.
;
; OPTIONAL KEYWORDS:
;   noplot     - toggles plotting
;
; OUTPUTS:  Files containing the polynomial coefficients of the fit to 
;           each fiber (red+blue combined). The naming convention is 
;           spFluxcorr-pppp-mmmm-s.fits where, pppp is the plateid, 
;           mmmmm is the mjd and s is the spectrograph ID.  Plots of
;           the smear correction vectors are also created for each
;           frame+spectrograph.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;   smearcorr is used to calculate and write to file a low order
;   polynomial function which registers the flux in the science exposures
;   to the flux in a fiducial exposure.  Point sources have their SED
;   tied to a "smear" exposure if one is available.  Galaxies have their
;   SED tied to the "calib" image which is typically the science frame
;   with the highest S/N.
;
;   The smear and calib images have different flux levels due to the very
;   different exposure times and observing strategies employed.  These 
;   images need to be put on the same scale, because the flux calibration
;   derived from the standard stars will be used to calibrate the galaxies.
;   This is accomplished by scaling the smear such that the median of
;   the smear image matches the median of the calib image. 
;
;   Unlike the flux correction, the smear correction varies from fiber
;   to fiber.  After a coarse rebinning along the dispersion direction,
;   the fiducial image is divided by the science image, and each fiber 
;   is fit with a low order function.  For the highest S/N fibers (typically 
;   10\%, including spectrophoto stds) the full 3rd order fit is done.
;
;   For medium S/N fibers, only a one parameter scaling of the best 
;   spectrophoto correction is produced.  For points sources the median
;   spectrophoto correction is created from the standard stars.  For 
;   galaxies the median correction is computed from all of the high S/N
;   galaxies.
;
;   For the lowest S/N fibers, the median spectrophoto correction is used 
;   (as above for stars and galaxies) but without any rescaling.  
;
; EXAMPLES:
;
; BUGS:
;
;  Blue wavelength region is hardwired: b1 = findgen(60)*4.0e-3 + 3.568
;  Red wavelength region is hardwired : r1 = findgen(54)*4.0e-3 + 3.756
;  Order of polynomial is hardwired:  3
; 
;  If no high S/N spectrophoto stars exist then the low and medium S/N point 
;  sources have their smear correction set to 1.0.  Galaxies may have a very 
;  different mean correction than this, resulting in a big star/galaxy
;  discrepancy in the calibration.  In general, this routine will work 
;  badly at low S/N.
;
;  Should probably use more mask bits to flag the 5(!!) different S/N cases
;  (high, medium-point-source, medium-galaxy, low-point-source, low-galaxy)
;
;  A difference in flats between frames could produce a different offsets
;  between the smear & science frames in the blue and red.  Since we fit
;  them together this could cause problems for the derived smear correction.
;  We attempt to corect for this by removing the median offset of the
;  blue & red frames in narrow passband (5000-5800, 6200-7500) that should be
;  relatively unaffected by differential refraction.  More testing is needed 
;  to know if this helps or just adds another source of noise.
;
; PROCEDURES CALLED:
;   djs_iterstat
;   mrdfits()
;   mwrfits
;   pixelmask_bits()
;   traceset2xy
;   xy2traceset
;
; INTERNAL SUPPORT ROUTINES:
;   median_rebin():  Used to rebin spectra in large wavelength blocks
;                    passed in parameter range 
;   rebin_exposure(): Used to rebin an entire frame (flux, ivar, mask, wave)
;                     to lower resolution  
;   smear_coeff(): Used to simplify the process of getting the correct
;                  polynomial coefficients for various S/N conditions.
;
; REVISION HISTORY:
;   17-Oct-2000  Formerly fluxcorr_new -- written by S. Burles
;   09-Oct-2002  Revised by C. Tremonti to calibrate point sources to the 
;                smear and galaxies to the calib images. Added rebin_exposure()
;                and smear_coeff() to streamline. 
;-
;------------------------------------------------------------------------------
function median_rebin, loglam, flux, ivar, color, mask=mask, sigrej=sigrej, $
         outwave = outwave, sn = sn, quality = quality

   if NOT keyword_set(sigrej) then sigrej = 20.0

   ;---------
   ; Define new wavelength sampling

   if color eq 'r' then begin
     w1 = findgen(54)*4.0e-3 + 3.756
     w2 = findgen(54)*4.0e-3 + w1[1]
   endif 
   if color eq 'b' then begin
     w1 = findgen(60)*4.0e-3 + 3.568
     w2 = findgen(60)*4.0e-3 + w1[1]
   endif
   range = transpose([[w1],[w2]])
   outwave = djs_median(range,1)

   ;--------------------------
   ; Set up vectors for output
   nr = (size(range))[2]
   ntrace = (size(flux))[2]
   fit = fltarr(nr, ntrace)
   mask = fltarr(nr, ntrace)

   ;--------------------------
   ; Do iterative rejection on each bin of each fiber
   for itrace=0,ntrace-1 do begin
     for irange = 0, nr -1 do begin
        inside = where(loglam[*,itrace] GE range[0,irange] $
                 AND loglam[*,itrace] LT range[1,irange], ninside)
        if ninside GT 0 then begin
           good = where(ivar[inside,itrace] GT 0, ngood)
       
           if ngood GT 1 then begin 
              djs_iterstat, flux[inside[good],itrace], median=md, sigma=sig
              fit[irange, itrace] = md
              if ngood GT 0.5 * ninside then $
                   mask[irange, itrace]  = total(ivar[inside[good],itrace])
              sn = sig * sqrt(mask[irange, itrace])
              if sn GT sigrej OR sn LE 0 then mask[irange, itrace] = 0.0
           endif
        endif
     endfor
   endfor

   sn = djs_median(flux * sqrt(ivar),1)

   badval = ( fibermask_bits('NOPLUG')   OR fibermask_bits('BADTRACE') $
       OR fibermask_bits('BADFLAT')  OR fibermask_bits('BADARC')   $
       OR fibermask_bits('NEARWHOPPER') )

   quality = (mask[0,*] AND badval)

   return, fit
end

;------------------------------------------------------------------------------
pro frame_flux_tweak, bloglam, rloglam, bflux, rflux, bivar, rivar, $
    best_exp, fibertag, corrfile

   expid = fibertag[uniq(fibertag.expid)].expid
   splog, 'Best Exposure =', best_exp
   indx = where(fibertag.expid eq best_exp)

   ;----------
   ; Create red & blue smoothed best images

   bfit = median_rebin(bloglam[*,indx], bflux[*,indx], bivar[*,indx], $
          'b', outwave = bwave, mask = bmask, sn = bsn, quality = bquality)

   rfit = median_rebin(rloglam[*,indx], rflux[*,indx], rivar[*,indx], $
          'r', outwave = rwave, mask = rmask, sn = rsn, quality = rquality)

   ;---------------------------------------------------------
   ; Now put blue and red together
    
   bestflux = [bfit,rfit]
   bestivar = [bmask, rmask]
   bestsnmed = transpose([[bsn],[rsn]])
   qbest =  (bquality OR rquality)

   npix = (size(bestflux,/dimens))[0]
   nfiber = (size(bestflux,/dimens))[1]
   wave = [bwave,rwave] # replicate(1,nfiber)

   ; --------------------------------------------------------
   ;  Create basic mask
    
   bestmask = lonarr(nfiber)

   ;----------
   ; Loop through the science images

   nfiles = n_elements(corrfile)

   for ifile=0, nfiles-1 do begin

      indx = where(fibertag.expid eq expid[ifile])
      splog, 'Fluxing science image =', expid[ifile]

      ;------------
      ; Create coeffients of fit and set them to a zero-order polynomaial
      ; with an amplitude of 1.0
 
      fitimg = bestflux*0.0 + 1.0
      xy2traceset, wave, fitimg, corrset, ncoeff=ncoeff, /silent
      corrset.coeff[0,*] = 1.0
      corrset.coeff[1:*,*] = 0.0
      thismask  = bestmask 

      ;--------------
      ; If the science image and best image are the same,
      ; force their ratio to be unity.
      if (expid[ifile] EQ best_exp) then begin
        mwrfits, corrset, corrfile[ifile], /create
      ; mwrfits, thismask, corrfile[ifile]
        continue
      endif 

      ;-------------
      ; Create rebinned blue and red science images

; Need to figure out the best way to set mask bits!!!

      bscifit = median_rebin(bloglam[*,indx], bflux[*,indx], bivar[*,indx], $
                'b', mask = bscimask, sn = bscisn, quality = bsciquality)

      rscifit = median_rebin(rloglam[*,indx], rflux[*,indx], rivar[*,indx], $
                'r', mask = rscimask, sn = rscisn, quality = rsciquality)
       
      ;----------------
      ; Combine red and blue

      sciflux = [bscifit,rscifit]
      sciivar = [bscimask,rscimask]
      scisnmed = transpose([[bscisn],[rscisn]])
      qsci = lonarr(nfiber)
      qsci = qsci OR bsciquality
      qsci = qsci OR rsciquality

      ;-------------------------------------------------------------
      ; Determine the S/N of the expsure and act accrodingly
      ; poly 3:  High S/N => 3 order polynomial fit (parabola)
      ; poly 2:  Medium S/N => 2 order polynomial fit (line)
      ; poly 1:  Low S/N => 1 order polynomial fit (zero point shift)

      fiber_coeff = lonarr(nfiber) + 1

      poly3 = where(scisnmed[0,*] GT 2.5 AND scisnmed[1,*] GT 5.0  $
               AND  bestsnmed[0,*] GT 2.5 AND  bestsnmed[1,*] GT 5.0 $
               AND  qsci EQ 0 AND qbest EQ 0, npoly3) 
      if poly3[0] NE -1 then fiber_coeff[poly3] = 3 

      poly2 = where(scisnmed[0,*] GT 1.0 AND scisnmed[1,*] GT 2.0  $
               AND  bestsnmed[0,*] GT 1.0 AND bestsnmed[1,*] GT 2.0 $
               AND  qsci EQ 0 AND qbest EQ 0 AND fiber_coeff NE 3, npoly2) 
      if poly2[0] NE -1 then fiber_coeff[poly2] = 2

      poly1 = where(fiber_coeff eq 1, npoly1)

      ;-----------------
      ; Calculate polynomial coeffs for the 3 cases
   
      if npoly3 gt 0 then begin
        xy2traceset, wave[*,poly3], bestflux[*,poly3], polyset, $
            invvar=bestivar[*,poly3], ncoeff=3, inputfunc=sciflux[*,poly3], $
            lower = 3, upper = 3
        corrset.coeff[*,poly3] = polyset.coeff
      endif

      if npoly2 gt 0 then begin
        xy2traceset, wave[*,poly2], bestflux[*,poly2], polyset, $
            invvar=bestivar[*,poly2], ncoeff=2, inputfunc=sciflux[*,poly2], $
            lower = 3, upper = 3
        corrset.coeff[0,poly2] = polyset.coeff[0,*]
        corrset.coeff[1,poly2] = polyset.coeff[1,*]
      endif
      if npoly1 gt 0 then begin
        xy2traceset, wave[*,poly1], bestflux[*,poly1], polyset, $
            invvar=bestivar[*,poly1], ncoeff=1, inputfunc=sciflux[*,poly1], $
            lower = 3, upper = 3
        corrset.coeff[0,poly1] = polyset.coeff
      endif

      ;----------
      ; Identify sky fibers from the plug map   
      sky = where(strmatch(fibertag[indx].objtype, '*SKY*'))
      corrset.coeff[0,sky] = 1.0
      corrset.coeff[1,sky] = 0.0
      corrset.coeff[2,sky] = 0.0

      ;---------------
      ; Reject any bad smear corrections and replace with low S/N solution
      ; solution.  Bad vectors are those with fit coefficients outside
      ; of pre-set boundaries.  The boundaries on coeffs 1 & 2 increase 
      ; with increasing spread in the coeff0 values -- this should help
      ; to accomodate really bad plates.

      ; Normalize coefficients 
      coef0 = corrset.coeff[0,*]  
      coef1 = corrset.coeff[1,*] / coef0 
      coef2 = corrset.coeff[2,*] / coef0 
 
      meanclip, coef0, coef0mean, coef0sig, clipsig = 5
      coef0 = coef0 / coef0mean
      coef0sig = (coef0sig / coef0mean) > 0.10
    
      ; Adjust boundary values for plate quality
      coef1bound = 0.35 * coef0sig / 0.10 
      coef2bound = 0.25 * coef0sig / 0.10

      isbad = coef0 LT  (1 - 5 * coef0sig) OR coef0 GT (1 + 5 * coef0sig) OR $
              coef1 LT -coef1bound OR coef1 GT coef1bound OR $
              coef2 LT -coef2bound OR coef2 GT coef2bound

      bad = where(isbad, nbad)
      if bad[0] NE -1 then begin
        splog, 'Warning: Large deviations in flux correction '
        splog, 'Warning: Replacing with zero point shift:', $
        string(bad + 1)

        xy2traceset, wave[*,bad], bestflux[*,bad], polyset, $
            invvar=bestivar[*,bad], ncoeff=1, inputfunc=sciflux[*,bad], $
            lower = 3, upper = 3
        corrset.coeff[0,bad] = polyset.coeff
        corrset.coeff[1,bad] = 0.0
        corrset.coeff[2,bad] = 0.0
      endif

      ;---------------
      ; Set maskbits and append to end of corrfile
         
      ;thismask = thismask OR (fibersn EQ 2) * pixelmask_bits('SMEARMEDSN')
      ;thismask = thismask OR (fibersn EQ 3) * pixelmask_bits('SMEARHIGHSN')

      ;------------
      ; Write out as FITS

      mwrfits, corrset, corrfile[ifile], /create
      ;mwrfits, thismask, corrfile[ifile]

      ;------------
      ; Plot smear correction vectors

      if keyword_set(noplot) then continue

      traceset2xy, corrset, wave, corrimage

      !P.MULTI = [0, 1, 2]
      std = where(strmatch(fibertag[indx].objtype, '*_STD*'), nstd)
      for ii = 0, nstd - 1 do begin
        istd = std[ii]
        plot, 10.0^wave[*,istd], corrimage[*,istd], yr=[0.5, 1.5], /nodata, $
                ytitle='Best Frame Flux / Science Frame Flux'
        oplot, 10.0^wave[*,istd], bestflux[*,istd]/sciflux[*,istd], psym=6, $
               syms=0.5, thick=3
        djs_oplot, 10.0^wave[*,istd], corrimage[*,istd], color='red', thick=3
      endfor
      !P.MULTI = 0

      djs_plot, 10.0^wave, corrimage, /nodata, yr=[0.5, 1.5], $
                xr=[min(10.0^wave)-100,max(10.0^wave)+100], $
                /xstyle, /ystyle, xtitle='\lambda [A]', $
                ytitle='Best Frame Flux / Science Frame Flux', $
                title= 'Flux Correction: ' + corrfile[ifile]

      hipts = where(fiber_coeff eq 3, nhipts)
      lowpts = where(fiber_coeff ne 3, nlowpts)
      std = where(strmatch(fibertag[indx].objtype, '*_STD*'), nstd)

      for iobj=0, nlowpts -1 do $
        djs_oplot, 10.0^wave[*,lowpts[iobj]], corrimage[*,lowpts[iobj]], $
                   color='blue', nsum=10
      for iobj=0, nhipts -1 do $
        djs_oplot, 10.0^wave[*,hipts[iobj]], corrimage[*,hipts[iobj]], $
                   color='green', nsum=10

      for iobj=0, nstd -1 do $
        djs_oplot, 10.0^wave[*,std[iobj]], corrimage[*,std[iobj]], $
                   color='red', nsum=10, thick=2

      djs_xyouts, 0.18, 0.95, 'High S/N Sources', color = 'green', /norm
      djs_xyouts, 0.18, 0.91, 'Low S/N Sources', color = 'blue', /norm
      djs_xyouts, 0.18, 0.87, 'Standard Stars', color = 'red', /norm

   endfor
   return
end
