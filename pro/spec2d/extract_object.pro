;+
; NAME:
;   extract_object
;
; PURPOSE:
;
;   Performs all object extraction tasks
;      0) Locate bright fibers, and test image background
;      1) 3 step Optimal extraction
;      2) Tweak to sky lines
;      3) Sky subtraction
;      4) Flux calibration
;      5) Telluric correction
;
; CALLING SEQUENCE:
;   extract_object, outname, objhdr, image, invvar, plugsort, wset, $
;    xarc, lambda, xtrace, fflat, fibermask, color=, $
;    [ widthset=, dispset=, skylinefile=, plottitle= ]
;
; INPUTS:
;   outname    - Name of outputs FITS file
;   objhdr     - Header cards from object image
;   image      - Object image [nx,ny]
;   invvar     - Inverse Variance of object image [nx,ny]
;   plugsort   - Plugmap structure for [ntrace] spectra
;   wset       - Wavelength solution from arc line spectra
;   xarc       - centroids measured in arc line spectra
;   lambda     - air wavelengths corresponding to xarc
;   xtrace     - spatial traces from flat field
;   fflat      - 1d flat field vectors
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   color      - ???
;   widthset   - ???
;   dispset    - ???
;   skylinefile- ???
;
; REQUIRED KEYWORDS:
;   color      - camera color (red or blue)
;
; OPTIONAL KEYWORDS:
;   plottitle  - Prefix for titles in QA plots.
;
; OUTPUTS:
;   A fits file is output in outname, which contains
;      FLOAT flambda [NX,NTRACE]
;      FLOAT flambda_invvar [NX,NTRACE]
;      LONG finalmask [NX,NTRACE]
;      STRUCT vacuum wavelengths
;      STRUCT wavelength dispersion
;      STRUCT plugmap [NTRACE]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   divideflat
;   djs_median()
;   djs_oplot
;   djs_plot
;   extract_boxcar()
;   extract_image
;   fibermask_bits()
;   calcscatimage()
;   fitvacset()
;   fluxcorr()
;   heliocentric()
;   locateskylines
;   mwrfits
;   pixelmask_bits()
;   qaplot_scatlight
;   qaplot_skydev
;   qaplot_skyline
;   qaplot_skyshift
;   qaplot_skysub
;   skysubtract
;   splog
;   sxaddpar
;   sxpar()
;   telluric_corr
;   tweaktrace
;
; INTERNAL SUPPORT ROUTINES:
;   find_whopping()
;
; REVISION HISTORY:
;   24-Jan-2000  Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------

function find_whopping, boxcar, thresh, whopct

    candidates = where(boxcar GT thresh, whopct)

    if (whopct LT 2) then return, candidates

    testc = [-20, candidates, n_elements(boxcar)+20]
    diff = testc[1:whopct+1] - testc[0:whopct]
    bottom = where(diff[0:whopct-1] NE 1, bottomct)
    top = where(diff[1:whopct] NE 1,topct)

    if (topct NE bottomct) then begin
      message, 'Bug introduced by Scott, look in find_whopping in extract_object.pro' ; ???
    endif

    whopping = lonarr(topct)
    for i=0, topct -1 do begin
       mmax = max(boxcar[bottom[i]:top[i]], place)
       whopping[i] = candidates[place + bottom[i]]
    endfor

    whopct = topct
    return, whopping
end

;------------------------------------------------------------------------------
pro extract_object, outname, objhdr, image, invvar, plugsort, wset, $
 xarc, lambda, xtrace, fflat, fibermask, color=color, $
 widthset=widthset, dispset=dispset, skylinefile=skylinefile, $
 plottitle=plottitle

   objname = strtrim(sxpar(objhdr,'OBJFILE'),2) 
   flavor  = strtrim(sxpar(objhdr,'FLAVOR'),2) 
 
   ;------------------
   ; Identify very bright objects
   ; Do a boxcar extraction, and look for fibers where the median
   ; counts are 10000 ADU per row.

   fextract = extract_boxcar(image, xtrace)

   scrunch = djs_median(fextract, 1) ; Find median counts/row in each fiber
   scrunch_sort = sort(scrunch)
   nfiber = (size(fextract,/dimens))[1]
   i5 = nfiber/20
   i95 = i5 * 19

   fullscrunch = djs_median(fextract) ; Find median counts/row in all fibers
   whopping = find_whopping(scrunch - fullscrunch, 10000.0, whopct)

   splog, 'Whopping fibers: ', whopping
   splog, 'Median counts in all fibers = ', fullscrunch
   splog, 'Number of bright fibers = ', whopct

   iskies = where(strtrim(plugsort.objtype,2) EQ 'SKY' $
    AND plugsort.fiberid GT 0 AND (fibermask EQ 0), nskies)

   if (nskies LT 2) then begin
      splog, 'ABORT: Only '+ string(nskies) + ' sky fibers found' 
      return
   endif 

   skymedian = djs_median(scrunch[iskies])
   splog, 'Sky fiber median '+string(skymedian)
;   if (skymedian GT 2000) then begin
;      splog, 'ABORT: Median sky flux is brighter than 2000 e-'
;      return
;   endif

   splog, '5% and 95% count levels ', scrunch[scrunch_sort[i5]], $
                                      scrunch[scrunch_sort[i95]]

   if (whopct GT 20) then begin
      splog, 'WARNING: Disable whopping terms ' + objname
      whopping = -1
      whopct = 0
   endif

   ;------------------------------------------------------------
   ;  Check for bad pixels within 3 pixels of trace

   badcheck = extract_boxcar((invvar LE 0), xtrace, radius=2.5)
   badplace = where(badcheck GT 0)

   nx = (size(fextract,/dim))[0] 
   ny = (size(fextract,/dim))[1] 
   pixelmask = lonarr(nx,ny)

   if (badplace[0] NE -1) then pixelmask[badplace] = $
                pixelmask[badplace] OR pixelmask_bits('NEARBADPIXEL')

   
   ;-----------------------------------------------------------------------
   ;  This is a kludge to fix first and last column
   ;-----------------------------------------------------------------------
   image[0,*] = image[0,*]*0.7
   image[2047,*] = image[2047,*]*0.7

   ;
   ;  First we should attempt to shift trace to object flexure
   ;

   xnow = match_trace(image, invvar, xtrace)
   bestlag = median(xnow-xtrace)

   splog, 'Shifting traces by match_trace ', bestlag

   if (abs(bestlag) GT 1.0) then begin
      splog, 'WARNING: pixel shift is large!'
   endif

   highrej = 5  ; just for first extraction steps
   lowrej = 5  ; just for first extraction steps
                ; We need to check npoly with new scattered light backgrounds
   npoly = 8 ; maybe more structure, lots of structure
   nrow = (size(image))[2]
   yrow = lindgen(nrow) 
   nfirst = n_elements(yrow)
   proftype = 3 ; Gaussian + Cubic

   splog, 'Extracting frame '+objname+' with 3 step process'

   traceset2xy, widthset, xx, sigma2
   ntrace = (size(sigma2,/dimens))[1]
   wfixed = [1,1]
   nterms = n_elements(wfixed)


   splog, 'Step 1: Initial Object extraction'

   extract_image, image, invvar, xnow, sigma2, tempflux, tempfluxivar, $
       proftype=proftype, wfixed=wfixed, yrow=yrow, $
       highrej=highrej, lowrej=lowrej, npoly=npoly, whopping=whopping, $
       ansimage=ansimage, chisq=firstchisq, ymodel=ym, /relative

     ; (2) Calculate scattered light

   splog, 'Step 2: Just find scattered light image'
   scatfit = calcscatimage(ansimage[ntrace*nterms:*,*], yrow)

   qaplot_scatlight, scatfit, yrow, $
    wset=wset, xcen=xtrace, fibermask=fibermask, $
    title=plottitle+'Scattered Light on '+objname

   ; (3) Calculate halo image
   splog, 'Step 3: Calculate Halo Image'
   smooth = smooth_halo(ym, wset)

   ;-----------------------------------------------------------------------
   ;  Now, subtract halo image and do final extraction with all rows
   ;-----------------------------------------------------------------------
   ; (4) Second and final extraction
   splog, 'Step 4: Final Object extraction'

   highrej = 4
   lowrej = 4
   wfixed = [1,1]
   nterms = n_elements(wfixed)
   reject = [0.2,0.6,0.6]


   extract_image, (image - smooth), invvar, xnow, sigma2, flux, fluxivar, $
    proftype=proftype, wfixed=wfixed, ansimage=ansimage, $
    highrej=highrej, lowrej=lowrej, npoly=npoly, whopping=whopping, $
    chisq=chisq, ymodel=ym, pixelmask=pixelmask, reject=reject, /relative

   ;----------------------------------------------------------------------
   ; can we find cosmic rays by looking for outlandish ansimage ratios?
   ;
   ; a = where(ansimage[lindgen(ntrace)*nterms, *] LT $
   ;           (-2*ansimage[lindgen(ntrace)*nterms+1, *])

   ;------------------
   ; QA chisq plot for fit calculated in extract image (make QAPLOT ???)

   xaxis = indgen(N_elements(chisq)) + 1
   djs_plot, xaxis, chisq, $
    xrange=[0,N_elements(chisq)], xstyle=1, $
    xtitle='Row number',  ytitle = '\chi^2', $
    title=plottitle+'Extraction chi^2 for '+objname

   djs_oplot, yrow, firstchisq[yrow], color='green'

   xyouts, 100, 0.05*!y.crange[0]+0.95*!y.crange[1], $
            'BLACK = Final chisq extraction'
   xyouts, 100, 0.08*!y.crange[0]+0.92*!y.crange[1], $
            'GREEN = Initial chisq extraction'

   ;------------------
   ; Flat-field the extracted object fibers with the global flat

   divideflat, flux, fluxivar, fflat, fibermask=fibermask

   lowflat = where(fflat LT 0.5)
   if (lowflat[0] NE -1) then pixelmask[lowflat] = $
                pixelmask[lowflat] OR pixelmask_bits('LOWFLAT')

   ;------------------
   ; Look for pixels where scattered light is dominating

   scatteredpix = where(extract_boxcar(scatfit, xnow, radius=2.0) GT flux)
   if (scatteredpix[0] NE -1) then pixelmask[scatteredpix] = $
                 pixelmask[scatteredpix] + pixelmask_bits('SCATTEREDLIGHT')

   ;------------------
   ; Tweak up the wavelength solution to agree with the sky lines.
   ; xshet contains polynomial coefficients to shift arc to sky line frame.

   locateskylines, skylinefile, flux, fluxivar, wset, $
    xarc, arcshift=arcshift, $
    xsky=xsky, skywaves=skywaves, skyshift=skyshift

   qaplot_skyshift, wset, xsky, skywaves, skyshift, $
    title=plottitle+'Sky Line Deviations for '+objname

   if (NOT keyword_set(arcshift)) then $
    splog, 'WARNING: Cannot shift to sky lines'

   ;------------------
   ; Apply heliocentric correction
   ; Correct LAMBDA, which is used to shift to vacuum wavelengths.

   helio=0.0
   ra = sxpar(objhdr,'RA')
   dec = sxpar(objhdr,'DEC')
   tai = sxpar(objhdr,'TAI')
   if (size(ra, /tname) NE 'INT' AND size(dec, /tname) NE 'INT' AND  $
    size(tai, /tname) NE 'INT') then begin
      helio = heliocentric(ra, dec, tai=tai)
      splog, 'Heliocentric correction = ', helio, ' km/s'
      sxaddpar, objhdr, 'HELIO_RV', helio, $
       ' Heliocentric correction (added to velocities)'
   endif else begin
      splog, 'WARNING: Header info not present to compute heliocentric correcoin'
   endelse

   ;------------------
   ; Shift to skylines and fit to vacuum wavelengths

   vacset = fitvacset(xarc, lambda, wset, arcshift, helio=helio)

   qaplot_skydev, flux, fluxivar, vacset, plugsort, color, $
    title=plottitle+objname

   sxaddpar, objhdr, 'VACUUM', 'T', 'Wavelengths are in vacuum'

   ;------------------
   ; Sky-subtract

   skystruct = skysubtract(flux, fluxivar, plugsort, vacset, $
    skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, $
    fibermask=fibermask, upper=3.0, lower=3.0, $
    relchi2struct=relchi2struct)

   ;-----------------------------------------------------------
   ;  Bad sky fibers???
   ;

   badskyfiber = where(djs_median(skysub[*,iskies]^2 * $
                skysubivar[*,iskies], 1) GT 2.0)               
   if badskyfiber[0] NE -1 then begin
       fibermask[iskies[badskyfiber]] = fibermask[iskies[badskyfiber]] OR $
          fibermask_bits('BADSKYFIBER')
       splog, 'WARNING: Calling Skysubtract again, masked skyfibers',$
            string(iskies[badskyfiber])
       skystruct = skysubtract(flux, fluxivar, plugsort, vacset, $
          skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, $
          fibermask=fibermask, upper=3.0, lower=3.0, $
          relchi2struct=relchi2struct)
   endif
 
   ;
   ; Sky-subtract again, this time with dispset (PSF subtraction)
   ; 

;   skystruct_psf = skysubtract(flux, fluxivar, plugsort, vacset, $
;     skysubpsf, skysubpsfivar, iskies=iskies, pixelmask=pixelmask, $
;     fibermask=fibermask, upper=3.0, lower=3.0, dispset=dispset)

   qaplot_skysub, flux, fluxivar, skysub, skysubivar, $
    vacset, iskies, title=plottitle+objname

   ;------------------
   ; QA for 2 skylines in the blue (specify vacuum wavelengths below)

   if (color EQ 'blue') then begin
      qaplot_skyline, 4359.5, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, iskies, fibermask=fibermask, dwave=4.0, $
       title=plottitle+objname
      qaplot_skyline, 5578.9, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, iskies, fibermask=fibermask, dwave=5.0, $
       title=plottitle+objname
   endif

   ;------------------
   ; QA for 2 skylines in the red (specify vacuum wavelengths below)

   if (color EQ 'red') then begin
      qaplot_skyline, 7343.0, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, iskies, fibermask=fibermask, dwave=7.0, $
       title=plottitle+objname
      qaplot_skyline, 8888.3, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, iskies, fibermask=fibermask, dwave=7.0, $
       title=plottitle+objname
   endif

   ;------------------------------------------
   ; Save the sky-subtracted flux values as is, and now modify flambda.

   flambda = skysub
   flambdaivar = skysubivar

   ;------------------------------------------
   ; Telluric correction called for 'red' side

   if (color EQ 'red')  then begin

      ; The following commented-out code essentially reproduces the
      ; telluric-correction code implemented from Oct 99 to Aug 00.
      ; However, if this is implemented, it should be done **after**
      ; the flux-calibration step below.
;      tellbands1 = { TELLBAND1, $
;       twave1: 6607., twave2: 8318., $
;       cwave1: 6607., cwave2: 8313. }
;      tellbands2 = { TELLBAND1, $
;       twave1: 8710., twave2: 9333., $
;       cwave1: 8710., cwave2: 9333. }
;      tellbands = [tellbands1, tellbands2]
;
;      ttt = telluric_corr(flambda, flambdaivar, vacset, plugsort, $
;       fibermask=fibermask, tellbands=tellbands, pixspace=100, $
;       upper=5, lower=5, /dospline, $
;       plottitle=plottitle+'Telluric correction for '+objname)

      telluricfactor = telluric_corr(flambda, flambdaivar, vacset, plugsort, $
       fibermask=fibermask, $
       plottitle=plottitle+'Telluric correction for '+objname)

      divideflat, flambda, flambdaivar, telluricfactor, minval=0.1

   endif

   ;------------------
   ; Flux calibrate to spectrophoto_std fibers

; Remove this code???
;   fluxfactor = fluxcorr(flambda, flambdaivar, vacset, plugsort, $
;    color=color, lower=3.0, upper=3.0, fibermask=fibermask)
;
;   minfluxfactor = median(fluxfactor) * 0.01
;   divideflat, flambda, flambdaivar, fluxfactor, minval=minfluxfactor

   ;----------
   ; Interpolate over masked pixels, just for aesthetic purposes

   flambda = djs_maskinterp(flambda, flambdaivar EQ 0, /const, iaxis=0 )

   ;----------
   ; Combine FIBERMASK and PIXELMASK to FINALMASK

   finalmask = pixelmask
   for itrace=0, ntrace-1 do $
    finalmask[*,itrace] = finalmask[*,itrace] OR fibermask[itrace]

   ;----------
   ; Add keywords to object header

   sxaddpar, objhdr, 'VERS2D', idlspec2d_version(), $
    'Version of idlspec2d for 2D reduction', after='VERSREAD'
   if (keyword_set(osigma)) then $
    sxaddpar, objhdr, 'OSIGMA',  sigma, $
     'Original guess at spatial sigma in pix'
   sxaddpar, objhdr, 'LOWREJ', lowrej, 'Extraction: low rejection'
   sxaddpar, objhdr, 'HIGHREJ', highrej, 'Extraction: high rejection'
   sxaddpar, objhdr, 'SCATPOLY', npoly, 'Extraction: Order of scattered light polynomial'
   sxaddpar, objhdr, 'PROFTYPE', proftype, 'Extraction profile: 1=Gaussian'
   sxaddpar, objhdr, 'NFITPOLY', nterms, 'Extraction: Number of parameters in each profile'

   ;----------
   ; Write extracted, lambda-calibrated, sky-subtracted spectra to disk

   mwrfits, flambda, outname, objhdr, /create
   mwrfits, flambdaivar, outname
   mwrfits, finalmask, outname
   mwrfits, vacset, outname
   mwrfits, dispset, outname
   mwrfits, plugsort, outname

; Save sky, fluxfactor, telluricfactor???
;   mwrfits, fluxfactor, outname
   if (color EQ 'red') then mwrfits, telluricfactor, outname

   return
end
;------------------------------------------------------------------------------
