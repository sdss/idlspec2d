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
;   fitansimage()
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
   if (skymedian GT 2000) then begin
      splog, 'ABORT: Median sky flux is brighter than 2000 e-'
      return
   endif

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

   maxshift = 2.0 ; ??? Need this for MJD=51579
   dims = size(image,/dimens)
   ncol = dims[0]
   nrow = dims[1]
   skiptrace = 20L
   ysample = lindgen(nrow) # replicate(1, nfiber - 2*skiptrace)
   xsample = xtrace[*,skiptrace:nfiber - skiptrace - 1]
   bestlag = shift_trace(image, xsample, ysample, lagrange=1.0, lagstep=0.1)

   splog, 'Shifting traces by pixel shift of ', bestlag

   if (abs(bestlag) GT 1.0) then begin
      splog, 'WARNING: pixel shift is large!'
   endif

   xnow = xtrace + bestlag

   highrej = 5  ; just for first extraction steps
   lowrej = 5  ; just for first extraction steps
                ; We need to check npoly with new scattered light backgrounds
   npoly = 8 ; maybe more structure, lots of structure
   skiprow = 8
   yrow = lindgen(nrow/skiprow) * skiprow + skiprow/2
   nfirst = n_elements(yrow)

   ;-----------------------------------------------------------------------
   ;  The fork in the road:
   ;    If we have widths from widthset, then just extract
   ;    otherwise determine width from object and extract
   ;-----------------------------------------------------------------------

   splog, 'Extracting frame '+objname+' with 3 step process'

   if (NOT keyword_set(widthset)) then begin 
     ;------------------
     ; Extract the object image

     ; Use the "whopping" terms
     ; We need to do 2 iteration extraction: 
     ;        1) Fit profiles in a subset of rows
     ;        2) Fit returned parameters with smooth functions
     ;        3) Extract all 2048 rows with new profiles given by
     ;              fitansimage

     sigma = 1.0
     proftype = 1 ; Gaussian
     wfixed = [1,1,1] ; gaussian term + centroid and  sigma terms
     nterms = 3
     sigmaterm = 1
     centerterm = 2

     ; (1) Extraction profiles in every 8th row

     splog, 'Object extraction: Step 1 (fit width)'
     extract_image, image, invvar, xnow, sigma, tempflux, tempfluxivar, $
       proftype=proftype, wfixed=wfixed, yrow=yrow, $
       highrej=highrej, lowrej=lowrej, npoly=npoly, whopping=whopping, $
       ansimage=ansimage, chisq=firstchisq

     ntrace = (size(tempflux,/dimens))[1]

     ; (2) Refit ansimage to smooth profiles


     splog, 'Answer Fitting: Step 2'

     ;---------------------------------------------------
     ;   Fitansimage is now hard wired for 320 fibers!!!!???

      fitans = fitansimage(ansimage, nterms, ntrace, npoly, yrow, $
          tempflux, fluxm=[1,1,0], scatfit=scatfit)

      fitans = fitans[0:nterms*ntrace-1,*]

      ; (3) Calculate new sigma and xtrace arrays
    
      sigmashift = transpose(fitans[lindgen(ntrace)*nterms + sigmaterm, *])
      centershift= transpose(fitans[lindgen(ntrace)*nterms + centerterm, *])

      splog, format='(a,3(f8.3))', 'Centershift ', min(centershift),  $
       median(centershift), max(centershift)

      splog, format='(a,3(f8.3))', 'Sigmashift ', min(sigmashift),  $
       median(sigmashift), max(sigmashift)

      sigma2 = sigma * (1.0 + sigmashift)

      if (max(abs(centershift)) GT maxshift OR $
       max(abs(sigmashift)) GT maxshift/3.0) then begin
         splog, 'ABORT: Shift terms are not well behaved!'
         return
      endif
     
   endif else begin

     traceset2xy, widthset, xx, sigma2
     ntrace = (size(sigma2,/dimens))[1]
     wfixed = [1,1]
     nterms = n_elements(wfixed)

     splog, 'Object extraction: Step 1 (use width from arcs)'
     extract_image, image, invvar, xnow, sigma2, tempflux, tempfluxivar, $
      proftype=proftype, wfixed=wfixed, yrow=yrow, $
      highrej=highrej, lowrej=lowrej, npoly=npoly, whopping=whopping, $
      ansimage=ansimage, chisq=firstchisq

     splog, 'Step 2: Just find scattered light image'
     junk = fitansimage(ansimage, nterms, ntrace, npoly, yrow, $
      tempflux, fluxm=[1,1], scatfit=scatfit)
   endelse

   ;-----------------------------------------------------------------------
   ;  Now, subtract scattered light and do final extraction with all rows
   ;-----------------------------------------------------------------------

   qaplot_scatlight, scatfit, yrow, $
    wset=wset, xcen=xtrace, fibermask=fibermask, $
    title=plottitle+'Scattered Light on '+objname

   ; (4) Second and final extraction
   splog, 'Object extraction: Step 3'

   highrej = 5
   lowrej = 5

   extract_image, (image - scatfit), invvar, xnow, sigma2, flux, $
    fluxivar, proftype=proftype, wfixed=wfixed, $
    highrej=highrej, lowrej=lowrej, npoly=0, whopping=whopping, $
    chisq=chisq, ymodel=ymodel2, pixelmask=pixelmask, $
    reject= [0.1,0.5,0.8]

   ;------------------
   ; QA chisq plot for fit calculated in extract image (make QAPLOT ???)

   xaxis = indgen(N_elements(chisq)) + 1
   djs_plot, xaxis, chisq, $
    xrange=[0,N_elements(chisq)], xstyle=1, $
    xtitle='Row number',  ytitle = '\chi^2', $
    title=plottitle+'Extraction chi^2 for '+objname

   djs_oplot, yrow, firstchisq[yrow], color='green', ps=4

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

   ;
   ; Sky-subtract again, this time with dispset (PSF subtraction)
   ; 

   skystruct_psf = skysubtract(flux, fluxivar, plugsort, vacset, $
    skysubpsf, skysubpsfivar, iskies=iskies, pixelmask=pixelmask, $
    fibermask=fibermask, upper=3.0, lower=3.0, dispset=dispset)

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

   fluxfactor = fluxcorr(flambda, flambdaivar, vacset, plugsort, $
    color=color, lower=3.0, upper=3.0, fibermask=fibermask)

   minfluxfactor = median(fluxfactor) * 0.01
   divideflat, flambda, flambdaivar, fluxfactor, minval=minfluxfactor

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
   sxaddpar, objhdr, 'SKIPROW', skiprow, 'Extraction: Number of rows skipped in step 1'
   sxaddpar, objhdr, 'LOWREJ', lowrej, 'Extraction: low rejection'
   sxaddpar, objhdr, 'HIGHREJ', highrej, 'Extraction: high rejection'
   sxaddpar, objhdr, 'SCATPOLY', npoly, 'Extraction: Order of scattered light polynomial'
   sxaddpar, objhdr, 'PROFTYPE', proftype, 'Extraction profile: 1=Gaussian'
   sxaddpar, objhdr, 'NFITPOLY', nterms, 'Extraction: Order of profile parameter fit'

   ;----------
   ; Write extracted, lambda-calibrated, sky-subtracted spectra to disk

   mwrfits, flambda, outname, objhdr, /create
   mwrfits, flambdaivar, outname
   mwrfits, finalmask, outname
   mwrfits, vacset, outname
   mwrfits, dispset, outname
   mwrfits, plugsort, outname

; Save sky, fluxfactor, telluricfactor???
   mwrfits, fluxfactor, outname
   if (color EQ 'red') then mwrfits, telluricfactor, outname

   return
end
;------------------------------------------------------------------------------
