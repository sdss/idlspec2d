;+
; NAME:
;   extract_object
;
; PURPOSE:
;
;   Performs all object extraction tasks
;      0) Locate bright fibers, and test image background
;      1) 6 step Optimal extraction
;      2) Tweak to sky lines
;      3) Sky subtraction
;      4) Flux calibration
;      5) Telluric correction
;
; CALLING SEQUENCE:
;   extract_object, outname, objhdr, image, invvar, plugsort, wset, $
;                    xarc, lambda, xtrace, fflat, fibermask, color=color
;
; INPUTS:
;   outname    - Name of outputs fits file
;   objhdr     - Header cards from object image
;   image      - Object image [nx,ny]
;   invvar     - Inverse Variance of object image [nx,ny]
;   plugsort   - Plugmap structure for [ntrace] spectra
;   wset       - wavelength solution from arc line spectra
;   xarc       - centroids measured in arc line spectra
;   lambda     - air wavelengths corresponding to xarc
;   xtrace     - spatial traces from flat field
;   fflat      - 1d flat field vectors
;   fibermask  - information bits on fiber status
;
; REQUIRED KEYWORDS:
;   color      - camera color (red or blue)
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   A fits file is output in outname, which contains
;      float flux [nx,ntrace]
;      float flux_invvar [nx,ntrace]
;      plugsort struct [ntrace]
;      vacuum wavelength coefficients 
;      integer pixelmask [nx,ntrace]
;      byte    fibermask [ntrace]
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
;   djs_median()
;   djs_plot
;   extract_boxcar()
;   extract_image
;   fibermask_bits()
;   fitansimage()
;   fitvacset()
;   fluxcorr()
;   locateskylines
;   mwrfits
;   pixelmask_bits()
;   qaplot_scatlight
;   qaplot_skyline
;   qaplot_skysub
;   qaplot_skydev
;   skysubtract
;   splog
;   telluric_corr
;   tweaktrace
;
; REVISION HISTORY:
;   24-Jan-2000  Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
pro extract_object, outname, objhdr, image, invvar, plugsort, wset, $
               xarc, lambda, xtrace, fflat, fibermask, color=color

      skylinefile = strtrim(sxpar(objhdr,'SKYLIST'),2) 
      objname = strtrim(sxpar(objhdr,'OBJFILE'),2) 
 
      ;------------------
      ; Identify very bright objects
      ; Do a boxcar extraction, and look for fibers where the median
      ; counts are 10000 ADU per row.

      fextract = extract_boxcar(image, xtrace)

      scrunch = djs_median(fextract, 1) ; Find median counts/row in each fiber
      scrunch_sort = sort(scrunch)
      nfiber = (size(fextract))[2]
      i5 = nfiber/20
      i95 = i5 * 19

      fullscrunch = djs_median(fextract) ; Find median counts/row in all fibers
      whopping = where(scrunch - fullscrunch GT 10000.0, whopct)
      splog, 'Median counts in all fibers = ', fullscrunch
      splog, 'Number of bright fibers = ', whopct

      iskies = where(plugsort.objtype EQ 'SKY' AND plugsort.fiberid GT 0 AND $
                     (fibermask EQ 0), nskies)

      if (nskies LE 1) then begin
              splog, 'ABORT: Only '+ string(nskies) + ' sky fibers found' 
              return
      endif 

      skymedian = djs_median(scrunch[iskies])
      splog, 'Sky fiber median '+string(skymedian)
      if (skymedian GT 3000.0) then begin
        splog, 'ABORT: Sky fibers are brighter than 3000 counts'
        return
      endif
          
      splog, '5% and 95% count levels ', scrunch[scrunch_sort[i5]], $
                                         scrunch[scrunch_sort[i95]]

;      if (scrunch[scrunch_sort[i5]] GT 5000.0) then begin
;         splog, 'ABORT: Fibers have '+ string(scrunch[scrunch_sort[i5]]) + $
;                ' at the 5% '
         ;
         ;	Here we just need to move to the next exposure, not return
         ;       completely.  This would argue for having a subroutine called
         ;       for each exposure, so a return would only skip the botched
         ;       exposure, and not all the rest
         ;
;	 return
;      endif


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
    

      ;------------------
      ; Extract the object image
      ; Use the "whopping" terms
      ; We need to do 2 iteration extraction: 
      ;        1) Fit profiles in a subset of rows
      ;        2) Fit returned parameters with smooth functions
      ;        3) Extract all 2048 rows with new profiles given by
      ;              fitansimage

      splog, 'Extracting frame '+objname+' with 6 step process'
      nrow = (size(image))[2]
      ncol = (size(image))[1]
      skiprow = 8
      yrow = lindgen(nrow/skiprow) * skiprow + skiprow/2
      nfirst = n_elements(yrow)
  
      sigma = 1.0
      proftype = 1 ; Gaussian
      highrej = 5  ; just for first extraction steps
      lowrej = 5  ; just for first extraction steps
      npoly = 8 ; maybe more structure, lots of structure
      wfixed = [1,1,1] ; gaussian term + centroid and  sigma terms
      nterms = 3
      sigmaterm = 1
      centerterm = 2
      xnow = xtrace; Is this modified when extracting???
      sigmanow = xtrace*0.0 + sigma
      maxshift = 1.5

      ; Kill or adjust first and last column ???
      ;invvar[0,*] = 0.0
      ;invvar[2047,*] = 0.0
      image[0,*] = image[0,*]*0.7
      image[2047,*] = image[2047,*]*0.7

      ; My fitansimage is not going to work well without good profiles ???
      ; Using it to tweak up traces, and then performing free fit to object

      for i = 0, 1 do begin
      
         ; (1) First extraction

         splog, 'Object extraction: Step', i*3+1
         extract_image, image, invvar, xnow, sigma, tempflux, tempfluxivar, $
          proftype=proftype, wfixed=wfixed, yrow=yrow, $
          highrej=highrej, lowrej=lowrej, npoly=npoly, whopping=whopping, $
          ansimage=ansimage

         ; (2) Refit ansimage to smooth profiles

         splog, 'Answer Fitting: Step', i*3+2
         nparams = 3
         ntrace = (size(tempflux))[2]
         fitans = fitansimage(ansimage, nparams, ntrace, npoly, nfirst, yrow, $
          fluxm = [1,1,0], crossfit=1, scatfit=scatfit, scatimage=scatimage)

         ; (3) Calculate new sigma and xtrace arrays
      
         sigmashift = transpose(fitans[lindgen(ntrace)*nterms + sigmaterm, *])
         centershift= transpose(fitans[lindgen(ntrace)*nterms + centerterm, *])

         splog, format='(a,3(f8.3))', 'Centershift ', min(centershift),  $
          median(centershift), max(centershift)

         splog, format='(a,3(f8.3))', 'Sigmashift ', min(sigmashift),  $
          median(sigmashift), max(sigmashift)

         if (max(abs(centershift)) GT maxshift OR $
             max(abs(sigmashift)) GT maxshift/3.0) then begin
              splog, 'ABORT: Shift terms are not well behaved!'
              return
         endif

         if (i EQ 0) then begin 
            splog, 'Trace Tweaking: Step', i*3+3, '    (Sigma not tweaked)'
	    ;----------------------------------------------------
            ;  This is essentially tweaktrace
            ;
            xnow = xnow - centershift * sigma
            sigmanow = (sigmashift + 1.0) * sigmanow

            ;tweaktrace, xnow, sigmanow, centershift, sigmashift
         endif
      endfor

      qaplot_scatlight, scatfit, yrow, $
       wset=wset, xcen=xtrace, fibermask=fibermask, filename=objname

      ; (4) Second and final extraction
      splog, 'Object extraction: Step 6'

      ; Using old sigma for now, which should be fine
      ; Different sigmas require a new profile for each trace, so will
      ; check timing in the future
      ; subtract off fit to scattered light and don't allow polynomial terms

;      extract_image, image-scatfit, invvar, xnow, sigma, flux, $
;       fluxivar, proftype=proftype, wfixed=wfixed, $
;       highrej=highrej, lowrej=lowrej, npoly=0, whopping=whopping, chisq=chisq

      highrej = 15
      lowrej = 15
      fitans = fitans[0:nparams*ntrace-1,*]
      extract_image, image-scatfit, invvar, xnow, sigma, flux, $
       fluxivar, proftype=proftype, wfixed=wfixed, fitans=fitans, $
       highrej=highrej, lowrej=lowrej, npoly=0, whopping=whopping, $
       chisq=chisq, ymodel=ymodel2, pixelmask=pixelmask
 

      ;------
      ; QA chisq plot for fit calculated in extract image (make QAPLOT ???)

      xaxis = indgen(N_elements(chisq)) + 1
      djs_plot, xaxis, chisq, $
         xtitle='Row number',  ytitle = '\chi^2', $
         title='Extraction chi^2 for '+objname

      ;------------------
      ; Flat-field the extracted object fibers with the global flat
      divideflat, flux, fluxivar, fflat, fibermask=fibermask

      lowflat = where(fflat LT 0.5)
      if (lowflat[0] NE -1) then pixelmask[lowflat] = $
                   pixelmask[lowflat] OR pixelmask_bits('LOWFLAT')

      ;--------------------------------------------------------------
      ;  Look for pixels where scattered light is dominating

      scatteredpix = where(extract_boxcar(scatfit, xnow, radius=2.0) GT flux)
      if (scatteredpix[0] NE -1) then pixelmask[scatteredpix] = $
                    pixelmask[scatteredpix] + pixelmask_bits('SCATTEREDLIGHT')

      ;------------------
      ; Tweak up the wavelength solution to agree with the sky lines.
      ; xshet contains polynomial coefficients to shift arc to sky line frame.
  
      locateskylines, skylinefile, flux, fluxivar, $
         wset, xsky, ysky, skywaves, xset=xset

      if (size(xset,/tname) EQ 'UNDEFINED') then begin
         splog, 'ABORT: Cannot shift to sky lines...'
         return
      endif 

      ; ----------------------------------------
      ;  fitvacset performs shift to skylines and fit to vacuum wavelengths
      
      vacset = fitvacset(xarc, lambda, wset, xset, ncoeff=arccoeff)

      sxaddpar, objhdr, 'VACUUM', 'WAVELENGTHS ARE IN VACUUM'
      sxaddpar, objhdr, 'AIR2VAC', systime()


      ;------------------
      ; Sky-subtract

      skystruct = skysubtract(flux, fluxivar, plugsort, vacset, $
        skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, $
        fibermask=fibermask, upper=3.0, lower=3.0, $
        relchi2struct=relchi2struct)

      ; plot, skystruct.wave, skystruct.flux, ps=3

      qaplot_skysub, flux, fluxivar, skysub, skysubivar, $
        vacset, iskies, filename=objname

      qaplot_skydev, flux, fluxivar, vacset, plugsort, color, $
              filename = objname

      ; QA for 2 skylines in the blue
      if (color EQ 'blue') then begin
        qaplot_skyline, 4359.5, flux, fluxivar, skysub, skysubivar, $
         plugsort, vacset, iskies, fibermask=fibermask, dwave=4.0, $
         filename=objname
        qaplot_skyline, 5578.9, flux, fluxivar, skysub, skysubivar, $
         plugsort, vacset, iskies, fibermask=fibermask, dwave=5.0, $
         filename=objname
      endif

      ; QA for 2 skylines in the red
      if (color EQ 'red') then begin
        qaplot_skyline, 7343.0, flux, fluxivar, skysub, skysubivar, $
         plugsort, vacset, iskies, fibermask=fibermask, dwave=7.0, $
         filename=objname
        qaplot_skyline, 8888.3, flux, fluxivar, skysub, skysubivar, $
         plugsort, vacset, iskies, fibermask=fibermask, dwave=7.0, $
         filename=objname
      endif

      ;------------------------------------------
      ; Flux calibrate to spectrophoto_std fibers

      fluxfactor = fluxcorr(skysub, skysubivar, vacset, plugsort, $
       color=color, lower=3.0, upper=3.0, fibermask=fibermask)

      flambda  = skysub
      flambdaivar = skysubivar

      minfluxfactor = median(fluxfactor)*0.01
      divideflat, flambda, flambdaivar, fluxfactor, minval=minfluxfactor

      ;------------------------------------------
      ; Telluric correction called for 'red' side
      ;
      ;  May want to move all of the telluric_corr and plotting into
      ;  new procedure: telluric_fit,flambda, flambdaivar, vacset, plugsort 

 
      if (color EQ 'red')  then begin

         ;-----------------------------------------------
         ;  Split into two regions, A,B bands first
         telluric1 = telluric_corr(flambda, flambdaivar, vacset, plugsort, $
            contwave=contwave1, contflux=contflux1, contivar=contivar1, $
            telluricbkpt=telluricbkpt1, telluriccoeff=telluriccoeff1, $
            minw=3.82, maxw=3.92, lower=5.0, upper=5.0, ncontbkpts=10, $
            fibermask=fibermask)
  
         ;-----------------------------------------------
         ;  9100 Ang absorption next?
         telluric2 = telluric_corr(flambda, flambdaivar, vacset, plugsort, $
            contwave=contwave2, contflux=contflux2, contivar=contivar2,    $
            telluricbkpt=telluricbkpt2, telluriccoeff=telluriccoeff2, $
            minw=3.94, maxw=3.97, lower=5.0, upper=5.0, ncontbkpts=5, $
            fibermask=fibermask)

         psave = !p.multi
 	 !p.multi = [0,1,3]
         djs_plot,10^contwave1,contflux1,ps=3,xr=10^[3.82,3.87],yr=[0.0,1.5], $
              ymargin=[2,4], charsize=1.6, xstyle=1, $ 
              xtitle='\lambda [A]', ytitle='Flux [electrons]', $
              title = 'Telluric correction for '+objname
         djs_oplot,10^contwave1,slatec_bvalu(contwave1,telluricbkpt1, $
                      telluriccoeff1),color='red'

         djs_plot,10^contwave1,contflux1,ps=3,xr=10^[3.87,3.92],yr=[0.0,1.5], $
              ymargin=[2,2], charsize=1.6, xstyle=1, $ 
              xtitle='\lambda [A]', ytitle='Flux [electrons]'
         djs_oplot,10^contwave1,slatec_bvalu(contwave1,telluricbkpt1, $
                      telluriccoeff1),color='red'

         djs_plot,10^contwave2,contflux2,ps=3,yr=[0.0,1.5], $
              ymargin=[4,2], charsize=1.6, xstyle=1, $ 
              xtitle='\lambda [A]', ytitle='Flux [electrons]'
         djs_oplot,10^contwave2,slatec_bvalu(contwave2,telluricbkpt2, $
                      telluriccoeff2),color='red'

 	 !p.multi = psave

         telluricfactor = telluric1 * telluric2
         divideflat, flambda, flambdaivar, telluricfactor, minval=0.1
  
      endif

      ;------
      ; Add keywords to object header

      sxaddpar, objhdr, 'OSIGMA',  sigma, $
       'Original guess at spatial sigma in pix'
      sxaddpar, objhdr, 'SKIPROW', skiprow, 'Number of rows skipped in step 1'
      sxaddpar, objhdr, 'LOWREJ', lowrej, 'Extraction, low rejection'
      sxaddpar, objhdr, 'HIGHREJ', highrej, 'Extraction, high rejection'
      sxaddpar, objhdr, 'SCATPOLY', npoly, 'Order of scattered light poly'
      sxaddpar, objhdr, 'PROFTYPE', proftype, '1=Gaussian'
      sxaddpar, objhdr, 'NFITPOLY', nparams, 'Order of profile parameter fit'

      ;------------------
      ; Write extracted, lambda-calibrated, sky-subtracted spectra to disk
  
      mwrfits, flambda, outname, objhdr, /create
      mwrfits, flambdaivar, outname
      mwrfits, plugsort, outname
      mwrfits, vacset, outname
      mwrfits, pixelmask, outname
      mwrfits, fibermask, outname

   return
end
;------------------------------------------------------------------------------
