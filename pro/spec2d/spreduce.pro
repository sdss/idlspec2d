;+
; NAME:
;   spreduce
;
; PURPOSE:
;   Extract, wavelength-calibrate, and flatten SDSS spectral frame(s).
;
; CALLING SEQUENCE:
;   spreduce, flatname, arcname, objname, pixflatname=pixflatname, $
;    plugfile=plugfile, lampfile=lampfile, $
;    indir=indir, plugdir=plugdir, outdir=outdir
;
; INPUTS:
;   flatname   - Name(s) of flat-field SDSS image(s)
;   arcname    - Name(s) of arc SDSS image(s)
;   objname    - Name(s) of object SDSS image(s)
;
; REQUIRED KEYWORDS:
;   plugfile   - Name of plugmap file (Yanny parameter file)
;
; OPTIONAL KEYWORDS:
;   pixflatname- Name of pixel-to-pixel flat, produced with SPFLATTEN.
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in the IDL path.
;   indir      - Input directory for FLATNAME, ARCNAME, OBJNAME;
;                default to '.'
;   plugdir    - Input directory for PLUGFILE; default to '.'
;   outdir     - Directory for output files; default to '.'
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Should test that arcs and flats are valid images with CHECKFLAVOR.
;
; PROCEDURES CALLED:
;   djs_locate_file()
;   djs_median()
;   djs_plot
;   extract_boxcar()
;   extract_image
;   fiberflat()
;   fit_skyset
;   fitarcimage
;   fluxcorr()
;   locateskylines
;   qaplot_arcline
;   qaplot_fflat
;   qaplot_scatlight
;   qaplot_skyline
;   qaplot_skysub
;   qaskylines
;   readcol
;   readfits()
;   sdssproc
;   skysubtract
;   splog
;   telluric_corr
;   tweaktrace
;   traceset2xy
;   writespectra
;   xy2traceset
;   yanny_free
;   yanny_read
;
; REVISION HISTORY:
;   12-Oct-1999  Written by D. Schlegel & S. Burles, APO
;-
;------------------------------------------------------------------------------

pro spreduce, flatname, arcname, objname, pixflatname=pixflatname, $
 plugfile=plugfile, lampfile=lampfile, $
 indir=indir, plugdir=plugdir, outdir=outdir

   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(plugdir)) then plugdir=indir
   if (NOT keyword_set(outdir)) then outdir = '.'

   stime0 = systime(1)

   ;---------------------------------------------------------------------------
   ; Locate skyline file for sky wavelength calibration
   ;---------------------------------------------------------------------------

   if (keyword_set(skylinefile)) then begin
      skyfilenames = (findfile(skylinefile, count=ct))[0]
      if (ct EQ 0) then message, 'No SKYLINEFILE found '+skylinefile
   endif else begin
      skydefault = getenv('IDLSPEC2D_DIR') + '/etc/skylines.dat'
      skyfilenames = (findfile(skydefault, count=ct))[0]
      if (skyfilenames EQ '') then message, 'No SKYLINEFILE found '+skydefault
   endelse

   skylinefile = skyfilenames[0]

   ;---------------------------------------------------------------------------
   ; Determine spectrograph ID and color from first object file
   ;---------------------------------------------------------------------------

   sdssproc, objname[0], indir=indir, spectrographid=spectrographid, color=color

   ;---------------------------------------------------------------------------
   ; Read PLUGMAP file and sort
   ;---------------------------------------------------------------------------
 
   plugpath = filepath(plugfile, root_dir=plugdir)
   plugfilename = (findfile(plugpath, count=ct))[0]
   if (ct NE 1) then $
    message, 'Cannot find plugMapFile ' + plugfile

   yanny_read, plugfilename, pstruct, hdr=hdrplug
   plugmap = *pstruct[0]
   yanny_free, pstruct

   ;-------------------------------------------------------------------------
   ; Plugsort will return mask of good (1) and bad (0) fibers too
   ;-------------------------------------------------------------------------

   plugsort = sortplugmap(plugmap, spectrographid, fibermask)
 
   ;---------------------------------------------------------------------------
   ; LOOP THROUGH FLAT+ARC IMAGE TO IDENTIFY THE BEST PAIR
   ;---------------------------------------------------------------------------

   nflat = N_elements(flatname)
   narc = N_elements(arcname)

   ibest = -1 ; Index number for best flat+arc pair
   bestcorr = -2.0 ; Any call to FITARCIMAGE returns values [-1,1]
   oldflatfile = ''

   for ifile=0, (nflat<narc)-1 do begin

      splog, ifile+1, (nflat<narc), $
       format='("Looping through flat+arc pair #",I3," of",I3)'

      ; Check if this is the same flat field as the last one read

      if (flatname[ifile] NE oldflatfile) then begin

         ;---------------------------------------------------------------------
         ; Read flat-field image
         ;---------------------------------------------------------------------

         splog, 'Reading flat ', flatname[ifile]
         sdssproc, flatname[ifile], image, invvar, indir=indir, $
          hdr=flathdr, pixflatname=pixflatname, nsatrow=nsatrow

         ;-----
         ; Decide if this flat is bad:
         ;   Reject if more than 1% of the pixels are marked as bad.
         ;   Reject if more than 10 rows are saturated.

         fbadpix = float(N_elements(where(invvar EQ 0))) / N_elements(invvar)

         qbadflat = 0
         if (fbadpix GT 0.01) then begin
            qbadflat = 1
            splog, 'Reject flat ' + flatname[ifile] + $
             ' (' + string(format='(i3)', fix(fbadpix*100)) + '% bad pixels)'
         endif
         if (nsatrow GT 10) then begin
            qbadflat = 1
            splog, 'Reject flat ' + flatname[ifile] + $
             ' (' + string(format='(i4)', nsatrow) + ' saturated rows)'
         endif

         if (NOT qbadflat) then begin
            ;------------------------------------------------------------------
            ; Create spatial tracing from flat-field image
            ;------------------------------------------------------------------

            splog, 'Tracing 320 fibers in ',  flatname[ifile]
            tmp_xsol = trace320crude(image, yset=ycen, maxdev=0.15)

            splog, 'Fitting traces in ',  flatname[ifile]
            xy2traceset, ycen, tmp_xsol, tset, ncoeff=5, maxdev=0.1
            traceset2xy, tset, ycen, tmp_xsol
         endif

         oldflatfile = flatname[ifile]
      endif

      if (NOT qbadflat) then begin

         ;---------------------------------------------------------------------
         ; Read the arc
         ;---------------------------------------------------------------------

         splog, 'Reading arc ', arcname[ifile]
         sdssproc, arcname[ifile], image, invvar, indir=indir, $
          hdr=archdr, pixflatname=pixflatname, nsatrow=nsatrow

         ;-----
         ; Decide if this arc is bad:
         ;   Reject if more than 40 rows are saturated.

         qbadarc = 0
         if (nsatrow GT 40) then begin
            qbadarc = 1
            splog, 'Reject arc ' + arcname[ifile] + $
             ' (' + string(format='(i4)', nsatrow) + ' saturated rows)'
         endif

         if (NOT qbadarc) then begin
            ;------------------------------------------------------------------
            ; Extract the arc image
            ;------------------------------------------------------------------

            splog, 'Extracting arc image with simple gaussian'
            sigma = 1.0
            proftype = 1 ; Gaussian
            highrej = 15
            lowrej = 15
            npoly = 1 ; maybe more structure
            wfixed = [1,1] ; Just fit the first gaussian term

            extract_image, image, invvar, tmp_xsol, sigma, flux, fluxivar, $
             proftype=proftype, wfixed=wfixed, $
             highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1

            ;-------------------------------------------------------------------
            ; Compute correlation coefficient for this arc image
            ;-------------------------------------------------------------------

            splog, 'Searching for wavelength solution'
            tmp_aset = 0
; FOR NOW, REVERT TO THE OLD CODE! ???
;            fitarcimage, flux, fluxivar, aset=tmp_aset, $
;             color=color, lampfile=lampfile, bestcorr=corr
            fitarcimage_old, flux, fluxivar, aset=tmp_aset, $
             color=color, lampfile=lampfile, bestcorr=corr

            ;-----
            ; Determine if this is the best flat+arc pair
            ; If so, then save the information that we need

            if (corr GT bestcorr) then begin
               ibest = ifile
               bestcorr = corr
               arcimg = flux
               arcivar = fluxivar
               xsol = tmp_xsol
               aset = tmp_aset
            endif
         endif
      endif

   endfor

   ;---------------------------------------------------------------------------
   ; Make sure that the best flat+arc pair is good enough
   ;---------------------------------------------------------------------------

   if (ibest EQ -1) then begin
      splog, 'ABORT: No good flats and/or no good arcs (saturated?)'
      return
   endif

   splog, 'Best flat = ', flatname[ibest]
   splog, 'Best arc = ', arcname[ibest]

   if ((color EQ 'blue' AND bestcorr LT 0.5) $
    OR (color EQ 'red'  AND bestcorr LT 0.5) ) then begin
      splog, 'ABORT: Best arc correlation = ', bestcorr
   endif else $
    if ((color EQ 'blue' AND bestcorr LT 0.7) $
    OR (color EQ 'red'  AND bestcorr LT 0.7) ) then begin
      splog, 'WARNING: Best arc correlation = ', bestcorr
      return
   endif

   ;---------------------------------------------------------------------------
   ; Compute wavelength calibration for arc lamp only
   ;---------------------------------------------------------------------------

   arccoeff = 5

   splog, 'Searching for wavelength solution'
; FOR NOW, REVERT TO THE OLD CODE! ???
;   fitarcimage, arcimg, arcivar, xpeak, ypeak, wset, ncoeff=arccoeff, $
;    aset=aset, $
;    color=color, lampfile=lampfile, lambda=lambda, xdif_tset=xdif_tset
   fitarcimage_old, arcimg, arcivar, xpeak, ypeak, wset, ncoeff=arccoeff, $
    aset=aset, $
    color=color, lampfile=lampfile, lambda=lambda, xdif_tset=xdif_tset

   if (NOT keyword_set(wset)) then begin
      splog, 'ABORT: Wavelength solution failed'
      return 
   endif

   qaplot_arcline, xdif_tset, lambda, filename=arcname[ibest], color=color

   ;---------------------------------------------------------------------
   ; Read best flat-field image (again)
   ;---------------------------------------------------------------------

   splog, 'Reading flat ', flatname[ibest]
   sdssproc, flatname[ibest], image, invvar, indir=indir, $
    hdr=flathdr, pixflatname=pixflatname

   ;---------------------------------------------------------------------------
   ; Extract the flat-field image
   ;---------------------------------------------------------------------------

   splog, 'Extracting flat-field image with simple gaussian'
   sigma = 1.0
   proftype = 1 ; Gaussian
   highrej = 15
   lowrej = 15
   npoly = 1  ; just fit flat background to each row
   wfixed = [1,1] ; Just fit the first gaussian term

   extract_image, image, invvar, xsol, sigma, flat_flux, flat_fluxivar, $
    proftype=proftype, wfixed=wfixed, $
    highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1

   junk = where(flat_flux GT 1.0e5, nbright)
   splog, 'Found ', nbright, ' bright pixels in extracted flat ', $
    flatname[ibest], format='(a,i7,a,a)'

   ;---------------------------------------------------------------------------
   ; Compute fiber-to-fiber flat-field variations
   ;---------------------------------------------------------------------------

;   fflat = fiberflat(flat_flux, flat_fluxivar, wset, fibermask=fibermask)
   fflat = fiberflat(flat_flux, flat_fluxivar, wset, fibermask=fibermask, $
    /dospline)

   qaplot_fflat, fflat, wset, filename=flatname[ibest]

   ;---------------------------------------------------------------------------
   ; LOOP THROUGH OBJECT FRAMES
   ;---------------------------------------------------------------------------

   for iobj=0, N_elements(objname)-1 do begin

      splog, camname=objname[iobj]

      ;------------------
      ; Read object image

      splog, 'Reading object ', objname[iobj]
      sdssproc, objname[iobj], image, invvar, indir=indir, hdr=objhdr, $
       pixflatname=pixflatname, spectrographid=spectrographid, color=color

      ;------------------
      ; Tweak up the spatial traces
      ; Currently doing this with 6 step extraction

      ;------------------
      ; Identify very bright objects
      ; Do a boxcar extraction, and look for fibers where the median
      ; counts are 10000 ADU per row.

      fextract = extract_boxcar(image, xsol)

      scrunch = djs_median(fextract, 1) ; Find median counts/row in each fiber
      fullscrunch = djs_median(fextract) ; Find median counts/row in all fibers
      whopping = where(scrunch - fullscrunch GT 10000.0, whopct)
      splog, 'Median counts in all fibers = ', fullscrunch
      splog, 'Number of bright fibers = ', whopct

      if (whopct GT 20) then begin
         splog, 'WARNING: Disable whopping terms ' + objname[iobj]
         whopping = -1
         whopct = 0
      endif

      ;------------------
      ; Extract the object image
      ; Use the "whopping" terms
      ; We need to do 2 iteration extraction: 
      ;        1) Fit profiles in a subset of rows
      ;        2) Fit returned parameters with smooth functions
      ;        3) Extract all 2048 rows with new profiles given by
      ;              fitansimage

      splog, 'Extracting frame '+objname[iobj]+' with 6 step process'
      nrow = (size(image))[2]
      ncol = (size(image))[1]
      skiprow = 8
      yrow = lindgen(nrow/skiprow) * skiprow + skiprow/2
      nfirst = n_elements(yrow)
  
      sigma = 1.0
      proftype = 1 ; Gaussian
      highrej = 5  ; just for first extraction steps
      lowrej = 5  ; just for first extraction steps
      npoly = 16 ; maybe more structure, lots of structure
      wfixed = [1,1,1] ; gaussian term + centroid and  sigma terms
      nterms = 3
      sigmaterm = 1
      centerterm = 2
      xnow = xsol ; Is this modified when extracting???
      sigmanow = xsol*0.0 + sigma

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

         ; (3) Calculate new sigma and xsol arrays
      
         if (i EQ 0) then begin 
            splog, 'Trace Tweaking: Step', i*3+3, '    (Sigma not tweaked)'
            sigmashift = $
             transpose(fitans[lindgen(ntrace)*nterms + sigmaterm, *])
            centershift = $
             transpose(fitans[lindgen(ntrace)*nterms + centerterm, *])
            tweaktrace, xnow, sigmanow, centershift, sigmashift
         endif

      endfor

      qaplot_scatlight, scatimage, scatfit, wset=wset, xcen=xsol, $
       fibermask=fibermask, filename=objname[iobj]

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
       chisq=chisq, ymodel=ymodel2

      ;------
      ; QA chisq plot for fit calculcated in extract image (make QAPLOT ???)

      xaxis = indgen(N_elements(chisq)) + 1
      djs_plot, xaxis, chisq, $
       xtitle='Row number',  ytitle = '\chi^2', $
       title='Extraction chi^2 for '+objname[iobj]

      ;------------------
      ; Flat-field the extracted object fibers with the global flat
      divideflat, flux, fluxivar, fflat, fibermask=fibermask

      ;------------------
      ; Tweak up the wavelength solution to agree with the sky lines.
      ; xshift contains polynomial coefficients to shift arc to sky lines.
  
      locateskylines, skylinefile, flux, fluxivar, $
       wset, xsky, ysky, skywaves, xset=xset

      ;------------------
      ; First convert lambda, and skywaves to log10 vacuum

      splog, 'Converting wavelengths to vacuum'
      vaclambda = lambda
      airtovac, vaclambda
      vacloglam = alog10(vaclambda)

      vacsky = skywaves
      airtovac, vacsky
      vaclogsky = alog10(vacsky)

      sxaddpar, hdr, 'VACUUM', 'WAVELENGTHS ARE IN VACUUM'
      sxaddpar, hdr, 'AIR2VAC', systime()

      splog, 'Tweaking to sky lines'

      ;------------------
      ; Fit to arc lines with sky shifts included

      ; This procedure produces numerical steps at 10^-6, 1 km/s

      xshift = xpeak * 0.0

      if (size(xset,/tname) NE 'UNDEFINED') then begin
         traceset2xy, xset, transpose(xpeak), xshift
         xshift = transpose(xshift)
      endif else begin
         splog, 'WARNING: Sky lines are too noisy! No shifting!'
      endelse

      ; Move this QA plot elsewhere ???
      plot, xpeak, xshift, ps=3, xtitle='Arc line position', /ynozero, $
       ytitle = 'Offset to Sky Line [pix]', $
       title = 'Offset to Sky Lines From Wavelength-Solution'

      vacset = wset
      fixabove = 2
      finalarcfit, xpeak+xshift, vacloglam, vacset, arccoeff, fixabove, $
       fibermask=fibermask, maxdev=0, maxiter=1, nsetcoeff=8, maxsig=2.0

;      fit_skyset, xpeak, ypeak, vacloglam, xsky, ysky, vaclogsky, skycoeff, $
;        goodlines, wset, ymin=ymin, ymax=ymax, func=func

      ;------------------
      ; Sky-subtract

      skysubtract, flux, fluxivar, plugsort, vacset, $
       skysub, skysubivar, fibermask=fibermask, upper=3.0, lower=3.0

      qaplot_skysub, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, fibermask=fibermask, filename=objname[iobj]

      ; QA for 2 skylines in the blue
      qaplot_skyline, 4359.5, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, fibermask=fibermask, dwave=3.0, $
       filename=objname[iobj]
      qaplot_skyline, 5578.9, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, fibermask=fibermask, dwave=5.0, $
       filename=objname[iobj]

      ; QA for 2 skylines in the red
      qaplot_skyline, 7343.0, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, fibermask=fibermask, dwave=5.0, $
       filename=objname[iobj]
      qaplot_skyline, 8888.3, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, fibermask=fibermask, dwave=5.0, $
       filename=objname[iobj]

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
  
      if (color EQ 'red')  then begin

         ;-----------------------------------------------
         ;  Split into two regions, A,B bands first
         telluric1 = telluric_corr(flambda, flambdaivar, vacset, plugsort, $
          minw=3.82, maxw=3.92, lower=5.0, upper=5.0, ncontbkpts=10)
  
         ;-----------------------------------------------
         ;  9100 Ang absorption next?
         telluric2 = telluric_corr(flambda, flambdaivar, vacset, plugsort, $
          minw=3.94, maxw=3.97, lower=5.0, upper=5.0)

         telluricfactor = telluric1 * telluric2
         divideflat, flambda, flambdaivar, telluricfactor, minval=0.1
  
         qaskylines, flambda, flambdaivar, vacset, plugsort
      endif

      ;------------------
      ; Write extracted, lambda-calibrated, sky-subtracted spectra to disk

      framenum = sxpar(objhdr, 'EXPOSURE')
  
      filebase = filepath( $
       's-'+string(format='(i1,a1,a,i4.4)',spectrographid, $
       color,'-',framenum), root_dir=outdir)
  
      ;------
      ; Add keywords to object header

      sxaddpar, objhdr, 'PLUGMAPF', plugfilename
      sxaddpar, objhdr, 'FLATFILE', flatname[ibest]
      sxaddpar, objhdr, 'ARCFILE', arcname[ibest]
      sxaddpar, objhdr, 'OBJFILE', objname[iobj]
      sxaddpar, objhdr, 'LAMPLIST', lampfile
      sxaddpar, objhdr, 'SKYLIST', skylinefile
      sxaddpar, objhdr, 'PIXFLAT', pixflatname
      sxaddpar, objhdr, 'OSIGMA',  sigma, $
       'Original guess at spatial sigma in pix'
      sxaddpar, objhdr, 'SKIPROW', skiprow, 'Number of rows skipped in step 1'
      sxaddpar, objhdr, 'LOWREJ', lowrej, 'Extraction, low rejection'
      sxaddpar, objhdr, 'HIGHREJ', highrej, 'Extraction, high rejection'
      sxaddpar, objhdr, 'SCATPOLY', npoly, 'Order of scattered light poly'
      sxaddpar, objhdr, 'PROFTYPE', proftype, '1=Gaussian'
      sxaddpar, objhdr, 'NFITPOLY', nparams, 'Order of profile parameter fit'

      writespectra, objhdr, flambda, flambdaivar, plugsort, vacset, $
       filebase=filebase

      heap_gc   ; Garbage collection for all lost pointers

   endfor

   splog, 'Elapsed time = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)', camname=''

   return
end
;------------------------------------------------------------------------------
