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
;    indir=indir, plugdir=plugdir, outdir=outdir, qadir=qadir, qa=qa
;
; INPUTS:
;   flatname   - Name of flat-field SDSS image
;   arcname    - Name of arc SDSS image
;   objname    - Name of object SDSS image(s)
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
;   qadir      - Directory for QA files; default to '.'
;   qa         - QA (quality assurance flag) 
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
;   Tweaking to sky lines + vacuum wavelengths commented out 15-nov-99 (DJS).
;
; PROCEDURES CALLED:
;   djs_locate_file()
;   djs_median()
;   extract_boxcar()
;   extract_image
;   fit_skyset
;   fitarcimage
;   fluxcorr()
;   locateskylines
;   qaplot_arcline
;   readcol
;   readfits()
;   sdssproc
;   skysubtract
;   splog
;   telluric_corr
;   traceset2xy
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
 indir=indir, plugdir=plugdir, outdir=outdir, qadir=qadir, qa=qa

   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(plugdir)) then plugdir=indir
   if (NOT keyword_set(outdir)) then outdir = '.'
   if (NOT keyword_set(qadir)) then qadir = outdir

   t_begin = systime(1)

   ;---------------------------------------------------------------------------
   ; Locate skyline file for sky wavelength calibration
   ;---------------------------------------------------------------------------

   if (keyword_set(skylinefile)) then begin
      skyfilenames = findfile(skylinefile, count=ct)
      if (ct EQ 0) then message, 'No SKYLINEFILE found '+skylinefile
   endif else begin
      skydefault = 'skylines.dat'
      skyfilenames = djs_locate_file(skydefault)
      if (skyfilenames EQ '') then message, 'No SKYLINEFILE found '+skydefault
   endelse

   skylinefile = skyfilenames[0]

   ;---------------------------------------------------------------------------
   ; Read PLUGMAP file
   ;---------------------------------------------------------------------------
 
   plugpath = filepath(plugfile, root_dir=plugdir)
   plugfilenames = findfile(plugpath, count=ct)
   if (ct NE 1) then $
    message, 'Cannot find plugMapFile ' + plugfile

   yanny_read, plugfilenames[0], pstruct, hdr=hdrplug
   plugmap = *pstruct[0]
   yanny_free, pstruct

   ;---------------------------------------------------------------------------
   ; Read flat-field image
   ;---------------------------------------------------------------------------

   flatpath = filepath(flatname, root_dir=indir)
   flatfilenames = findfile(flatpath, count=ct)
   if (ct NE 1) then $
    message, 'Cannot find flat image ' + flatname

   splog, 'Reading in flat ', flatfilenames[0]
   sdssproc, flatfilenames[0], image, invvar, hdr=flathdr, $
    pixflatname=pixflatname, spectrographid=spectrographid, color=color

   ;-------------------------------------------------------------------------
   ; Plugsort will return mask of good (1) and bad (0) fibers too
   ;-------------------------------------------------------------------------
   plugsort = sortplugmap(plugmap, spectrographid, fibermask)
 
   ;---------------------------------------------------------------------------
   ; Create spatial tracing from flat-field image
   ;---------------------------------------------------------------------------

   xsol = trace320crude(image, yset=ycen, maxdev=0.15)

   xy2traceset, ycen, xsol, tset, ncoeff=5, maxdev=0.1
   traceset2xy, tset, ycen, xsol

   ;---------------------------------------------------------------------------
   ; Extract the flat-field image
   ;---------------------------------------------------------------------------

   splog, 'Extracting flat-field with simple gaussian'
   sigma = 1.0
   proftype = 1 ; Gaussian
   highrej = 20
   lowrej = 25
   nPoly = 6
   wfixed = [1] ; Just fit the first gaussian term

   extract_image, image, invvar, xsol, sigma, flux, fluxivar, $
    proftype=proftype, wfixed=wfixed, $
    highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1

   highpixels = where(flux GT 1.0e5, numhighpixels)

   splog, 'Found ', numhighpixels, ' highpixels in extracted flat ', $
    flatname

   ;---------------------------------------------------------------------------
   ; Compute fiber-to-fiber flat-field variations
   ;---------------------------------------------------------------------------

   fflat = fiberflat(flux, fluxivar, fibermask)

   ;---------------------------------------------------------------------------
   ; Read the arc
   ;---------------------------------------------------------------------------

   arcpath = filepath(arcname, root_dir=indir)
   arcfilenames = findfile(arcpath, count=ct)
   if (ct NE 1) then $
    message, 'Cannot find arc image ' + arcname

   splog, 'Reading in arc ', arcfilenames[0]
   sdssproc, arcfilenames[0], image, invvar, hdr=archdr, $
    pixflatname=pixflatname, spectrographid=spectrographid, color=color

     
  ;--------------------------------------------------------------------------
  ; Extract the arc image
  ;--------------------------------------------------------------------------

  splog, 'Extracting arc-lamp with simple gaussian'
  sigma = 1.0
  proftype = 1 ; Gaussian
  highrej = 10
  lowrej = 15
  nPoly = 6 ; maybe more structure
  wfixed = [1] ; Just fit the first gaussian term

  extract_image, image, invvar, xsol, sigma, flux, fluxivar, $
   proftype=proftype, wfixed=wfixed, $
   highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1

  ;------------------
  ; Flat-field the extracted arcs with the global flat
  ; Hmmm.... would be circular if we need a wavelength calibration before
  ; making that flat.  We don't at the moment.

  divideflat, flux, fluxivar, fflat, fibermask

;  flux = flux / fflat
;  fluxivar = fluxivar * fflat^2

  ;-------------------------------------------------------------------------
  ; Compute wavelength calibration for arc lamp only
  ;-------------------------------------------------------------------------

   splog, 'Searching for wavelength solution with fitarcimage'
   fitarcimage, flux, fluxivar, xpeak, ypeak, wset, $
    color=color, lampfile=lampfile, lambda=lambda, $
    xdif_lfit=xdif_lfit, xdif_tset=xdif_tset

   qaplot_arcline, xdif_tset, lambda, arcname

; Plot flat-field ???
plot,fflat[*,0], yr=[0,2], /ystyle, $
 xtitle='Wavelength',ytitle='Intensity', $
 title='Flat Vectors (every 20th fib)', $
 charsize=2, charthick=2
for i=0,16 do oplot,fflat[*,i*19]

   ;---------------------------------------------------------------------------
   ; LOOP THROUGH OBJECT FRAMES
   ;---------------------------------------------------------------------------

   for iobj=0, N_elements(objname)-1 do begin

      splog, 'Start time ',systime(1)-t_begin, ' seconds so far', $
       format='(A,F8.2,A)'

      ;------------------
      ; Read object image

      objpath = filepath(objname[iobj], root_dir=indir)
      objfilenames = findfile(objpath, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find object image ' + objname[iobj]

      objfile = objfilenames[0]
      splog, 'Reading in object ', objfile
      sdssproc, objfile, image, invvar, hdr=objhdr, $
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
      whopping = where(scrunch GT 10000.0, whopct)
      splog, 'Number of bright fibers = ', whopct

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
      highrej = 10
      lowrej = 15
      nPoly = 6 ; maybe more structure
      wfixed = [1,1,1] ; gaussian term + centroid and  sigma terms
      nTerms = 3
      sigmaterm = 1
      centerterm = 2
      xnow = xsol
      sigmanow = xsol*0.0 + sigma

      for i = 0, 1 do begin
      
        ; 1) First extraction
        splog, 'Object extraction: Step', i*3+1
        extract_image, image, invvar, xnow, sigma, tempflux, tempfluxivar, $
         proftype=proftype, wfixed=wfixed, yrow=yrow, $
         highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping, $
         ansimage=ansimage

        ; 2) Refit ansimage to smooth profiles

        splog, 'Answer Fitting: Step', i*3+2
        nparams = 3
        nTrace = (size(flux))[2]
        fitans = fitansimage(ansimage, nparams, nTrace, nPoly, nfirst, yrow, $
         fluxm = [1,1,0], crossfit=1-i)

        ; 3) Calculate new sigma and xsol arrays
      
        if (i EQ 0) then begin 
           splog, 'Trace Tweaking: Step', i*3+3, '    (Sigma not tweaked)'
           sigmashift = transpose(fitans[lindgen(nTrace)*nTerms + sigmaterm, *])
           centershift = $
            transpose(fitans[lindgen(nTrace)*nTerms + centerterm, *])
           tweaktrace, xnow, sigmanow, centershift, sigmashift
        endif
      endfor

      ; 4) Second and final extraction
      splog, 'Object extraction: Step 6'

      ; Using old sigma for now, which should be fine
      ; Different sigmas require a new profile for each trace, so will
      ; check timing in the future

      extract_image, image, invvar, xnow, sigma, flux, $
       fluxivar, proftype=proftype, wfixed=wfixed, fitans=fitans, $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping 
;       ymodel=ymodel2

      ;------------------
      ; Flat-field the extracted object fibers with the global flat
  divideflat, flux, fluxivar, fflat, fibermask
;      flux = flux / fflat
;      fluxivar = fluxivar * fflat^2

      ;------------------
      ; Tweak up the wavelength solution to agree with the sky lines.

      locateskylines, skylinefile, flux, fluxivar, $
       wset, xsky, ysky, skywaves, lambda=skylambda

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
      skycoeff = 2
      if (n_elements(vaclogsky) GT 3) then skycoeff = 3

      fit_skyset, xpeak, ypeak, vacloglam, xsky, ysky, vaclogsky, skycoeff, $
        goodlines, wset, ymin=ymin, ymax=ymax, func=func

      locateskylines, skylinefile, flux, fluxivar, $
       wset, xsky, ysky, skywaves, lambda=vaclogsky

      ;------------------
      ; Sky-subtract

      skysubtract, flux, fluxivar, plugsort, wset, $ 
       skysub, skysubivar

      ;------------------------------------------
      ; Flux calibrate to spectrophoto_std fibers

      fluxfactor = fluxcorr(skysub, skysubivar, wset, plugsort, $
                             lower=1.5, upper=5)

      flux = skysub * fluxfactor
      fluxivar = skysubivar / (fluxfactor^2)

      ;------------------------------------------
      ; Telluric correction called for 'red' side

      if (color EQ 'red')  then begin

         telluricfactor = telluric_corr(flux, fluxivar, wset, plugsort)
         flux = flux / telluricfactor
         fluxivar = fluxivar * (telluricfactor^2)

      endif

      ;------------------
      ; Write extracted, lambda-calibrated, sky-subtracted spectra to disk

      framenum = sxpar(objhdr, 'EXPOSURE')

      filebase = filepath( $
       's-'+string(format='(i1,a1,a,i4.4)',spectrographid, $
       color,'-',framenum), root_dir=outdir)

      ;------
      ; Add everything we can think of to object header

      sxaddpar, objhdr, 'PLUGMAPF', plugfilenames[0]
      sxaddpar, objhdr, 'FLATFILE', flatfilenames[0]
      sxaddpar, objhdr, 'ARCFILE',  arcfilenames[0]
      sxaddpar, objhdr, 'OBJFILE',  objfile
      sxaddpar, objhdr, 'LAMPLIST',  lampfile
      sxaddpar, objhdr, 'SKYLIST',  skylinefile
      sxaddpar, objhdr, 'PIXFLAT',  pixflatname
      sxaddpar, objhdr, 'OSIGMA',  sigma, $
           'Original guess at sigma of spatial profiles'
      sxaddpar, objhdr, 'SKIPROW', skiprow, 'Number of rows skipped in step 1'
      sxaddpar, objhdr, 'LOWREJ', lowrej, 'Extraction, low rejection'
      sxaddpar, objhdr, 'HIGHREJ', highrej, 'Extraction, high rejection'
      sxaddpar, objhdr, 'SCATPOLY', nPoly, 'Order of scattered light poly'
      sxaddpar, objhdr, 'PROFTYPE', proftype, '1 is Gaussian'
      sxaddpar, objhdr, 'NFITPOLY', nparams, 'order of profile parameter fit'

      writespectra, objhdr, plugsort, flux, fluxivar, wset, $
       filebase=filebase

      heap_gc   ; Garbage collection for all lost pointers
   endfor

   splog, 'End time ', systime(1)-t_begin, ' seconds TOTAL', $
    format='(A,F8.2,A)'

end
;------------------------------------------------------------------------------
