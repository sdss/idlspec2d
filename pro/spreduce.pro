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
;   flatname   - Name of flat-field SDSS image(s)
;   arcname    - Name of arc SDSS image(s)
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
;                default to './'
;   plugdir    - Input directory for PLUGFILE; default to './'
;   outdir     - Directory for output files; default to './'
;   qadir      - Directory for QA files; default to './'
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
;
; PROCEDURES CALLED:
;   djs_locate_file()
;   djs_median()
;   extract_boxcar()
;   extract_image
;   readcol
;   readfits()
;   sdssproc
;   traceset2xy
;   xy2traceset
;   yanny_read
;
; REVISION HISTORY:
;   12-Oct-1999  Written by D. Schlegel & S. Burles, APO
;-
;------------------------------------------------------------------------------

pro spreduce, flatname, arcname, objname, pixflatname=pixflatname, $
 plugfile=plugfile, lampfile=lampfile, $
 indir=indir, plugdir=plugdir, outdir=outdir, qadir=qadir, qa=qa

   if (NOT keyword_set(indir)) then indir = './'
   if (NOT keyword_set(plugdir)) then plugdir=indir
   if (NOT keyword_set(outdir)) then outdir = './'
   if (NOT keyword_set(qadir)) then qadir = outdir

   t_begin = systime(1)

   ;---------------------------------------------------------------------------
   ; Read LAMPLIST file for wavelength calibration
   ;---------------------------------------------------------------------------

   if (keyword_set(lampfile)) then begin
      lampfilenames = findfile(lampfile, count=ct)
      if (ct EQ 0) then message, 'No LAMPFILE found '+lampfile
   endif else begin
      lampdefault = 'lamphgcdne.dat'
      lampfilenames = djs_locate_file(lampdefault)
      if (lampfilenames EQ '') then message, 'No LAMPFILE found '+lampdefault
   endelse
         
   readcol, lampfilenames[0], lampwave, lampinten, lampquality, format='d,f,a'
   lamplist = [[lampwave], [lampinten], [(lampquality EQ 'GOOD')]]

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
   ptr_free, pstruct  ;;!!!I don't know if this destroys the pointer totally

   ;---------------------------------------------------------------------------
   ; Read flat-field image
   ;---------------------------------------------------------------------------

   iflat = 0
   iarc = 0
   arcdone = 0
   numarcfiles = (size(arcname))[1]
   numflatfiles = (size(flatname))[1]
   while (arcdone EQ 0 AND iarc LT numarcfiles AND $
         iflat LT numflatfiles) do begin

     flatpath = filepath(flatname[iflat], root_dir=indir)
     flatfilenames = findfile(flatpath, count=ct)
     if (ct NE 1) then $
      message, 'Cannot find flat image ' + flatname[iflat]

   print, 'Reading in Flat ', flatfilenames[0]
   sdssproc, flatfilenames[0], flatimg, flativar, hdr=flathdr, $
    pixflatname=pixflatname
 
   ;---------------------------------------------------------------------------
   ; Create spatial tracing from flat-field image
   ;---------------------------------------------------------------------------

   xcen = trace320crude(flatimg, yset=ycen, maxdev=0.15)

   xy2traceset, ycen, xcen, tset, ncoeff=5, maxdev=0.1
   traceset2xy, tset, ycen, xsol

   ;---------------------------------------------------------------------------
   ; Extract the flat-field image
   ;---------------------------------------------------------------------------

   print, 'Extracting flat-field with simple gaussian'
   sigma = 1.0
   proftype = 1 ; Gaussian
   highrej = 20
   lowrej = 25
   nPoly = 6
   wfixed = [1] ; Just fit the first gaussian term

   extract_image, flatimg, flativar, xsol, sigma, flat_flux, flat_fluxivar, $
    proftype=proftype, wfixed=wfixed, $
    highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1

   highpixels = where(flat_flux GT 1.0e5,numhighpixels)

   print, 'Found ', numhighpixels, ' highpixels in extracted flatfield ', $
       flatname[iflat]

   flatimg = 0
   flativar = 0

   ;---------------------------------------------------------------------------
   ; Compute fiber-to-fiber flat-field variations
   ;---------------------------------------------------------------------------

   fflat = fiberflat(flat_flux, flat_fluxivar)

   flat_flux=0
   flat_fluxivar=0
   ;---------------------------------------------------------------------------
   ; Read the arc
   ;---------------------------------------------------------------------------

     arcpath = filepath(arcname[iarc], root_dir=indir)
     arcfilenames = findfile(arcpath, count=ct)
     if (ct NE 1) then $
      message, 'Cannot find arc image ' + arcname[iarc]

     print, 'Reading in Arc ', arcfilenames[0]
     sdssproc, arcfilenames[0], arcimg, arcivar, hdr=archdr, $
      pixflatname=pixflatname, spectrographid=spectrographid, color=color

     ;------------------------------------------------------------------------
     ; Do stats here to make sure it's an arc


     
     ;--------------------------------------------------------------------------
     ; Extract the arc image
     ;--------------------------------------------------------------------------

     print, 'Extracting arc-lamp with simple gaussian'
     sigma = 1.0
     proftype = 1 ; Gaussian
     highrej = 10
     lowrej = 15
     nPoly = 6 ; maybe more structure
     wfixed = [1] ; Just fit the first gaussian term

     extract_image, arcimg, arcivar, xsol, sigma, arc_flux, arc_fluxivar, $
      proftype=proftype, wfixed=wfixed, $
      highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1

     arcimg = 0
     arcivar = 0
     ;------------------
     ; Flat-field the extracted arcs with the global flat
     ; Hmmm.... would be circular if we need a wavelength calibration before
     ; making that flat.  We don't at the moment.

     arc_flux = arc_flux / fflat
     arc_fluxivar = arc_fluxivar * fflat^2

     ;-------------------------------------------------------------------------
     ; Compute wavelength calibration for arc lamp only
     ;-------------------------------------------------------------------------

     print, 'Searching for wavelength solution with fitarcimage'
     arcstatus = fitarcimage(arc_flux, arc_fluxivar, color, lamplist, $
      xpeak, ypeak, wset, invset, lambda=lambda, $
      xdif_lfit=xdif_lfit, xdif_tset=xdif_tset)

     if (arcstatus) then begin
       iarc = iarc + 1
       iflat = iflat + 1
       print, 'Trying the next arc ', arcname[iarc]
       print, '   and next flat ', flatname[iflat]
     endif else arcdone = 1
     
   endwhile

   if (arcdone NE 1) then $
     message, 'Did not find a good arc solution'

   qaplot_arcline, xdif_tset, lambda, arcname[iarc]

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

      print,'> SPREDUCE: ',systime(1)-t_begin, ' seconds so far', $
       format='(A,F8.2,A)'

      ;------------------
      ; Read object image

      objpath = filepath(objname[iobj], root_dir=indir)
      objfilenames = findfile(objpath, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find object image ' + objname[iobj]

      objfile = objfilenames[0]
     print, 'Reading in Object ', objfile
      sdssproc, objfile, objimg, objivar, hdr=objhdr, $
       pixflatname=pixflatname, spectrographid=spectrographid, color=color

      ;------------------
      ; Tweak up the spatial traces

      ; ???

      ;------------------
      ; Identify very bright objects
      ; Do a boxcar extraction, and look for fibers where the median
      ; counts are 10000 ADU per row.

      fextract = extract_boxcar(objimg, xsol, ycen)
      scrunch = djs_median(fextract, 1) ; Find median counts/row in each fiber
      whopping = where(scrunch GT 10000.0, whopct)
      print, 'Number of bright fibers = ', whopct

; QUICK EXTRACTION...
;extract_image, objimg, objivar, xsol, sigma, obj_flux, obj_fluxivar, $
; proftype=proftype, wfixed=[1,1,1], $
; highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping

      ;------------------
      ; Extract the object image
      ; Use the "whopping" terms
      ; We need to do 2 iteration extraction: 
      ;        1) Fit profiles in a subset of rows
      ;        2) Fit returned parameters with smooth functions
      ;        3) Extract all 2048 rows with new profiles given by
      ;              fitansimage

      print, 'Extracting frame '+objname[iobj]+' with 6 step process'
      nrow = (size(objimg))[2]
      ncol = (size(objimg))[1]
      skiprow = 8
      yrow = lindgen(nrow/skiprow)*skiprow + skiprow/2
      nfirst = n_elements(yrow)

      osigma = 1.0
      sigma = xsol*0.0 + osigma
      proftype = 1 ; Gaussian
      highrej = 10
      lowrej = 15
      nPoly = 6 ; maybe more structure
      wfixed = [1,1,1] ; gaussian term + centroid and  sigma terms
      nTerms = 3
      sigmaterm = 1
      centerterm = 2
      xnow = xsol
      sigmanow = sigma

      for i = 0, 1 do begin
      
        ; 1) First extraction
        print, 'Object extraction: Step', i*3+1
        extract_image, objimg, objivar, xnow, sigma, obj_flux,obj_fluxivar, $
         proftype=proftype, wfixed=wfixed, yrow=yrow, $
         highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping, $
         ansimage=ansimage

        ; 2) Refit ansimage to smooth profiles

        print, 'Answer Fitting: Step', i*3+2
        nparams = 3
        nTrace = (size(obj_flux))[2]
        fitans = fitansimage(ansimage, nparams, nTrace, nPoly, nfirst, yrow, $
         fluxm = [1,1,0], crossfit=1-i)

        ; 3) Calculate new sigma and xsol arrays
      
        if (i EQ 0) then begin 
           print, 'Trace Tweaking: Step', i*3+3, '    (Sigma not tweaked)'
           sigmashift = transpose(fitans[lindgen(nTrace)*nTerms + sigmaterm, *])
           centershift = $
              transpose(fitans[lindgen(nTrace)*nTerms + centerterm, *])
           tweaktrace, xnow, sigmanow, centershift, sigmashift
        endif
      endfor

      ; 4) Second and final extraction
      print, 'Object extraction: Step 6'
;
;	Using old sigma for now, which should be fine
;	Different sigmas require a new profile for each trace, so will
;	check timing in the future
;
      extract_image, objimg, objivar, xnow, sigma, obj_flux, $
       obj_fluxivar, proftype=proftype, wfixed=wfixed, fitans=fitans, $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping 
;       ymodel=ymodel2

      ;------------------
      ; Flat-field the extracted object fibers with the global flat
      obj_flux = obj_flux / fflat
      obj_fluxivar = obj_fluxivar * fflat^2

      ;------------------
      ; Tweak up the wavelength solution to agree with the sky lines.

      locateskylines, skylinefile, obj_flux, obj_fluxivar, $
       wset, invset, xsky, ysky, skywaves

; DO NOT TWEAK THE SKY LINES -- THIS ROUTINE MESSES UP NEAR ROW 145 ???
; COMMENT OUT locateskylines, fit_skyset

       wset_tweak = wset
       invset_tweak = invset

      ;
      ;	First convert lambda, and skywaves to log10 vacuum
      ;
      print, 'converting wavelengths to vacuum'
	vaclambda = 10^lambda
        airtovac, vaclambda
	vaclambda = alog10(vaclambda)

	vacsky = skywaves
        airtovac, vacsky
	vacsky = alog10(vacsky)

	sxaddpar, hdr, 'VACUUM', 'WAVELENGTHS ARE IN VACUUM'
	sxaddpar, hdr, 'AIR2VAC', systime()
      print, 'now tweaking to sky lines'
      skycoeff = 2
      fit_skyset, xpeak, ypeak, vaclambda, xsky, ysky, vacsky, skycoeff, $
        goodlines, wset_tweak, invset_tweak, ymin=ymin, ymax=ymax, func=func

      ;------------------
      ; Sky-subtract

      plugsort = sortplugmap(plugmap, spectrographid)

      skysubtract, obj_flux, obj_fluxivar, plugsort, wset_tweak, $ 
       skysub, skysubivar

      obj_flux =0
      obj_fluxivar =0
      ;------------------------------------------
      ; Flux calibrate to spectrophoto_std fibers

      fluxfactor = fluxcorr(skysub, skysubivar, wset_tweak, plugsort, $
                             lower=2.5, upper=10)

      fluxout = skysub * fluxfactor
      fluxoutivar = skysubivar / (fluxfactor^2)

      ;------------------------------------------
      ; Telluric correction called for 'red' side

      if (color EQ 'red')  then begin

        telluricfactor = telluric_corr(fluxout, fluxoutivar, plugsort)
	fluxout = fluxout / telluricfactor
	fluxoutivar = fluxoutivar * (telluricfactor^2)

      endif

      ;------------------
      ; Write extracted, lambda-calibrated, sky-subtracted spectra to disk

      framenum = sxpar(objhdr, 'EXPOSURE')

      filebase = filepath( $
       's-'+string(format='(i1,a1,a,i4.4)',spectrographid, $
       color,'-',framenum), root_dir=outdir)

;
;	Add everything we can think of to object header
;
	sxaddpar, objhdr, 'PLUGMAPF', plugfilenames[0]
	sxaddpar, objhdr, 'FLATFILE', flatfilenames[0]
	sxaddpar, objhdr, 'ARCFILE',  arcfilenames[0]
	sxaddpar, objhdr, 'OBJFILE',  objfile
	sxaddpar, objhdr, 'LAMPLIST',  lampfilenames[0]
	sxaddpar, objhdr, 'SKYLIST',  skylinefile
	sxaddpar, objhdr, 'PIXFLAT',  pixflatname
	sxaddpar, objhdr, 'OSIGMA',  osigma, $
            'Original guess at sigma of spatial profiles'
	sxaddpar, objhdr, 'SKIPROW', skiprow, 'Number of rows skipped in step 1'
	sxaddpar, objhdr, 'LOWREJ', lowrej, 'Extraction, low rejection'
	sxaddpar, objhdr, 'HIGHREJ', highrej, 'Extraction, high rejection'
	sxaddpar, objhdr, 'SCATPOLY', nPoly, 'Order of scattered light poly'
	sxaddpar, objhdr, 'PROFTYPE', proftype, '1 is Gaussian'
	sxaddpar, objhdr, 'NFITPOLY', nparams, 'order of profile parameter fit'


      writespectra, objhdr, plugsort, fluxout, fluxoutivar, wset_tweak, $
       filebase=filebase

;
;	Clear out variables for memory considerations
;
     objimg = 0
     objivar = 0
     skysub = 0
     skysubivar = 0
     fluxout = 0
     fluxoutivar = 0
     fluxfactor = 0
     telluricfactor = 0
     wset_tweak = 0
     invset_tweak = 0
     sigma = 0
     xnow = 0
     ansimage = 0
     fitans = 0
     
     heap_gc   ; Garbage collection for all lost pointers
   endfor

   print,'> SPREDUCE: ', systime(1)-t_begin, ' seconds TOTAL', $
    format='(A,F8.2,A)'

end
;------------------------------------------------------------------------------
