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
;   Should test that arcs and flats are valid images with CHECKFLAVOR.
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
   bestcorr = -1.0

   for ifile=0, (nflat<narc)-1 do begin

      splog, ifile+1, (nflat<narc), $
       format='("Looping through flat+arc pair #",I3," of",I3)'

      ;------------------------------------------------------------------------
      ; Read flat-field image
      ;------------------------------------------------------------------------

      splog, 'Reading flat ', flatname[ifile]
      sdssproc, flatname[ifile], image, invvar, indir=indir, $
       hdr=flathdr, pixflatname=pixflatname

      ;------------------------------------------------------------------------
      ; Create spatial tracing from flat-field image
      ;------------------------------------------------------------------------

      splog, 'Tracing 320 fibers in ',  flatname[ifile]
      tmp_xsol = trace320crude(image, yset=ycen, maxdev=0.15)

      splog, 'Fitting traces in ',  flatname[ifile]
      xy2traceset, ycen, tmp_xsol, tset, ncoeff=5, maxdev=0.1
      traceset2xy, tset, ycen, tmp_xsol

      ;------------------------------------------------------------------------
      ; Extract the flat-field image
      ;------------------------------------------------------------------------

      splog, 'Extracting flat-field image with simple gaussian'
      sigma = 1.0
      proftype = 1 ; Gaussian
      highrej = 20
      lowrej = 25
      nPoly = 6
      wfixed = [1] ; Just fit the first gaussian term

      extract_image, image, invvar, tmp_xsol, sigma, flux, fluxivar, $
       proftype=proftype, wfixed=wfixed, $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1

      highpixels = where(flux GT 1.0e5, numhighpixels)

      splog, 'Found ', numhighpixels, ' highpixels in extracted flat ', $
       flatname[ifile]

      ;------------------------------------------------------------------------
      ; Compute fiber-to-fiber flat-field variations
      ;------------------------------------------------------------------------

      tmp_fflat = fiberflat(flux, fluxivar, fibermask)

      ;------------------------------------------------------------------------
      ; Read the arc
      ;------------------------------------------------------------------------

      splog, 'Reading arc ', arcname[ifile]
      sdssproc, arcname[ifile], image, invvar, indir=indir, $
       hdr=archdr, pixflatname=pixflatname

      ;------------------------------------------------------------------------
      ; Extract the arc image
      ;------------------------------------------------------------------------

      splog, 'Extracting arc image with simple gaussian'
      sigma = 1.0
      proftype = 1 ; Gaussian
      highrej = 10
      lowrej = 15
      nPoly = 6 ; maybe more structure
      wfixed = [1] ; Just fit the first gaussian term

      extract_image, image, invvar, tmp_xsol, sigma, flux, fluxivar, $
       proftype=proftype, wfixed=wfixed, $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1

      ;------------------
      ; Flat-field the extracted arcs with the global flat
      ; Hmmm.... would be circular if we need a wavelength calibration before
      ; making that flat.  We don't at the moment.

      divideflat, flux, fluxivar, tmp_fflat, fibermask

      ;-------------------------------------------------------------------------
      ; Compute correlation coefficient for this arc image
      ;-------------------------------------------------------------------------

      splog, 'Searching for wavelength solution'
      tmp_aset = 0
      fitarcimage, flux, fluxivar, aset=tmp_aset, $
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
         fflat = tmp_fflat
         aset = tmp_aset
      endif

   endfor

   ;---------------------------------------------------------------------------
   ; Make sure that the best flat+arc pair is good enough
   ;---------------------------------------------------------------------------

   splog, 'Best flat = ', flatname[ibest]
   splog, 'Best arc = ', arcname[ibest]

   if (bestcorr LT 0.7) then begin
      splog, 'Best arc correlation = ', bestcorr
      splog, 'Abort!'
      return
   endif

   ;---------------------------------------------------------------------------
   ; Compute wavelength calibration for arc lamp only
   ;---------------------------------------------------------------------------

   arccoeff = 5

   splog, 'Searching for wavelength solution'
   fitarcimage, arcimg, arcivar, xpeak, ypeak, wset, ncoeff=arccoeff, $
    aset=aset, $
    color=color, lampfile=lampfile, lambda=lambda, xdif_tset=xdif_tset
     
   wsave = wset

   qaplot_arcline, xdif_tset, lambda, arcname[ibest]

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

      ;
      ;	My fitansimage is not going to work well without good profiles
      ; Using it to tweak up traces, and then performing free fit to
      ; object
      ;
      for i = 0, 0 do begin
      
        ; 1) First extraction
        splog, 'Object extraction: Step', i*3+1
        extract_image, image, invvar, xnow, sigma, tempflux, tempfluxivar, $
         proftype=proftype, wfixed=wfixed, yrow=yrow, $
         highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping, $
         ansimage=ansimage

        ; 2) Refit ansimage to smooth profiles

        splog, 'Answer Fitting: Step', i*3+2
        nparams = 3
        nTrace = (size(tempflux))[2]
        fitans = fitansimage(ansimage, nparams, nTrace, nPoly, nfirst, yrow, $
         fluxm = [1,1,0], crossfit=1-i, scatfit=scatfit, scatimage=scatimage)

        ; 3) Calculate new sigma and xsol arrays
      
        if (i EQ 0) then begin 
           splog, 'Trace Tweaking: Step', i*3+3, '    (Sigma not tweaked)'
           sigmashift = transpose(fitans[lindgen(nTrace)*nTerms + sigmaterm, *])
           centershift = $
            transpose(fitans[lindgen(nTrace)*nTerms + centerterm, *])
           tweaktrace, xnow, sigmanow, centershift, sigmashift
        endif
      endfor

      splog, 'Skipping steps 4 and 5'
      splog, 'Scattered light: median ', median(scatfit), ' electrons'

      ; 4) Second and final extraction
      splog, 'Object extraction: Step 6'

      ; Using old sigma for now, which should be fine
      ; Different sigmas require a new profile for each trace, so will
      ; check timing in the future
      ; subtract off fit to scattered light and don't allow any polynomial terms

      extract_image, image-scatfit, invvar, xnow, sigma, flux, $
       fluxivar, proftype=proftype, wfixed=wfixed, $
       highrej=highrej, lowrej=lowrej, nPoly=0, whopping=whopping, chisq=chisq


;      fitans = fitans[0:nparams*nTrace-1,*]
;      extract_image, image-scatfit, invvar, xnow, sigma, flux, $
;       fluxivar, proftype=proftype, wfixed=wfixed, fitans=fitans, $
;       highrej=highrej, lowrej=lowrej, nPoly=0, whopping=whopping, chisq=chisq
;       ymodel=ymodel2

      plot, chisq

      ;------------------
      ; Flat-field the extracted object fibers with the global flat
      divideflat, flux, fluxivar, fflat, fibermask

      ;------------------
      ; Tweak up the wavelength solution to agree with the sky lines.
      ; xshift contains polynomial coefficients to shift arc to sky lines.

      locateskylines, skylinefile, flux, fluxivar, $
       wset, xsky, ysky, skywaves, xcoeff=xcoeff

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

      ;
      ;	Fit to arc lines with sky shifts included
      ;


    
      xshift = double(xpeak) 
      wold = wset
      for i=0,nTrace - 1 do xshift[i,*] =  poly(xpeak[i,*],xcoeff[i,*])

       plot, xpeak, xshift, ps=3 

      vacset = wset
      fixabove = 2
      finalarcfit, xpeak+xshift, vacloglam, vacset, arccoeff, fixabove, $
             maxdev=0, maxiter=1, nsetcoeff=8, maxsig=2.0


;      fit_skyset, xpeak, ypeak, vacloglam, xsky, ysky, vaclogsky, skycoeff, $
;        goodlines, wset, ymin=ymin, ymax=ymax, func=func

      ;------------------
      ; Sky-subtract

      skysubtract, flux, fluxivar, plugsort, vacset, $ 
       skysub, skysubivar, fibermask=fibermask

      ;------------------------------------------
      ; Flux calibrate to spectrophoto_std fibers

      fluxfactor = fluxcorr(skysub, skysubivar, vacset, plugsort, $
                             lower=1.5, upper=5, fibermask=fibermask)

      flux = skysub * fluxfactor
      fluxivar = skysubivar / (fluxfactor^2)

      ;------------------------------------------
      ; Telluric correction called for 'red' side

      if (color EQ 'red')  then begin

         telluricfactor = telluric_corr(flux, fluxivar, vacset, plugsort)
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

      sxaddpar, objhdr, 'PLUGMAPF', plugfilename
      sxaddpar, objhdr, 'FLATFILE', flatname[ibest]
      sxaddpar, objhdr, 'ARCFILE', arcname[ibest]
      sxaddpar, objhdr, 'OBJFILE', objname[iobj]
      sxaddpar, objhdr, 'LAMPLIST', lampfile
      sxaddpar, objhdr, 'SKYLIST', skylinefile
      sxaddpar, objhdr, 'PIXFLAT', pixflatname
      sxaddpar, objhdr, 'OSIGMA',  sigma, $
           'Original guess at sigma of spatial profiles'
      sxaddpar, objhdr, 'SKIPROW', skiprow, 'Number of rows skipped in step 1'
      sxaddpar, objhdr, 'LOWREJ', lowrej, 'Extraction, low rejection'
      sxaddpar, objhdr, 'HIGHREJ', highrej, 'Extraction, high rejection'
      sxaddpar, objhdr, 'SCATPOLY', nPoly, 'Order of scattered light poly'
      sxaddpar, objhdr, 'PROFTYPE', proftype, '1 is Gaussian'
      sxaddpar, objhdr, 'NFITPOLY', nparams, 'order of profile parameter fit'

      writespectra, objhdr, flux, fluxivar, plugsort, vacset, $
       filebase=filebase

      heap_gc   ; Garbage collection for all lost pointers
   endfor

   splog, 'End time ', systime(1)-t_begin, ' seconds TOTAL', $
    format='(A,F8.2,A)'

   return
end
;------------------------------------------------------------------------------
