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
      tempname = findfile(lampfile, count=ct)
      if (ct EQ 0) then message, 'No LAMPFILE found '+lampfile
   endif else begin
      lampdefault = 'lamphgcdne.dat'
      tempname = djs_locate_file(lampdefault)
      if (tempname EQ '') then message, 'No LAMPFILE found '+lampdefault
   endelse
         
   readcol, tempname[0], lampwave, lampinten, lampquality, format='d,f,a'
   lamplist = [[lampwave], [lampinten], [(lampquality EQ 'GOOD')]]

   ;---------------------------------------------------------------------------
   ; Locate skyline file for sky wavelength calibration
   ;---------------------------------------------------------------------------

   if (keyword_set(skylinefile)) then begin
      tempname = findfile(skylinefile, count=ct)
      if (ct EQ 0) then message, 'No SKYLINEFILE found '+lampfile
   endif else begin
      skydefault = 'skylines.dat'
      tempname = djs_locate_file(skydefault)
      if (tempname EQ '') then message, 'No SKYLINEFILE found '+skydefault
   endelse

   skylinefile = tempname[0]

   ;---------------------------------------------------------------------------
   ; Read PLUGMAP file
   ;---------------------------------------------------------------------------
 
   fullpath = filepath(plugfile, root_dir=plugdir)
   fullname = findfile(fullpath, count=ct)
   if (ct NE 1) then $
    message, 'Cannot find plugMapFile ' + plugMapFile

   yanny_read, fullname[0], pstruct, hdr=hdrplug
   plugmap = *pstruct[0]

   ;---------------------------------------------------------------------------
   ; Read pixel-to-pixel flat-field
   ;---------------------------------------------------------------------------

   if (keyword_set(pixflatname)) then begin
      fullpath = filepath(pixflatname, root_dir=indir)
      fullname = findfile(fullpath, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find pixflat image ' + pixflatname
      pixflat = readfits(fullname[0])
   endif

   ;---------------------------------------------------------------------------
   ; Read flat-field image
   ;---------------------------------------------------------------------------

   fullpath = filepath(flatname, root_dir=indir)
   fullname = findfile(fullpath, count=ct)
   if (ct NE 1) then $
    message, 'Cannot find flat image ' + flatname
   sdssproc, fullname[0], flatimg, flativar, hdr=flathdr
 
   ; Flat-field the flat image
   if (keyword_set(pixflatname)) then $
    flatimg = flatimg / pixflat

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
   nPoly = 4
   wfixed = [1] ; Just fit the first gaussian term

   extract_image, flatimg, flativar, xsol, sigma, flat_flux, flat_fluxivar, $
    proftype=proftype, wfixed=wfixed, $
    highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1, ymodel=ymodel

   ;---------------------------------------------------------------------------
   ; Compute fiber-to-fiber flat-field variations
   ;---------------------------------------------------------------------------

   fflat = fiberflat(flat_flux, flat_fluxivar)

   ;---------------------------------------------------------------------------
   ; Read the arc
   ;---------------------------------------------------------------------------

   fullpath = filepath(arcname, root_dir=indir)
   fullname = findfile(fullpath, count=ct)
   if (ct NE 1) then $
    message, 'Cannot find arc image ' + arcname
   sdssproc, fullname[0], arcimg, arcivar, hdr=archdr, $
    spectrographid=spectrographid, color=color

   ; Flat-field the arc image
   if (keyword_set(pixflatname)) then begin
      arcimg = arcimg / pixflat
      arcivar = arcivar * pixflat^2
   endif

   ;---------------------------------------------------------------------------
   ; Extract the arc image
   ;---------------------------------------------------------------------------

   print, 'Extracting arc-lamp with simple gaussian'
   sigma = 1.0
   proftype = 1 ; Gaussian
   highrej = 10
   lowrej = 15
   nPoly = 4
   wfixed = [1] ; Just fit the first gaussian term

   extract_image, arcimg, arcivar, xsol, sigma, arc_flux, arc_fluxivar, $
    proftype=proftype, wfixed=wfixed, $
    highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1

   ;------------------
   ; Flat-field the extracted arcs with the global flat
   ; Hmmm.... would be circular if we need a wavelength calibration before
   ; making that flat.  We don't at the moment.

   arc_flux = arc_flux / fflat
   arc_fluxivar = arc_fluxivar * fflat^2

   ;---------------------------------------------------------------------------
   ; Compute wavelength calibration for arc lamp only
   ;---------------------------------------------------------------------------

   print, 'Searching for wavelength solution with fitarcimage'

   fitarcimage, arc_flux, arc_fluxivar, color, lamplist, xpeak, ypeak, wset, $
    invset, ans=wavesolution, lambda=lambda, $
    xdif_lfit=xdif_lfit, xdif_tset=xdif_tset, errcode=errcode

   if (errcode NE 0) then begin
      message, 'Fitarcimage failed'
   endif 

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

      print,'> SPREDUCE: ',systime(1)-t_begin, ' seconds so far', $
       format='(A,F8.2,A)'

      ;------------------
      ; Read object image

      fullpath = filepath(objname[iobj], root_dir=indir)
      fullname = findfile(fullpath, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find object image ' + objname[iobj]
      sdssproc, fullname[0], objimg, objivar, hdr=objhdr, $
       spectrographid=spectrographid, color=color

      if (keyword_set(pixflatname)) then begin
         objimg = objimg / pixflat
         objivar = objivar * pixflat^2
      endif

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

      print, 'Extracting frame '+objname[iobj]+' with 3 step process'
      nrow = (size(objimg))[2]
      ncol = (size(objimg))[1]
      skiprow = 8
      yrow = lindgen(nrow/skiprow)*skiprow + skiprow/2
      nfirst = n_elements(yrow)

; COMMENT OUT FOR NOW SINCE extract_image crashes on step 3 ???
      ; 1) First extraction
      print, 'Object extraction: Step 1'
      extract_image, objimg, objivar, xsol, sigma, obj_flux, obj_fluxivar, $
       proftype=proftype, wfixed=[1,1,1], yrow=yrow, $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping, $
       ansimage=ansimage

      ; 2) Refit ansimage to smooth profiles

      print, 'Object extraction: Step 2'
      nparams = 3
      nTrace = (size(obj_flux))[2]
      fitans = fitansimage(ansimage, nparams, nTrace, nPoly, nfirst, yrow, $
       fluxm = [1,1,0])

;
;	shiftfit, fitans, nTrace, xsol, sigma, xsolout, sigmaout
;      

      ; 3) Second and final extraction
      print, 'Object extraction: Step 3'
      extract_image, objimg, objivar, xsol, sigma, obj_flux, obj_fluxivar, $
       proftype=proftype, wfixed=[1,1,1], fitans=fitans, $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping, $
       ymodel=ymodel2

      ;------------------
      ; Flat-field the extracted object fibers with the global flat
      obj_flux = obj_flux / fflat
      obj_fluxivar = obj_fluxivar * fflat^2

      ;------------------
      ; Tweak up the wavelength solution to agree with the sky lines.

      locateskylines, skylinefile, obj_flux, obj_fluxivar, $
       wset, invset, wset_tweak, invset_tweak, $
       xsky, ysky, skywaves

      ;------------------
      ; Sky-subtract

      plugsort = sortplugmap(plugmap, spectrographid)

;      skysubtract, flat1, flat1ivar, plugsort, wset, skysub, skysubivar, $
;                     allwave=allwave, allsky=allsky, allfit=allskyfit

      skysubtract, obj_flux, obj_fluxivar, plugsort, wset_tweak, $ 
       skysub, skysubivar

      ;------------------
      ; Write extracted, lambda-calibrated, sky-subtracted spectra to disk

      framenum = sxpar(objhdr, 'EXPOSURE')

      filebase = filepath( $
       's-'+string(format='(i1,a1,a,i4.4)',spectrographid, $
       color,'-',framenum), root_dir=outdir)

      writespectra, objhdr, plugsort, skysub, skysubivar, wset_tweak, $
       filebase=filebase

   endfor

   print,'> SPREDUCE: ', systime(1)-t_begin, ' seconds TOTAL', $
    format='(A,F8.2,A)'

stop
end
;------------------------------------------------------------------------------
