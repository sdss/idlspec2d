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
   if (NOT keyword_set(qadir)) then outdir = './'

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

;   pixflat = spflatten(flatname, indir=indir, xsol=xsol, ycen=ycen, $
;         flat_flux=flat_flux)

   sigma = 1.0
   proftype = 1 ; Gaussian
   highrej = 20
   lowrej = 25
   nPoly = 4

   extract_image, flatimg, flativar, xsol, sigma, flat_flux, flat_fluxivar, $
    proftype=proftype, wfixed=[1,1,1], $
    highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1, ymodel=ymodel

   ;---------------------------------------------------------------------------
   ; Read the arc
   ;---------------------------------------------------------------------------

   fullpath = filepath(arcname, root_dir=indir)
   fullname = findfile(fullpath, count=ct)
   if (ct NE 1) then $
    message, 'Cannot find arc image ' + arcname
   sdssproc, fullname[0], arcimg, arcivar, hdr=archdr

   ; Flat-field the arc image
   if (keyword_set(pixflatname)) then begin
      arcimg = arcimg / pixflat
      arcivar = arcivar * pixflat^2
   endif

   ;---------------------------------------------------------------------------
   ; Extract the arc image
   ;---------------------------------------------------------------------------

   sigma = 1.0
   proftype = 1 ; Gaussian
   highrej = 10
   lowrej = 15
   nPoly = 4

   extract_image, arcimg, arcivar, xsol, sigma, arc_flux, arc_fluxivar, $
    proftype=proftype, wfixed=[1,1,1], $
    highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1

   ;------------------
   ; Flat-field the extracted arcs with the global flat
   ; Hmmm.... Circular here, since we need a wavelength calibration before
   ; making that flat
;   arc_flux = arc_flux / fflat
;   arc_fluxivar = arc_fluxivar * fflat^2

   ;---------------------------------------------------------------------------
   ; Compute wavelength calibration for arc lamp only
   ;---------------------------------------------------------------------------

   ; Decide if this is a red or blue spectrograph based upon CCDCOL
   ; in the arc header, which was written to the header by SDSSPROC.
   camcol = sxpar(archdr, 'CAMCOL')
   if (camcol EQ 1 OR camcol EQ 3) then color='blue' $
    else if (camcol EQ 4 OR camcol EQ 2) then color='red' $
    else message, 'No CAMCOL keyword in arc header'

   fitarcimage, arc_flux, arc_fluxivar, color, lamplist, xpeak, ypeak, wset, $
    invset, ans=wavesolution, lambda=lambda, goodlines=goodlines

   ;---------------------------------------------------------------------------
   ; Compute fiber-to-fiber flat-field variations
   ;---------------------------------------------------------------------------

   fflat = fiberflat(flat_flux, wset)

   ;---------------------------------------------------------------------------
   ; LOOP THROUGH OBJECT FRAMES
   ;---------------------------------------------------------------------------

   for iobj=0, N_elements(objname)-1 do begin

      ;------------------
      ; Read object image

      fullpath = filepath(objname[iobj], root_dir=indir)
      fullname = findfile(fullpath, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find object image ' + objname[iobj]
      sdssproc, fullname[0], objimg, objivar, hdr=objhdr

      ; What about invvar???  Is below okay???
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

      ;------------------
      ; Extract the object image
      ; Use the "whopping" terms
      ; We need to do 2 iteration extraction: 
      ;        1) Fit profiles in a subset of rows
      ;        2) Fit returned parameters with smooth functions
      ;        3) Extract all 2048 rows with new profiles given by
      ;              fitansimage

      nrow = (size(objimg))[2]
      ncol = (size(objimg))[1]
      skiprow = 8
      yrow = lindgen(nrow/skiprow)*skiprow + skiprow/2
      nfirst = n_elements(yrow)

      ; 1) First extraction
      extract_image, objimg, objivar, xsol, sigma, obj_flux, obj_fluxivar, $
       proftype=proftype, wfixed=[1,1,1], yrow=yrow, $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping, $
       ansimage=ansimage

      ; 2) Refit ansimage to smooth profiles

      nparams = 3
      nTrace = (size(obj_flux))[2]
      fitans = fitansimage(ansimage, nparams, nTrace, nPoly, nfirst, yrow, $
             fluxm = [1,1,0])

      ; 3) Second and final extraction
      extract_image, objimg, objivar, xsol, sigma, obj_flux, obj_fluxivar, $
       proftype=proftype, wfixed=[1,1,1], fitans=fitans, $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping, $
       ymodel=ymodel2

      ;------------------
      ; Flat-field the extracted object fibers with the global flat
      obj_flux = obj_flux / fflat
      obj_fluxivar = obj_fluxivar * fflat^2
stop

      ;------------------
      ; Tweak up the wavelength solution to agree with the sky lines.

      ; ???

      ;------------------
      ; Sky-subtract

      ; ???

      ;------------------
      ; Write extracted, lambda-calibrated, sky-subtracted spectra to disk

      ; ???

   endfor

stop
end
;------------------------------------------------------------------------------
