;+
; NAME:
;   spreduce
;
; PURPOSE:
;   Extract, wavelength-calibrate, and flatten SDSS spectral frame(s).
;
; CALLING SEQUENCE:
;   spreduce, flatname, arcname, objname, pixflatname, $
;    plugfile=plugfile, lampfile=lampfile, $
;    indir=indir, plugdir=plugdir, outdir=outdir
;
; INPUTS:
;   flatname   - Name of flat-field SDSS image
;   arcname    - Name of arc SDSS image
;   objname    - Name of object SDSS image(s)
;
; OPTIONAL INPUT:
;   pixflatname- Name of pixel-to-pixel flat, produced with SPFLATTEN.
;
; REQUIRED KEYWORDS:
;   plugfile   - Name of plugmap file (Yanny parameter file)
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in the IDL path.
;
; OPTIONAL KEYWORDS:
;   indir      - Input directory for FLATNAME, ARCNAME, OBJNAME;
;                default to './'
;   plugdir    - Input directory for PLUGFILE; default to './'
;   outdir     - Directory for output files; default to './'
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

pro spreduce, flatname, arcname, objname, pixflatname, $
 plugfile=plugfile, lampfile=lampfile, $
 indir=indir, plugdir=plugdir, outdir=outdir

   if (NOT keyword_set(indir)) then indir = './'
   if (NOT keyword_set(plugdir)) then plugdir=indir
   if (NOT keyword_set(outdir)) then outdir = './'

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
         
   readcol, tempname[0], lampwave, lampinten, lampquality, format='f,f,a'
   lamplist = [[lampwave], [lampinten], [(lampquality EQ 'GOOD')]]

   ;---------------------------------------------------------------------------
   ; Read PLUGMAP file
   ;---------------------------------------------------------------------------
 
   fullpath = filepath(plugfile, root_dir=plugdir)
   fullname = findfile(fullpath, count=ct)
   if (ct NE 1) then $
    message, 'Cannot find plugMapFile' + plugMapFile

   yanny_read, fullname[0], pstruct, hdr=hdrplug
   plugmap = *pstruct[0]

   ;---------------------------------------------------------------------------
   ; Read pixel-to-pixel flat-field
   ;---------------------------------------------------------------------------

   if (keyword_set(pixflatname)) then begin
      fullpath = filepath(pixflatname, root_dir=indir)
      fullname = findfile(fullpath, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find pixflat image' + pixflatname
      pixflat = readfits(pixflatname)
   endif

   ;---------------------------------------------------------------------------
   ; Read flat-field image
   ;---------------------------------------------------------------------------

   fullpath = filepath(flatname, root_dir=indir)
   fullname = findfile(fullpath, count=ct)
   if (ct NE 1) then $
    message, 'Cannot find flat image' + flatname
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
    message, 'Cannot find arc image' + arcname
   sdssproc, fullname[0], arcimg, arcivar, hdr=archdr

   ; Flat-field the arc image
   if (keyword_set(pixflatname)) then $
    arcimg = arcimg / pixflat

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
; Now do something with fiberflat !!!???

   ; Construct a vector for each fiber to take out global variations of
   ; the fibers relative to each other.  This needs the wavelength calibration.
   medval = median(flat_flux)
   medsize = 35
   flatglobal = 0 * flat_flux
   ntrace=(size(flat_flux))[2]
   for i=0, ntrace-1 do $
      flatglobal[*,i] = median(flat_flux[*,i], medsize) / medval

   ; Compute the wavelengths for the extracted spectra
   traceset2xy,wset,xx,yy
   waves = 10^yy

   ; Find the wavelength range [WAVEMIN,WAVEMIN] that is in common for
   ; all fibers.  Also find the minimum wavelength separation between pixels.
   ny = (size(waves))[1]
   if (waves[0,0] LT waves[ny-1,0]) then begin
      ; Ascending wavelengths
      wavemin = max(waves[0,*])
      wavemax = min(waves[ny-1,*])
      deltaw = min( waves[1:ny-1,*] - waves[0:ny-2,*] )
   endif else begin
      ; Descending wavelengths
      wavemin = max(waves[ny-1,*])
      wavemax = min(waves[0,*])
      deltaw = min( waves[0:ny-2,*] - waves[1:ny-1,*] )
   endelse

   ; Linearly interpolate all of the flat-field vectors onto a common
   ; wavelength scale
   wtemp = wavemin + findgen(fix(wavemax-wavemin)/deltaw)
   ntrace=(size(flat_flux))[2]
   flattemp = fltarr(N_elements(wtemp),ntrace)
   for i=0, ntrace-1 do $
      flattemp[*,i] = interpol(flat_flux[*,i], waves[*,i], wtemp)

; mm = djs_median(flattemp,2)
; mm = mm / median(mm)
; junk = 0*flattemp
; for i=0, 319 do $
;  junk[*,i]= median( flattemp[*,i]/(median(flattemp[*,i])*mm), 25)
; writefits,'flat_b1.fits',junk
; dfpsplot, 'flat_b1.ps', /square
; djs_plot, wtemp, [0], yr=[0.8,3.0],/ystyle, xtitle='Lambda', $
;  ytitle='Flat-field vectors + const', charsize=2
; for i=0, 300, 20 do $
;  djs_oplot, wtemp, $
;   median( flattemp[*,i]/(median(flattemp[*,i])*mm), 25) +i/200.
; dfpsclose

stop
   ;---------------------------------------------------------------------------
   ; LOOP THROUGH OBJECT FRAMES
   ;---------------------------------------------------------------------------

   for iobj=0, N_elements(objname)-1 do begin

      ;------------------
      ; Read object image

      fullpath = filepath(objname[iobj], root_dir=indir)
      fullname = findfile(fullpath, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find object image' + objname[iobj]
      sdssproc, fullname[0], objimg, objivar, hdr=objhdr

      ; Flat-field the object image
      if (keyword_set(pixflatname)) then $
       objimg = objimg / pixflat

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

      extract_image, objimg, objivar, xsol, sigma, obj_flux, obj_fluxivar, $
       proftype=proftype, wfixed=[1,1,1], $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, whopping=whopping

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
