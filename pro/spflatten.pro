;+
; NAME:
;   spflatten
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   spflatten, flatname, [ pixflat, sigrej=, maxiter=, $
;    outfile=, indir=, outdir=, tmpdir= ]
;
; INPUTS:
;   flatname   - Name of flat-field SDSS image(s)
;
; OPTIONAL INPUTS:
;   sigrej     - Sigma rejection level; default to 1, 1, 1.1, 1.3, 1.6 or 1.9
;                for 1,2,3,4,5 or 6 flats.  For more then 6 flats, default
;                to 2.0.
;   maxiter    - Number of rejection iterations; default to 2.
;   outfile    - Write the image PIXFLAT to this file.
;   indir      - Input directory for FLATNAME; default to './'
;   outdir     - Output directory for OUTFILE; default to './'
;   tmpdir     - Directory for temporary files; default to same as OUTDIR
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   pixflat    - Image containing all the information about pixel-to-pixel
;                variations.  Illumination variations are removed.
;
; COMMENTS:
;   This program writes 2*nflat temporary files to disk to save internal memory.
;   But it still creates an two arrays (FLATARR (float) and OUTMASK (byte))
;   that is as large as all of the input flats.
;   Thus, if you are using ten 16 MB flat images,
;   that array will be 10*(16+4)=200 MB (in addition to a few other images).
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_avsigclip()
;   extract_image
;   readfits()
;   sdssproc
;   writefits
;
; REVISION HISTORY:
;   13-Oct-1999  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------

pro spflatten, flatname, pixflat, sigrej=sigrej, maxiter=maxiter, $
 outfile=outfile, indir=indir, outdir=outdir, tmpdir=tmpdir

   if (NOT keyword_set(indir)) then indir = './'
   if (NOT keyword_set(outdir)) then outdir = './'
   if (NOT keyword_set(tmpdir)) then tmpdir = outdir

   nflat = N_elements(flatname)
   ngrow = 2

   if (NOT keyword_set(sigrej)) then begin
      if (nflat LE 2) then sigrej = 1.0 $ ; Irrelevant for only 1 or 2 flats
       else if (nflat EQ 3) then sigrej = 1.1 $
       else if (nflat EQ 4) then sigrej = 1.3 $
       else if (nflat EQ 5) then sigrej = 1.6 $
       else if (nflat EQ 6) then sigrej = 1.9 $
       else sigrej = 2.0
   endif
   if (NOT keyword_set(maxiter)) then maxiter = 2

   tmpname1 = tmpdir+'/tmp.flatimg.'+strtrim(string(indgen(nflat)),2)+'.fits'
   tmpname2 = tmpdir+'/tmp.ymodel.'+strtrim(string(indgen(nflat)),2)+'.fits'

   for iflat=0, nflat-1 do begin

      ;----------------------
      ; Read flat-field image

      fullpath = filepath(flatname[iflat], root_dir=indir)
      fullname = findfile(fullpath, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find flat image' + flatname[iflat]
      sdssproc, fullname[0], flatimg, flativar, hdr=flathdr
;flatimg = flatimg[*,0:49]   ; TEST ???
;flativar = flativar[*,0:49] ; TEST ???
      dims = size(flatimg, /dimens)
      nx = dims[0]
      ny = dims[1]

      ;----------------------
      ; Create the array of preliminary flats
      if (iflat EQ 0) then pixflatarr = fltarr(nx,ny,nflat)

      ;----------------------
      ; Create spatial tracing from flat-field image

;      xcen = trace320crude(flatimg, yset=ycen, maxdev=0.15)
;      ntrace = (size(xcen, /dimens))[1]

;      xy2traceset, ycen, xcen, tset, ncoeff=5, maxdev=0.1
;      traceset2xy, tset, ycen, xsol

      ;----------------------
      ; Extract the flat-field image

      sigma = 1.0
      proftype = 1 ; Gaussian
      highrej = 10
      lowrej = 15
      nPoly = 4

      ; Determine YMODEL image
print, 'Begin extract ', fullname[0]
;      extract_image, flatimg, flativar, xsol, sigma, flat_flux, flat_fluxivar, $
;       proftype=proftype, wfixed=[1,1,1], $
;       highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1, $
;       ymodel=ymodel

      ; There are still systematics in the result from EXTRACT_IMAGE,
      ; so instead fit down each column of the image (which is slowly
      ; varying) by just doing a median filter.
      ymodel = 0.0 * flatimg
      for i=0, nx-1 do $
       ymodel[i,*] = median(transpose(flatimg[i,*]), 25)

      pixflatarr[*,*,iflat] = (flatimg > 1) / (ymodel > 1)

      ;----------------------
      ; Write FLATIMG and YMODEL to disk
      writefits, tmpname1[iflat], flatimg
flatimg = 0
      writefits, tmpname2[iflat], ymodel
ymodel = 0

   endfor

   ; Find deviant pixels in each pixflat
   meanimg = djs_avsigclip(pixflatarr, 3, sigrej=sigrej, maxiter=maxiter, $
    outmask=outmask)
meanimg = 0
pixflatarr = 0
   outmask = temporary(1-outmask) ; Change to 0=bad, 1=good

   flatimgsum = 0
   ymodelsum = 0
   for iflat=0, nflat-1 do begin

      flatimgsum = flatimgsum + outmask[*,*,iflat] * readfits(tmpname1[iflat])
      ymodelsum = ymodelsum + outmask[*,*,iflat] * readfits(tmpname2[iflat])

      rmfile, tmpname1[iflat]
      rmfile, tmpname2[iflat]

   endfor

   pixflat = (flatimgsum > 1) / (ymodelsum > 1)

   if (keyword_set(outfile)) then $
    writefits, outdir+outfile, pixflat

   return
end
;------------------------------------------------------------------------------
