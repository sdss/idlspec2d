;+
; NAME:
;   spflatten
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   spflatten, flatarr
;
; INPUTS:
;   flatname   - Name of flat-field SDSS image(s)
;
; OPTIONAL KEYWORDS:
;   indir      - Input directory for FLATNAME; default to './'
;   tmpdir     - Directory for temporary files; default to './'
;
; OUTPUTS:
;   pixflat    - Image containing all the information about pixel-to-pixel
;                variations.  Illumination variations are removed.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Need to remove cosmics!!! ???
;   Look for very devient frames and exclude them.
;
; PROCEDURES CALLED:
;   readfits()
;   writefits
;
; REVISION HISTORY:
;   13-Oct-1999  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------

function spflatten, flatname, indir=indir

   if (NOT keyword_set(indir)) then indir = './'
   if (NOT keyword_set(tmpdir)) then tmpdir = './'
   nflat = N_elements(flatname)
   ngrow = 2

   for iflat=0, nflat-1 do begin

      ;----------------------
      ; Read flat-field image

      fullpath = filepath(flatname[iflat], root_dir=indir)
      fullname = findfile(fullpath, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find flat image' + flatname[iflat]
      sdssproc, fullname[0], flatimg, flativar, hdr=flathdr
      dims = size(flatimg, /dimens)
      nx = dims[0]
      ny = dims[1]

      ;----------------------
      ; Create the array of preliminary flats
      if (iflat EQ 0) then pixflatarr = fltarr(nx,ny,nflat)

      ;----------------------
      ; Create spatial tracing from flat-field image

      xcen = trace320crude(flatimg, yset=ycen, maxdev=0.15)
      ntrace = (size(xcen, /dimens))[1]

      xy2traceset, ycen, xcen, tset, ncoeff=5, maxdev=0.1
      traceset2xy, tset, ycen, xsol

      ;----------------------
      ; Extract the flat-field image

      sigma = 1.0
      proftype = 1 ; Gaussian
      highrej = 10
      lowrej = 15
      nPoly = 4

      ; Extract, setting mask=0 for bad pixels, =1 for good
      mask = 0
      extract_image, flatimg, flativar, xsol, sigma, flat_flux, flat_fluxivar, $
       proftype=proftype, wfixed=[1,1,1], $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1, $
       ymodel=ymodel, mask=mask

      pixflatarr[*,*,iflat] = (flatimg > 1) / (ymodel > 1)

      ;----------------------
      ; Write FLATIMG and YMODEL to disk
      writefits, tmpdir+'/tmp.flatimg.'+strtrim(string(iflat),2)+'.fits', $
       flatimg
flatimg = 0
      writefits, tmpdir+'/tmp.ymodel.'+strtrim(string(iflat),2)+'.fits', $
       ymodel
ymodel = 0

   endfor

   ; Find deviant pixels in each pixflat
   sigrej = 3.0
   maxiter = 2
   meanimg = djs_avsigclip(pixflatarr, 3, sigrej=sigrej, maxiter=maxiter, $
    outmask=outmask)
;meanimg = 0
;pixflatarr = 0
   outmask = temporary(1-outmask) ; Change to 0=bad, 1=good

   flatimgsum = 0
   ymodelsum = 0
   for iflat=0, nflat-1 do begin

      flatimgsum = flatimgsum + outmask[*,*,iflat] * $
       readfits(tmpdir+'/tmp.flatimg.'+strtrim(string(iflat),2)+'.fits')
      ymodelsum = ymodelsum + outmask[*,*,iflat] * $
       readfits(tmpdir+'/tmp.ymodel.'+strtrim(string(iflat),2)+'.fits')

   endfor

   pixflat = (flatimgsum > 1) / (ymodelsum > 1)
stop

   return, pixflat
end
;------------------------------------------------------------------------------
