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
;
; OUTPUTS:
;   pixflat    - Image containing all the information about pixel-to-pixel
;                variations.  Illumination variations are removed.
;
; OPTIONAL OUTPUTS:
;   xsol       - From last image
;   flat_flux  - From last image
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Need to remove cosmics!!! ???
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   13-Oct-1999  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------

function spflatten, flatname, indir=indir, xsol=xsol, flat_flux=flat_flux

   if (NOT keyword_set(indir)) then indir = './'
   nflat = N_elements(flatname)

   flatimgsum = 0
   ymodelsum = 0

   for iflat=0, nflat-1 do begin

      ;----------------------
      ; Read flat-field image

      fullpath = filepath(flatname[iflat], root_dir=indir)
      fullname = findfile(fullpath, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find flat image' + flatname[iflat]
      sdssproc, fullname[0], flatimg, flativar, hdr=flathdr

      ;----------------------
      ; Create spatial tracing from flat-field image

      xcen = trace320crude(flatimg, yset=ycen, maxdev=0.15)

      xy2traceset, ycen, xcen, tset, ncoeff=5, maxdev=0.1
      traceset2xy, tset, ycen, xsol

      ;----------------------
      ; Extract the flat-field image

      sigma = 1.0
      proftype = 1 ; Gaussian
      highrej = 20
      lowrej = 25
      nPoly = 4

      extract_image, flatimg, flativar, xsol, sigma, flat_flux, flat_fluxivar, $
       proftype=proftype, wfixed=[1,1,1], $
       highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1, ymodel=ymodel

      flatimgsum = flatimgsum + flatimg
      ymodelsum = ymodelsum + ymodel

   endfor

   pixflat = (flatimgsum > 1) / (ymodelsum > 1)

   return, pixflat
end
;------------------------------------------------------------------------------
