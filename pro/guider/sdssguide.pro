;+
; NAME:
;   sdssguide
;
; PURPOSE:
;   Sloan telescope guider code
;
; CALLING SEQUENCE:
;   sdssguide, [ plugdir=, guidedir= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   plugdir    - Directory containing plPlugMapM files;
;                default to '/data/plugMapM'
;   guidedir   - Top-level directory containing guider images, with the
;                largest MJD subdirectory chosen;
;                default to '/data/gcam'
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The guide camera must already be in an independent loop taking images.
;   Pausing the guider does not pause whatever routine is independently
;     generating guider images.  When the guider is re-started, it will
;     begin by analyzing the most recent image, which may mean never
;     analyzing some older images.
;   Stars are fit to a symmetric Gaussian.
;
; EXAMPLES:
;
; BUGS:
;   Implement GUI interface with IDL widgets.
;   Send TCC offset commands.
;   If we can communicate with the TCC, then verify cartridge ID is mounted.
;   Decide upon gain for telescope offsets, and whether to apply rotation
;     and scale offsets (perhaps only once they become significantly bad).
;   Reject outlier stars in telescope offset fits.
;   Subtract dark frames, and mask bad pixels.
;   Should write the SOP-like guider output file.
;   Apply flat-fielding.
;   Add S/N ceiling to images.
;   Compute expected throughput, given the plugmap magnitudes.
;   Store history of all previous guider info, for use in any
;      more complex telescope offset decisions.
;   Line plot with history of offsets, FWHM, throughput.
;   Figure out how to change the ATV stretch levels dynamically, or use DS9.
;   Offset guiding option.
;   Need replacement script for spiralpattern when lost on the sky.
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/cartridgeInfo.par
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   24-Aug-2008  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
function sdssguide_guideobj, nfiber

   guideobj = replicate(create_struct( $
    'timestamp', 0L, $
    'plateid', 0, $
    'fiber', 0, $
    'enable', 0, $
    'goodfiber', 0, $
    'xfiber', 0., $
    'yfiber', 0., $
    'rfiber', 0., $
    'rot', 0., $
    'xstar', 0., $
    'ystar', 0., $
    'xstar_err', 0., $
    'ystar_err', 0., $
    'counts', 0., $
    'counts_err', 0., $
    'fwhm', 0., $
    'fwhm_err', 0., $
    'skyval', 0., $
    'skyval_err', 0., $
    'dx', 0., $
    'dy', 0., $
    'dra', 0., $
    'ddec', 0., $
    'dra_err', 0., $
    'ddec_err', 0. ), nfiber)

   return, guideobj
end

;------------------------------------------------------------------------------
; Generate an image of radial distance from a position xcen, ycen.

function sdssguide_radimg, nx, ny, xcen, ycen

   xaxis = djs_laxisgen([nx,ny], iaxis=0) - xcen
   yaxis = djs_laxisgen([nx,ny], iaxis=1) - ycen
   return, sqrt(xaxis^2 + yaxis^2)
end

;------------------------------------------------------------------------------
pro sdssguide_select_plugmap, pluginfo, plughdr, plugdir=plugdir
   common com_sdssguide_pluginfo, allplug

   if (NOT keyword_set(allplug)) then begin
      print, 'Locating most recent plPlugMapM files...'
      qcart = lonarr(17)
      files = file_search(djs_filepath('plPlugMapM-????-?????-??.par', $
       root_dir=plugdir))
      if (keyword_set(files)) then begin
         finfo = file_info(files)
         files = files[reverse(sort(finfo.mtime))]
         nfile = n_elements(files)
         ifile = 0L
         while (total(1-qcart) NE 0 AND ifile LT nfile) do begin
            yanny_read, files[ifile], hdr=plughdr
            ; Deal with cases where cartid is a string!?
            cartid = yanny_par(plughdr, 'cartridgeId')
            if (NOT keyword_set(cartid)) then cartid = 0
            if (qcart[(cartid-1)>0] EQ 0 AND cartid NE 0) then begin
               ; Find all pointings for this plugmap...
               file1 = fileandpath(files[ifile], path=path1)
               file1 = repstr(file1,'.par','*.par')
               filelist = file_search(djs_filepath(file1, root_dir=path1))
               for j=0L, n_elements(filelist)-1L do begin
                  pointing = strmid(fileandpath(filelist[j]), 24, 1)
                  pointing = repstr(pointing,'.','')
                  thisplug = {cartid: strtrim(cartid,2)+pointing, $
                   file: filelist[j]}
                  if (keyword_set(allplug)) then allplug = [allplug,thisplug] $
                   else allplug = thisplug
               endfor
               qcart[cartid-1] = 1
            endif
            ifile++
         endwhile
         if (keyword_set(allplug)) then allplug = allplug[sort(allplug.cartid)]
      endif
   endif

   if (keyword_set(allplug)) then begin
      print
      for i=0L, n_elements(allplug)-1 do $
       print, allplug[i].cartid, fileandpath(allplug[i].file), $
        format='(a4,x,a)'
      ct = 0
      while (ct EQ 0) do begin
         val = ''
         read, prompt='Select cartridge from first column above: ', val
         iselect = where(allplug.cartid EQ strupcase(strtrim(val,2)), ct)
      endwhile

      pluginfo = yanny_readone(allplug[iselect].file, hdr=plughdr, /anon)
      if (NOT keyword_set(pluginfo)) then begin
         splog, 'Error reading plugmap '+allplug[iselect].file
         return
      endif
      pluginfo = pluginfo[ where(pluginfo.holetype EQ 'GUIDE') ]
      pluginfo = pluginfo[ sort(pluginfo.fiberid) ]
   endif

   return
end

;------------------------------------------------------------------------------
function sdssguide_read_image, thisfile, hdr=imhdr

   if (fits_wait(thisfile, deltat=2, tmax=10) EQ 0) then return, 0
   image = mrdfits(thisfile, 0, imhdr, /fscale, /silent)

   ; Make the image look like the old camera binned 2x2
   dims = size(image,/dimens)
   case dims[0] of
      384: begin
         camera = 'old'
         camoffsets = [0., 0.]
         pixscale = 0.149 * 2
         rdnoise = 100. ; electrons
         gain = 10. ; electrons/ADU
         end
      768: begin
         camera = 'old'
         camoffsets = [0., 0.]
         pixscale = 0.149 * 2
         rdnoise = 100. ; electrons
         gain = 10. ; electrons/ADU
         image = rebin(image, dims/2) ; old camera unbinned
         end
      512: begin
         camera = 'new'
         camoffsets = [0., 0.]
         pixscale = 0.214 * 2
         rdnoise = 8.1 ; electrons
         gain = 1.4 ; electrons/ADU
         end
      524: begin ; changed to this dimension on Aug 17, 2010 = MJD 55425
         camera = 'new'
         camoffsets = [0., 0.]
         pixscale = 0.214 * 2
         rdnoise = 8.1 ; electrons
         gain = 1.4 ; electrons/ADU
         end
      1024: begin
         camera = 'new'
         camoffsets = [0., 0.]
         pixscale = 0.214 * 2
         rdnoise = 8.1 ; electrons
         gain = 1.4 ; electrons/ADU
         image = rebin(image, dims/2) ; new camera unbinned
         end
      1048: begin ; changed to this dimension on Aug 17, 2010 = MJD 55425
         camera = 'new'
         camoffsets = [0., 0.]
         pixscale = 0.214 * 2
         rdnoise = 8.1 ; electrons
         gain = 1.4 ; electrons/ADU
         image = rebin(image, dims/2) ; new camera unbinned
         end
   endcase

   biasval = median(image)
   image -= biasval ; bias-subtract
   var = (image>0) / gain + rdnoise^2 / gain^2
   invvar = 1. / var
   img = create_struct( $
    'filename', thisfile, $
    'camera', camera, $
    'camoffsets', camoffsets, $
    'pixscale', pixscale, $
    'biasval', biasval, $
    'imhdr', imhdr, $
    'image', image, $
    'invvar', invvar, $
    'immask', intarr(size(image,/dimens)) )

   return, img
end

;------------------------------------------------------------------------------
pro sdssguide_disable_fiber, cartinfo, guideobj

   read, prompt='Select fiber to disable: ', num
   i = where(guideobj.fiber EQ num, ct)
   if (ct GT 0) then begin
      cartinfo[i].exists = 'F'
      guideobj[i].enable = 0
   endif else begin
      print, 'Fiber not found'
   endelse

   return
end

;------------------------------------------------------------------------------
pro sdssguide_enable_fiber, cartinfo, guideobj

   read, prompt='Select fiber to enable: ', num
   i = where(guideobj.fiber EQ num, ct)
   if (ct GT 0) then begin
      cartinfo[i].exists = 'T'
      guideobj[i].enable = 1
   endif else begin
      print, 'Fiber not found'
   endelse

   return
end

;------------------------------------------------------------------------------
function sdssguide_is_valid, img
   if (n_elements(img.image) GT 1) then return, 1
   return, 0
end

;------------------------------------------------------------------------------
function sdssguide_cartinfo, cartid, img
   common com_sdssguide_cartinfo, allcart

   cartinfo = 0

   if (NOT keyword_set(allcart)) then begin
      cartfile = djs_filepath('cartridgeInfo.par', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdir='examples')
      allcart = yanny_readone(cartfile, /anon)
   endif

   if (keyword_set(allcart)) then begin
      i = where(allcart.cartridgeid EQ cartid)
      if (i[0] NE -1) then cartinfo = allcart[i] $
       else return, 0

      ; Convert to positions for this image's pixel scale and CCD offsets
;      cartinfo.xcen *= 0.149 / img.pixscale
;      cartinfo.ycen *= 0.149 / img.pixscale
;      cartinfo.radius *= 0.149 / img.pixscale
;      cartinfo.xcen += (0.149 / img.pixscale) * img.camoffsets[0]
;      cartinfo.ycen += (0.149 / img.pixscale) * img.camoffsets[1]
      if (cartinfo[0].cartridgeid LE 9) then begin
         cartinfo.xcen = 0.345 * allcart[i].ycen + 168
         cartinfo.ycen = 0.345 * allcart[i].xcen + 122
         cartinfo.radius = 0.345 * allcart[i].radius
      endif
   endif

   return, cartinfo
end

;------------------------------------------------------------------------------
; Parameters are as follows: XCEN, YCEN, FWHM, FLUX, SKYVAL
function sdssguide_psf_fn, pp, subimg=subimg, subivar=subivar

   dims = size(subimg, /dimens)
   modelimg = pp[3] * psf_gaussian(npixel=dims, fwhm=pp[2], centroid=pp[0:1], $
    /normalize) + pp[4]
   chivec = reform( (subimg - modelimg) * sqrt(subivar), n_elements(subimg))

   return, chivec
end

;------------------------------------------------------------------------------
function sdssguide_mask_fn, pp, subimg=subimg, radius=radius

   dims = size(subimg,/dimens)
   djs_photfrac, pp[0], pp[1], radius, $
    xdimen=dims[0], ydimen=dims[1], pixnum=pixnum, fracs=fracs
   mask = 0*subimg
   mask[pixnum] = fracs
   maxval = max(subimg)
   vdiff = reform(2 - subimg * mask / maxval, n_elements(subimg))

   return, sqrt(vdiff) ; This will maximize on simply the total counts
end

;------------------------------------------------------------------------------
; Centroid the circular fiber regions by maximizing the counts that
; fall within a circle.  Also decide if a fiber sees enough light
; above the bias level to be considered a good fiber.

pro sdssguide_centroid_fibers, cartinfo, img, guideobj, minval=minval

   nfiber = n_elements(cartinfo)

   for i=0L, nfiber-1L do begin
      xcen = cartinfo[i].xcen
      ycen = cartinfo[i].ycen
      radius = cartinfo[i].radius
      dims = size(img.image,/dimens)

      ; Define a sub-region within to work
      x1 = (floor(xcen - 2.0 * radius)) > 0
      x2 = (ceil(xcen + 2.0 * radius)) < (dims[0]-1)
      y1 = (floor(ycen - 2.0 * radius)) > 0
      y2 = (ceil(ycen + 2.0 * radius)) < (dims[1]-1)
      subimg = img.image[x1:x2,y1:y2]

      parinfo = {value: 0.D, limited: [0b,0b], limits: [0.d0,0.d0], $
       step: 0.05D, mpmaxstep: 1.0D}
      parinfo = replicate(parinfo, 2)
      parinfo[0].value = xcen - x1
      parinfo[0].limits = [0.5*radius,3.5*radius]
      parinfo[1].value = ycen - y1
      parinfo[1].limits = [0.5*radius,3.5*radius]

      functargs = {subimg: subimg, radius: radius}
      fitval = mpfit('sdssguide_mask_fn', parinfo=parinfo, $
       functargs=functargs, perror=perror, niter=niter, status=status, $
       /quiet)

      ; Leave these centroids unchanged if the fit failed
      if (n_elements(fitval) EQ 2) then begin
         guideobj[i].xfiber = fitval[0] + x1
         guideobj[i].yfiber = fitval[1] + y1
         ii = where(sdssguide_radimg(x2-x1+1,y2-y1+1,fitval[0],fitval[1]) $
          LT radius)
         medval = median(subimg[ii])
         guideobj[i].goodfiber = medval GT minval
      endif else begin
         guideobj[i].xfiber = cartinfo[i].xcen
         guideobj[i].yfiber = cartinfo[i].ycen
         guideobj[i].goodfiber = 0
      endelse
      guideobj[i].fiber = cartinfo[i].gprobeid
      guideobj[i].enable = cartinfo[i].exists EQ 'T'
      guideobj[i].rfiber = cartinfo[i].radius
      guideobj[i].rot = cartinfo[i].rot
   endfor

   return
end

;------------------------------------------------------------------------------
pro sdssguide_centroid_stars, img, guideobj, minsn=minsn

   nfiber = n_elements(guideobj)

   for i=0L, nfiber-1L do begin
      if (guideobj[i].enable AND guideobj[i].goodfiber GT 0) then begin
         xcen = guideobj[i].xfiber
         ycen = guideobj[i].yfiber
         radius = guideobj[i].rfiber
         dims = size(img.image,/dimens)

         ; Define a sub-region within to work
         x1 = (floor(xcen - 1.0 * radius)) > 0
         x2 = (ceil(xcen + 1.0 * radius)) < (dims[0]-1)
         y1 = (floor(ycen - 1.0 * radius)) > 0
         y2 = (ceil(ycen + 1.0 * radius)) < (dims[1]-1)
         subimg = img.image[x1:x2,y1:y2]
         radimg = sdssguide_radimg(x2-x1+1,y2-y1+1, xcen-x1, ycen-y1)
         subivar = img.invvar[x1:x2,y1:y2] * (radimg LT radius)
         qgood = subivar GT 0
         ii = where(qgood, ngood)
         ii = sort(subimg[ii])
         skyval = subimg[ii[ 0.25*ngood ]] ; select the 25th-percentile
         skyval = median(subimg[ii])

         ; Parameters are as follows: XCEN, YCEN, FWHM, FLUX, SKYVAL
         parinfo = {value: 0.D, limited: [0b,0b], limits: [0.d0,0.d0], $
          step: 0.D, mpmaxstep: 0.D}
         parinfo = replicate(parinfo, 5)

         ; Initial guesses for XCEN,YCEN will be flux-weighted means
         xvec = total((subimg - skyval) * qgood,2)
         yvec = total((subimg - skyval) * qgood,1)
         xthis = total(xvec * findgen(x2-x1+1)) / total(xvec)
         ythis = total(yvec * findgen(y2-y1+1)) / total(yvec)

         parinfo[0].value = xthis
         parinfo[0].limits = [-0.1*radius,+2.1*radius]
         parinfo[0].limited = [1b,1b]
         parinfo[0].step = 0.05
         parinfo[0].mpmaxstep = 1.0

         parinfo[1].value = ythis
         parinfo[1].limits = [-0.1*radius,+2.1*radius]
         parinfo[1].limited = [1b,1b]
         parinfo[1].step = 0.05
         parinfo[1].mpmaxstep = 1.0

         parinfo[2].value = 2.0/img.pixscale ; initial guess is 2 arcsec
         parinfo[2].limits = [0.75, 4.0]/img.pixscale ; seeing range in pixel
         parinfo[2].limited = [1b,1b]
         parinfo[2].step = 0.05
         parinfo[2].mpmaxstep = 1.0

         parinfo[3].value = total(subimg * qgood) - skyval * ngood

         parinfo[4].value = skyval

         functargs = {subimg: subimg, subivar: subivar}
         fitval = mpfit('sdssguide_psf_fn', parinfo=parinfo, $
          functargs=functargs, perror=perror, niter=niter, status=status, $
          /quiet)
      endif else begin
         fitval = 0
      endelse

      qgood = keyword_set(fitval)
      if (qgood) then qgood = status NE 0
      if (qgood) then qgood = fitval[3] GT minsn*perror[3] $
       AND perror[3] GT 0 ; AND perror[2] GT 0
      if (qgood) then begin
         guideobj[i].xstar = x1 + fitval[0]
         guideobj[i].ystar = y1 + fitval[1]
         guideobj[i].xstar_err = perror[0]
         guideobj[i].ystar_err = perror[1]
         guideobj[i].fwhm = fitval[2]
         guideobj[i].fwhm_err = perror[2]
         guideobj[i].counts = fitval[3]
         guideobj[i].counts_err = perror[3]
         guideobj[i].skyval = fitval[4]
         guideobj[i].skyval_err = perror[4]
      endif else begin
         guideobj[i].xstar = 0
         guideobj[i].ystar = 0
         guideobj[i].xstar_err = 0
         guideobj[i].ystar_err = 0
         guideobj[i].fwhm = 0
         guideobj[i].fwhm_err = 0
         guideobj[i].counts = 0
         guideobj[i].counts_err = 0
         guideobj[i].skyval = 0
         guideobj[i].skyval_err = 0
         guideobj[i].goodfiber = -1
      endelse
   endfor

   return
end

;------------------------------------------------------------------------------
pro sdssguide_compute_offsets, guideobj, img

   pixscale = img.pixscale
   sinphi = sin(guideobj.rot / !radeg)
   cosphi = cos(guideobj.rot / !radeg)
   dx = guideobj.xstar - guideobj.xfiber
   dy = guideobj.ystar - guideobj.yfiber
   guideobj.dra = pixscale * (dx * cosphi + dy * sinphi)
   guideobj.ddec = pixscale * (dx * cosphi + dy * sinphi)
   guideobj.dra_err = pixscale $
    * sqrt((guideobj.xstar_err * cosphi)^2 + (guideobj.ystar_err * sinphi)^2)
   guideobj.ddec_err = pixscale $
    * sqrt((guideobj.xstar_err * cosphi)^2 + (guideobj.ystar_err * sinphi)^2)

   return
end

;------------------------------------------------------------------------------
pro sdssguide_compute_teloffset, plughdr, img, guideobj, guideiter, pluginfo

   asecpermm = 16.5245

   guideiter = create_struct( $
    'filename', '', $
    'timestamp', 0L, $
    'exptime', 0, $
    'nguide', 0, $
    'errors', fltarr(4), $
    'offsets', fltarr(4), $
    'lst', 0., $
    'ra', 0., $
    'dec', 0.)
   guideiter.filename = img.filename
   guideiter.timestamp = 24.0D * 3600.0D * current_mjd()
   guideiter.exptime = sxpar(img.imhdr, 'EXPTIME')
   guideiter.nguide = n_elements(guideobj)
   guideiter.ra = yanny_par(plughdr, 'raCen')
   guideiter.dec = yanny_par(plughdr, 'decCen')

   nfiber = n_elements(guideobj)
   bvec = [guideobj.dra, guideobj.ddec]
   berr = [guideobj.dra_err, guideobj.ddec_err] * (guideobj.goodfiber GT 0)
   sqivar = (berr GT 0) / (berr + (berr LE 0))
   theta = atan(pluginfo.yfocal, pluginfo.xfocal)
   sint = sin(theta)
   cost = cos(theta)
   amatrix = fltarr(2*nfiber,4)
   amatrix[*,0] = [fltarr(nfiber)+1, fltarr(nfiber)+0]
   amatrix[*,1] = [fltarr(nfiber)+0, fltarr(nfiber)+1]
   amatrix[*,2] = [-pluginfo.xfocal*asecpermm, pluginfo.yfocal*asecpermm]
   amatrix[*,3] = $
    [-pluginfo.xfocal*asecpermm*cost, pluginfo.yfocal*asecpermm*sint]
   chi2 = computechi2(bvec, sqivar, amatrix, acoeff=acoeff, var=var, dof=dof)
   print, 'Chi^2 of fit = ', chi2, ' for ', long(dof), ' DOF'
   print, 'Computed RA offset = ', acoeff[0], ' arsec'
   print, 'Computed Dec offset = ', acoeff[1], ' arcsec'
   print, 'Computed rotator offset = ', acoeff[2]*!radeg*3600., ' arcsec'
   print, 'Computed scale offset = ', acoeff[3]*100., ' percent'

   guideiter.errors = acoeff

   return
end

;------------------------------------------------------------------------------
pro sdssguide_display, img, guideobj

   skyval = median(guideobj.skyval)
   atv, img.image, min=0.5*skyval-200, max=3.0*(skyval+200)

   nfiber = n_elements(guideobj)

   ; Draw a circle around each good fiber (but not around bad ones)
   for i=0L, nfiber-1L do begin
      if (guideobj[i].enable EQ 0) then text = 'DISABLE' $
       else if (guideobj[i].goodfiber EQ 0) then text = 'NO SKY' $
       else if (guideobj[i].goodfiber EQ -1) then text = 'FAINT' $
       else text = string(guideobj[i].fwhm*img.pixscale,format='(f4.2)') $
        +string(34b)
      rplot = guideobj[i].rfiber
      npt = 2*!pi*rplot
      xplot = cos(2*!pi*findgen(npt+1)/npt) * rplot + guideobj[i].xfiber
      yplot = sin(2*!pi*findgen(npt+1)/npt) * rplot + guideobj[i].yfiber
      atvplot, xplot, yplot, color='red'
      atvxyouts, guideobj[i].xfiber, $
       guideobj[i].yfiber+guideobj[i].rfiber, $
       strtrim(guideobj[i].fiber,2)+' ', color='red', charsize=2, align=1
      atvxyouts, guideobj[i].xfiber, $
       guideobj[i].yfiber-guideobj[i].rfiber-16, $
       text, color='green', charsize=1.5, align=0.5

      if (guideobj[i].xstar_err GT 0 AND guideobj[i].ystar_err GT 0) then $
       atvplot, guideobj[i].xstar, guideobj[i].ystar, psym=1, color='red', $
        symsize=2
   endfor

   dims = size(img.image,/dimens)
   atvxyouts, dims[0]/2, dims[1]/2-220, $
    'Cart #'+strtrim(sxpar(img.imhdr,'FLATCART'),2) $
    +'    '+strmid(fileandpath(img.filename),0,9) $
    +'    '+strtrim(sxpar(img.imhdr,'DATE-OBS'),2), $
    charsize=2, color='red', align=0.5

   return
end

;------------------------------------------------------------------------------
pro sdssguide, plugdir=plugdir1, guidedir=guidedir1

;plugdir1='/Users/schlegel/guider/plugmap' ; ???
;guidedir1 ='/Users/schlegel/rawdata/51820/guider' ; ???
;guidedir1 ='/Users/schlegel/rawdata/gcam/55123' ; ???
;guidedir1 ='/Users/schlegel/rawdata/gcam' ; ???

   ; Set defaults
   if (keyword_set(plugdir1)) then plugdir = plugdir1 $
    else plugdir = '/data/plugMapM'
   if (keyword_set(guidedir1)) then guidedir = guidedir1 $
    else guidedir = '/data/gcam'
   minval = 0. ; Minimum median value within fiber
   minsn = 10. ; Minimum star S/N
   quiet = !quiet
   !quiet = 1
   except = !except
   !except = 0

   ; Guider loop
   lastfile = ''
   cc = 'M'
   while (cc NE 'Q') do begin
;      if (cc EQ 'M') then begin
;         sdssguide_select_plugmap, pluginfo, plughdr, plugdir=plugdir
;         cartid = yanny_par(plughdr, 'cartridgeId')
;         cartinfo = 0
;      endif

      if (cc EQ 'D') then sdssguide_disable_fiber, cartinfo, guideobj
      if (cc EQ 'E') then sdssguide_enable_fiber, cartinfo, guideobj
      if (cc EQ 'P') then begin
         junk = ''
         read, 'Paused -- Press return to continue: ', junk
      endif

      mjd_dir = file_search(djs_filepath('5[0-9][0-9][0-9][0-9]', $
       root_dir=guidedir))
      if (keyword_set(mjd_dir) EQ 0) then begin
         splog, 'No MJD directories found in GUIDEDIR='+guidedir
         return
      endif
      junk = max(long(fileandpath(mjd_dir)), imax)
      mjd_dir = mjd_dir[imax]
      files = file_search(djs_filepath('g*.fits*', root_dir=mjd_dir))
      finfo = file_info(files)
      junk = max(finfo.mtime, i) ; Find the most recently modified file
      thisfile = files[i]
;thisfile ='/Users/schlegel/rawdata/51820/guider/gimg0052.fits.gz' ; ???
;thisfile ='/Users/schlegel/rawdata/51820/guider/gimg0127.fits.gz' ; ???
;thisfile ='/Users/schlegel/rawdata/gcam/55123/gimg-1000.fits' ; ???
      if (thisfile NE lastfile) then begin
         lastfile = thisfile
         print, ' '
         print, 'File = '+thisfile
         img = sdssguide_read_image(thisfile, hdr=imhdr)
         if (sdssguide_is_valid(img)) then begin
            ; Read the cartridge info with the first good image
            cartid = sxpar(imhdr, 'FLATCART')
;if (NOT keyword_set(cartid)) then cartid = 10 ; ???
            if (NOT keyword_set(cartinfo)) then $
             cartinfo = sdssguide_cartinfo(cartid, img)
            if (NOT keyword_set(cartinfo)) then begin
               splog, 'No cartridge info found for cartridge ID ', cartid
               return
            end

            guideobj = sdssguide_guideobj(n_elements(cartinfo))
            sdssguide_centroid_fibers, cartinfo, img, guideobj, minval=minval
            sdssguide_centroid_stars, img, guideobj, minsn=minsn
            sdssguide_compute_offsets, guideobj, img
;            sdssguide_compute_teloffset, plughdr, img, guideobj, $
;             guideiter, pluginfo
            sdssguide_display, img, guideobj
         endif else begin
            wait, 1 ; Wait 1 sec if no new image found yet
         endelse
      endif
;lastfile = '' ; ???

      print
      print, 'Options:'
      print, '         d - Disable fiber'
      print, '         e - Enable fiber'
;      print, '         m - Load plugmap'
      print, '         p - Pause'
      print, '         q - Quit'

      cc = strupcase(get_kbrd(0))
   endwhile

   !quiet = quiet
   !except = except

   return
end
;------------------------------------------------------------------------------
