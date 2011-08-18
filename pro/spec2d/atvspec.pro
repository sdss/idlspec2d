
;+
; NAME:
;   atvspec
;
; PURPOSE:
;   Display cutouts of all raw images around specified wavelengths
;
; CALLING SEQUENCE:
;   atvspec, plate, fiberid, wave=, [ mjd=, xsize=, ysize=, $
;    psym=, symsize=, color=, _EXTRA=, /verbose ]
;
; INPUTS:
;   plate      - Plate (scalar)
;   fiberid    - Fiber ID (scalar)
;   wave       - Wavelength(s) at which to display raw images (vacuum Ang)
;
; OPTIONAL INPUTS:
;   mjd        - MJD (scalar)
;   xsize      - X dimension of cutout of raw images; default to 50 pix
;   ysize      - Y dimension of cutout of raw images; default to 50 pix
;   psym       - PSYM keyword for ATVPLOT; default to 6 (squares)
;   symsize    - SYMSIZE keyword for ATVPLOT; default to 3
;   color      - COLOR keyword for ATVPLOT; default to 'red'
;   verbose    - If set, then print X,Y locations on CCD
;   _EXTRA     - Keywords for ATV, such as /ALIGN or /STRETCH
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Display cutouts from the raw images, with each row of cutouts
;   corresponding to each WAVE and each column corresponding to the
;   exposure numbers.  Assume wavelengths < 6000 Ang should display
;   the blue images, and > 6000 Ang should dislpay the red images.
;
; EXAMPLES:
;   Display raw images around two bright sky lines
;     IDL> atvspec, 3686, 1, mjd=55268, wave=[5578.8876,6302.0466]
;
; BUGS:
;
; REVISION HISTORY:
;   08-Jul-2010  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro atvspec, plate, fiberid, mjd=mjd, wave=wavecen, $
 xsize=xsize1, ysize=ysize1, psym=psym1, symsize=symsize1, color=color1, $
 verbose=verbose, _EXTRA=KeywordsForATV

   if (n_params() LT 2 OR keyword_set(wavecen) EQ 0) then begin
      print, 'Must set PLATE, FIBERID, WAVECEN'
      return
   endif

   ; Set defaults
   nwave = n_elements(wavecen)
   if (keyword_set(xsize1)) then xsize = xsize1 else xsize = 50L
   if (keyword_set(ysize1)) then ysize = ysize1 else ysize = 50L
   if (keyword_set(psym1)) then psym = psym1 else psym = 6
   if (keyword_set(symsize1)) then symsize = symsize1 else symsize = 3
   if (keyword_set(color1)) then color = color1 else color = 'red'

   readonespec, plate, fiberid, mjd=mjd, wave=wave, ximg=ximg, $
    framehdr=framehdr, expnum=expnum
   nfile = n_elements(expnum)
   if (nfile EQ 0) then begin
      print, 'No files found'
      return
   end
   cameras = strarr(nfile)
   for i=0, nfile-1 do cameras[i] = strtrim(sxpar(*framehdr[i], 'CAMERAS'),2)
   mjds = strarr(nfile)
   for i=0, nfile-1 do mjds[i] = strtrim(sxpar(*framehdr[i], 'MJD'),2)
   filename = strarr(nfile)
   for i=0, nfile-1 do begin
      filename[i] = 'sdR-'+cameras[i] $
       +'-'+string(expnum[i],format='(i8.8)')+'.fit*'
      filename[i] = filepath(filename[i], root_dir=getenv('BOSS_SPECTRO_DATA'), $
       subdir=mjds[i])
      filename[i] = findfile(filename[i], count=ct)
   endfor

   explist = expnum[uniq(expnum,sort(expnum))]
   nexp = n_elements(explist)

   npix = (size(wave,/dimens))[0]
   xvec = lindgen(npix)

   bigimg = fltarr(nexp*xsize,nwave*ysize)

   xlist = fltarr(nexp,nwave)
   ylist = fltarr(nexp,nwave)
   lastfile = ''
   for iexp=0, nexp-1 do begin
      for iwave=0, nwave-1 do begin
         thiscam = wavecen[iwave] LT 6000 ? 'b' : 'r'
         j = (where(strmid(cameras,0,1) EQ thiscam $
          AND expnum EQ explist[iexp], ct))[0]
         if (ct GT 0 AND keyword_set(filename[j>0])) then begin
            if (filename[j] NE lastfile) then $
             sdssproc, filename[j], thisimg, /silent $
            else $
             lastfile = filename[j]
            thisdim = size(thisimg,/dimens)

            ; Find the x,y location on this particular image
            linterp, wave[*,j], xvec, wavecen[iwave], ycen
            linterp, xvec, ximg[*,j], ycen, xcen
            if (keyword_set(verbose)) then $
             print, fileandpath(filename[j]), ' X=', xcen, ' Y=', ycen
            x0 = ((round(xcen) - xsize/2) > 0) < (thisdim[0]-xsize)
            x1 = x0 + xsize-1
            y0 = ((round(ycen) - ysize/2) > 0) < (thisdim[1]-ysize)
            y1 = y0 + ysize-1
            bigimg[iexp*xsize:(iexp+1)*xsize-1,iwave*ysize:(iwave+1)*ysize-1] $
             = thisimg[x0:x1,y0:y1]
            xlist[iexp,iwave] = xcen - x0 + iexp*xsize
            ylist[iexp,iwave] = ycen - y0 + iwave*ysize
         endif
      endfor
   endfor

   atv, bigimg, _EXTRA=KeywordsForATV
   atvplot, xlist, ylist, psym=psym, symsize=symsize, color=color
   for iexp=1, nexp-1 do $
    atvplot, iexp*xsize+0.5+[0,0], [-0.5,nwave*ysize-0.5], color=color
   for iwave=1, nwave-1 do $
    atvplot, [-0.5,nexp*xsize-0.5], iwave*ysize+0.5+[0,0], color=color

   return
end
;------------------------------------------------------------------------------
