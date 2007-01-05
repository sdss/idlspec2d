;+
; NAME:
;   spframe_read
;
; PURPOSE:
;   Read data from either an spFrame or spCFrame file.
;
; CALLING SEQUENCE:
;   spframe_read, filename, [ indx, objflux=, objivar=, mask=, $
;    wset=, loglam=, dispset=, dispimg=, ximg=, $
;    plugmap=, skyflux=, superflat=, hdr=, adderr= ]
;
; INPUTS:
;   filename   - Input file name
;
; OPTIONAL INPUTS:
;   indx       - Optional 0-indexed row numbers; default to all
;   adderr     - Additional error to add to the formal errors, as a
;                fraction of the flux; default to none
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   objflux    - Object flux
;   objivar    - Object inverse variance (units of 1/OBJFLUX^2)
;   mask       - Pixel bit mask
;   wset       - Trace-set for wavelength solution
;   loglam     - Wavelength image (vacuum log-10 Ang)
;   dispset    - Trace-set for dispersion solution
;   dispimg    - Dispersion image (per native pixel)
;   ximg       - X position on CCD image
;   skyflux    - Sky flux (same units as OBJFLUX)
;   superflat  - Superflat vector from quartz lamps
;   hdr        - FITS header for HDU#0
;
; COMMENTS:
;   The spFrame and spCFrame files contain nearly identical HDUs,
;   except for HDU #3 and #4 which are trace-sets in spFrame, and
;   more easily-interpreted 2-dimensional images in spCFrame:
;
;            spFrame   spCFrame
;   HDU #0:  Flux      Flux
;   HDU #1:  Invvar    Invvar
;   HDU #2:  mask      mask
;   HDU #3:  wset      loglam
;   HDU #4:  dispset   dispimg
;   HDU #5:  plugmap   plugmap
;   HDU #6:  sky       sky
;   HDU #7:  ximg      ximg
;   HDU #8:  superflat
;   HDU #9:  skystruct
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   fileandpath()
;   headfits()
;   lookforgzip()
;   mrdfits()
;   traceset2xy
;   traceset_trim()
;
; REVISION HISTORY:
;   05-Feb-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro spframe_read, filename, indx, objflux=objflux, objivar=objivar, $
 mask=mask, wset=wset, loglam=loglam, dispset=dispset, dispimg=dispimg, $
 ximg=ximg, plugmap=plugmap, skyflux=skyflux, superflat=superflat, $
 hdr=hdr, adderr=adderr

   qtrim = n_elements(indx) GT 0

   if (NOT keyword_set(filename)) then $
    message, 'Must specify FILENAME'

   ; Is this a flux-caliibrated frame file?
   qcframe = strmatch(fileandpath(filename),'spCFrame*')
   if (qcframe AND (arg_present(wset) OR arg_present(dispset) $
    OR arg_present(superflat))) then $
    message, 'Cannot request WSET or DISPSET from a spCFrame file'

   thisfile = lookforgzip(filename[0])

   if (arg_present(hdr)) then hdr = headfits(thisfile[0])

   if (arg_present(objflux) $
    OR (arg_present(objivar) AND keyword_set(adderr))) then begin
      objflux = mrdfits(thisfile[0], 0, /silent)
      if (qtrim) then objflux = objflux[*,indx]
   endif

   if (arg_present(objivar)) then begin
      objivar = mrdfits(thisfile[0], 1, /silent)
      if (qtrim) then objivar = objivar[*,indx]
      if (keyword_set(adderr)) then begin
         gmask = objivar NE 0 ; =1 for good points
         objivar = 1.0 / ( 1.0/(objivar + (1-gmask)) $
          + (adderr * (objflux>0))^2 ) * gmask
      endif
   endif

   if (arg_present(mask)) then begin
      mask = mrdfits(thisfile[0], 2, /silent)
      if (qtrim) then mask = mask[*,indx]
   endif

   if (qcframe) then begin
      if (arg_present(loglam)) then begin
         loglam = mrdfits(thisfile[0], 3, /silent)
         if (qtrim) then loglam = loglam[*,indx]
      endif
   endif else begin
      if (arg_present(wset) OR arg_present(loglam)) then begin
         wset = mrdfits(thisfile[0], 3, /silent)
         if (qtrim) then wset = traceset_trim(wset, indx)
         if (arg_present(loglam) AND keyword_set(wset[0])) then $
          traceset2xy, wset, xtmp, loglam
         xtmp = 0
      endif
   endelse

   if (qcframe) then begin
      if (arg_present(dispimg)) then begin
         dispimg = mrdfits(thisfile[0], 4, /silent)
         if (qtrim) then dispimg = dispimg[*,indx]
      endif
   endif else begin
      if (arg_present(dispset) OR arg_present(dispimg)) then begin
         dispset = mrdfits(thisfile[0], 4, /silent)
         if (qtrim) then dispset = traceset_trim(dispset, indx)
         if (arg_present(dispimg)) then traceset2xy, dispset, xtmp, dispimg
         xtmp = 0
      endif
   endelse

   if (arg_present(plugmap)) then begin
      plugmap = mrdfits(thisfile[0], 5, /silent)
      if (qtrim) then plugmap = plugmap[indx]
   endif

   if (arg_present(skyflux)) then begin
      skyflux = mrdfits(thisfile[0], 6, /silent)
      if (qtrim) then skyflux = skyflux[*,indx]
   endif

   if (arg_present(ximg)) then begin
      ximg = mrdfits(thisfile[0], 7, /silent)
      if (qtrim) then ximg = ximg[*,indx]
   endif

   if (arg_present(superflat)) then begin
      superflat = mrdfits(thisfile[0], 8, /silent)
      if (qtrim) then superflat = superflat[*,indx]
   endif

   return
end
;------------------------------------------------------------------------------
