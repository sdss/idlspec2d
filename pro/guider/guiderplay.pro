;+
; NAME:
;   guiderplay
;
; PURPOSE:
;   Replay guider images from a night
;
; CALLING SEQUENCE:
;   guidermovie, [ mjd=, wtime= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - MJD subdirectory; if not set, then use largest number
;                in $GUIDE_DIR directory
;   wtime      - Wait time between images; default to 0 sec
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
;
; REVISION HISTORY:
;   01-May-2010  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
forward_function sdssguide_read_image, sdssguide_is_valid, $
 sdssguide_cartinfo, sdssguide_guideobj
;------------------------------------------------------------------------------
pro guiderplay, mjd=mjd1, wtime=wtime1

   sdssguide, guidedir='test' ; Compile these routines

   ; Set defaults
   guidedir = getenv('GUIDE_DIR')
   if (NOT keyword_set(guidedir)) then guidedir = '/data/gcam'
   if (keyword_set(mjd1)) then begin
      mjd = get_mjd_dir(guidedir, mjd=long(mjd1[0]))
   endif else begin
      mjd = get_mjd_dir(guidedir,mjstart=1,mjend=99999)
      if (keyword_set(mjd)) then mjd = mjd[(reverse(sort(mjd)))[0]]
   endelse
   if (NOT keyword_set(mjd)) then begin
      splog, 'No MJD directory found in '+guidedir
      return
   endif
   if (n_elements(wtime1) EQ 0) then wtime = 0. $
    else wtime = wtime1[0]

   minval = 0. ; Minimum median value within fiber
   minsn = 10. ; Minimum star S/N
   quiet = !quiet
   !quiet = 1
   except = !except
   !except = 0

   allfile = findfile(filepath('gimg-????.fits*', root_dir=guidedir, $
    subdir=strtrim(mjd,2)), count=nfile)
   if (nfile EQ 0) then begin
      splog, 'No gimg files found in '+guidedir+'/'+strtrim(mjd,2)
      return
   endif

   for iexp=0L, nfile-1L do begin

      thisfile = allfile[iexp++]
      print, 'File = '+thisfile
      img = sdssguide_read_image(thisfile, hdr=imhdr)
      if (sdssguide_is_valid(img)) then begin
         ; Read the cartridge info with the first good image
         cartid = sxpar(imhdr, 'FLATCART')
         if (NOT keyword_set(cartid)) then cartid = 10 ; ???
         cartinfo = sdssguide_cartinfo(cartid, img)
         if (NOT keyword_set(cartinfo)) then begin
            splog, 'No cartridge info found for cartridge ID ', cartid
            return
         end

         guideobj = sdssguide_guideobj(n_elements(cartinfo))
         sdssguide_centroid_fibers, cartinfo, img, guideobj, minval=minval
         sdssguide_centroid_stars, img, guideobj, minsn=minsn
         sdssguide_compute_offsets, guideobj, img
;         sdssguide_compute_teloffset, plughdr, img, guideobj, $
;          guideiter, pluginfo
         sdssguide_display, img, guideobj
         wait, wtime ; Wait 1 sec if no new image found yet
      endif
   endfor

   !quiet = quiet
   !except = except

   return
end
;------------------------------------------------------------------------------
