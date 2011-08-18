
;+
; NAME:
;   daytime_test
;
; PURPOSE:
;   Look for spectro exposures taken during the day and not marked as test.
;
; CALLING SEQUENCE:
;   daytime_test
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Find all the files matching '$BOSS_SPECTRO_DATA/$MJD/sdR-b1-????????.fit*',
;   and print a warning message for any science, flat, or arc exposures
;   that are not marked as test exposures, and were taken during the
;   daytime.  The warnings are also printed to the file 'daytime.log'.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   get_mjd_dir()
;   get_tai
;   headfits()
;   sdsshead()
;   splog
;   sunpos
;   sxpar()
;   zenpos
;
; INTERNAL SUPPORT PROCEDURES:
;   daytime_test1
;
; REVISION HISTORY:
;   27-Apr-2003  Written by David Schlegel, Princeton (not checked in then)
;-
;------------------------------------------------------------------------------
pro daytime_test1, filename

   COMMON SITE, lat, lng, tzone
   lat = 32.780361d0
   lng = 360. - 105.820417d0
   tzone = 7 ; ???

   DRADEG = 180.d0/!dpi

   hdr = sdsshead(filename)
   if (strtrim(sxpar(hdr,'QUALITY'),2)  NE 'excellent') then return
   flavor = strtrim(sxpar(hdr,'FLAVOR'),2)
   if (flavor NE 'science' AND flavor NE 'flat' AND flavor NE 'arc') $
    then return

   get_tai, hdr, tai_beg, tai_mid, tai_end
   jd = 2400000.5D + tai_mid / (24.D*3600.D)
   sunpos, jd, ra1, dec1 ; returns degrees
   zenpos, jd, ra2, dec2 ; returns radians
   ra2 = ra2 * DRADEG
   dec2 = dec2 * DRADEG
   adiff = djs_diff_angle(ra1,dec1,ra2,dec2)

   if (adiff LT 90) then begin
     infstring  = string(sxpar(hdr,'PLATEID'), sxpar(hdr,'MJD'), $
           sxpar(hdr,'EXPOSURE'), format='(i4," ",i5," ",i8.8)')
      splog, 'Warning: Sun above the horizon ' + string(90.-adiff) + ' deg ' + $
       ' Plate=' + infstring + ' ' + flavor
   endif

end
;------------------------------------------------------------------------------
pro daytime_test

   splog, filename='daytime.log'

   ; Get the list of MJDs
   rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   mjdlist = get_mjd_dir(rawdata_dir, mjstart=1, mjend=99999, mjd='?????')

   for imjd=0L, n_elements(mjdlist)-1 do begin
      files = findfile(filepath('sdR-b1-????????.fit*', $
       root_dir=rawdata_dir, subdir=mjdlist[imjd]), count=nfile)
      splog, 'Working on MJD=', mjdlist[imjd], ' nfile=', nfile, prelog=''
      for ifile=0, nfile-1 do begin
         print, format='("File ",i5," of ",i5,a1,$)', $
          ifile+1, nfile, string(13b)
         daytime_test1, files[ifile]
      endfor
   endfor

   splog, /close

   return
end
;------------------------------------------------------------------------------
