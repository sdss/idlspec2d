;+
; NAME:
;   warn_daytime
;
; PURPOSE:
;   Print a warning if an excellent-quality exposure taken with the Sun up.
;
; CALLING SEQUENCE:
;   warn_daytime, hdr
;
; INPUTS:
;   hdr        - FITS header from SDSS spectroscopic exposure
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The only output from this procedure is to print (using SPLOG) a
;   warning message if an "excellent"-quality exposure for a flat, arc,
;   or science was taken with the Sun above the horizon.  This should
;   basically never happen -- those should almost always be test data.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   get_tai
;   splog
;   sunpos
;   sxpar()
;   zenpos
;
; REVISION HISTORY:
;   28-Apr-2003  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro warn_daytime, hdr

   COMMON SITE, lat, lng, tzone
   lat = 32.780361d0
   lng = 360. - 105.820417d0
   tzone = 7

   DRADEG = 180.d0/!dpi

   ;----------
   ; Don't trigger any warnings if not an "excellent"-quality exposure,
   ; and if not a flat, arc, or science exposure.

   if (strtrim(sxpar(hdr,'QUALITY'),2)  NE 'excellent') then return
   flavor = strtrim(sxpar(hdr,'FLAVOR'),2)
   if (flavor NE 'science' AND flavor NE 'flat' AND flavor NE 'arc') $
    then return

   ;----------
   ; Compute the mid-point of the exposure

   get_tai, hdr, tai_beg, tai_mid, tai_end, /silent

   ;----------
   ; Compute where the Sun is relative to the horizon.
   ; Specifically, ADIFF is the angle of the Sun *above* the horizon.

   jd = 2400000.5D + tai_mid / (24.D*3600.D)
   sunpos, jd, ra1, dec1 ; returns degrees
   zenpos, jd, ra2, dec2 ; returns radians
   ra2 = ra2 * DRADEG
   dec2 = dec2 * DRADEG
   adiff = 90. - djs_diff_angle(ra1,dec1,ra2,dec2)

   ;----------
   ; Print a warning if indeed the Sun is above the horizon.

   if (adiff GT 0) then begin
      splog, 'WARNING: Sun above the horizon by ' $
       + string(adiff,format='(f5.1)') + ' deg for non-test exposure'
   endif

   return
end
;------------------------------------------------------------------------------
