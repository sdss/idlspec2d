;+
; NAME:
;   heliocentric
;
; PURPOSE:
;   Compute correction term to add to velocities to convert to heliocentric.
;
; CALLING SEQUENCE:
;   vhelio = heliocentric, ra, dec, [ epoch, jd=, tai=, $
;    longitude=, latitude=, altitude= ]
;
; INPUTS:
;   ra             - Right ascension [degrees]
;   dec            - Declination [degrees]
;   epoch          - Epoch of observation for RA, DEC; default to 2000.
;
; OPTIONAL KEYWORDS:
;   jd             - Decimal Julian date.
;   tai            - Number of seconds since Nov 17 1858; either JD or TAI
;                    must be specified.
;   longitude      - Longitude of observatory; default to 105.820417 deg for APO
;   latitute       - Latitude of observatory; default to 32.780361 deg for APO
;   altitude       - Altitude of observatory; default to 2788 m for APO
;
; OUTPUTS:
;   vhelio         - Velocity correction term, in km/s, to add to measured
;                    radial velocity to convert it to the heliocentric frame.
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
;   baryvel
;   zenpos
;
; REVISION HISTORY:
;   09-May-2000  Written by S. Burles & D. Schlegel
;-
;------------------------------------------------------------------------------

function heliocentric, ra, dec, epoch, jd=jd, tai=tai, $
 longitude=longitude, latitude=latitude, altitude=altitude 

   if (NOT keyword_set(epoch)) then epoch = 2000.0

   ; Default to location of Apache Point Observatory
   if (NOT keyword_set(latitude)) then latitude = 105.820417
   if (NOT keyword_set(longitude)) then longitude = 32.780361
   if (NOT keyword_set(altitude)) then altitude = 2788.

   if (NOT keyword_set(jd)) then begin
      if (keyword_set(tai)) then begin
         jd = 2400000.D + tai / (24.*60.*60.)
      endif else begin
         message, 'Must specify either JD or TAI', /cont
         return, 0
      endelse
   endif

   ; Set up common block for ZENPOS
   COMMON SITE, lat, lng, tzone
   lat = latitude
   lng = longitude
   tzone = 7 ; Not used!

   DRADEG = 180.d0 / !DPI

   ;----------
   ; Compute baryocentric velocity

   baryvel, jd, epoch, dvelh, dvelb

   ; Project velocity toward star
   vbarycen = dvelb[0]*cos(dec)*cos(ra) + $
            dvelb[1]*cos(dec)*sin(ra) + dvelb[2]*sin(dec) 

   ;----------
   ; Compute rotational velocity of observer on the Earth

   ; LAT is the latitude in radians.
   lat = latitude / DRADEG

   ; Reduction of geodetic latitude to geocentric latitude (radians).
   ; DLAT is in arcseconds.

   dlat = -(11.d0 * 60.d0 + 32.743000d0) * sin(2.d0 * lat) + $
            1.163300d0 * sin(4.d0 * lat) -0.002600d0 * sin(6.d0 * lat)
   lat  = lat + (dlat / 3600.d0) / DRADEG

   ; R is the radius vector from the Earth's center to the observer (meters).
   ; VC is the corresponding circular velocity
   ; (meters/sidereal day converted to km / sec).
   ; (sidereal day = 23.934469591229 hours (1986))

   r = 6378160.0d0 * (0.998327073d0 + 0.00167643800d0 * cos(2.d0 * lat) - $
       0.00000351d0 * cos(4.d0 * lat) + 0.000000008d0 * cos(6.d0 * lat)) + $
         altitude
   vc = 2.d0 * !DPI * (r / 1000.d0)  / (23.934469591229d0 * 3600.d0)

   ; Compute the zenith position
   zenpos, jd, zrarad, zdecrad ; returns RA,DEC of zenith in radians
   LST = zrarad * DRADEG
   HA = LST - ra

   ; Project the velocity onto the line of sight to the star.
   vrotate = vc * cos(lat) * cos(dec/DRADEG) * sin(HA/DRADEG)

   return, (-vbarycen - vrotate)
end

