function heliocentric, ra, dec, mjd, $
      longitude=longitude, latitude=latitude, altitude=altitude 

   ; ra, dec, longitude and latitude are in degrees
   ; MJD is probably not the best thing to use

   ; Longitude for APO
  if (NOT keyword_set(latitude)) then latitude=105.82042

   ; Latitude for APO
  if (NOT keyword_set(latitude)) then latitude= 32.780361

   ; Altitude for APO
  if (NOT keyword_set(altitude)) then altitude=4200.0
;   first get baryocentric velocity:

   jd = mjd + 0.5d
   baryvel, jd, 2000.0, vh, vb

;   Now estimate rotation of observer

; LAT is the latitude in radians.
   lat = latitude * !Pi/180.d0

; Reduction of geodetic latitude to geocentric latitude (radians).
; Dlat is in arcseconds.

   dlat = -(11.d0 * 60.d0 + 32.743000d0) * sin(2.d0 * lat) + $
            1.163300d0 * sin(4.d0 * lat) -0.002600d0 * sin(6.d0 * lat)
   lat  = lat + (dlat / 3600.d0) * !Pi / 180.d0

; R is the radius vector from the Earth's center to the observer (meters).
; Vc is the corresponding circular velocity
; (meters/sidereal day converted to km / sec).
; (sidereal day = 23.934469591229 hours (1986))

   r = 6378160.0d0 * (0.998327073d0 + 0.00167643800d0 * cos(2.d0 * lat) - $
       0.00000351d0 * cos(4.d0 * lat) + 0.000000008d0 * cos(6.d0 * lat)) + $
         altitude
   vc = 2.d0 * !Pi * (r / 1000.d0)  / (23.934469591229d0 * 3600.d0)

; Project the velocity onto the line of sight to the star.
      lmst = ast_mst(epoch,longitude)
      v = vc * cos(lat) * cos(dec*!Pi/180.d0) * sin((ra-lmst)*!Pi/180.d0)

      return
      end



