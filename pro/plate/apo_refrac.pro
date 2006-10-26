;+
; NAME:
;   apo_refrac
; PURPOSE:
;   Take RA and DEC and apply refraction at APO for a specific LST
; CALLING SEQUENCE:
;   apo_refrac, ra, dec, racen, deccen, ra_refrac, dec_refrac [, $
;      lst=, airtemp= ]
; INPUTS:
;   ra,dec     - [N] arrays of true locations (J2000 deg)
;   racen      - RA center for tile (J2000 deg)
;   deccen     - DEC center for tile (J2000 deg)
; OPTIONAL INPUTS:
;   lst        - anticipated LST of observation (in deg; default to racen)
;   airtemp    - air temperature (in C; default to 5)
; OUTPUTS:
;   ra_refrac,dec_refrac - [N] arrays of observed locations (J2000 deg)
; COMMENTS:
;   Without setting LST explicitly, it assumes that you will observe
;   at an hour angle of 0. 
;   Uses Lat=32.7803 deg, and Height=2788 m (relevant for APO)
; REVISION HISTORY:
;   26-Oct-2006  Written by MRB, NYU (cribbed from PRIMUS code by Burles)
;-
;------------------------------------------------------------------------------
pro apo_refrac, ra, dec, racen, deccen, ra_refrac, dec_refrac, $
                lat=lat, airtemp=airtemp, lst=lst

if n_elements(lat) EQ 0 then lat = 32.7803D
if n_elements(lst) EQ 0 then lst = double(racen)
if n_elements(height) EQ 0 then height = 2788.D
if n_elements(airtemp) EQ 0 then airtemp = 5.D
airtemp_k=airtemp+273.155  ; C to Kelvin

ha = lst - ra
hacen = lst - racen

hadec2altaz, ha, dec, lat, alt, az
new_alt = co_refract_list(alt, altitude=height, $
                          epsilon=0.00001, temperature=airtemp_k, $
                          /to_observed, pressure=1013.25)

hadec2altaz, hacen, deccen, lat, altcen, azcen
new_altcen = co_refract_list(altcen, altitude=height, $
                             epsilon=0.00001, temperature=airtemp_k, $
                             /to_observed, pressure=1013.25)

alt_refrac = new_alt - (new_altcen - altcen)
altaz2hadec, alt_refrac, az, lat, ha_refrac, dec_refrac

ha_refrac = 180.D/!DPi * asin(sin(ha_refrac*!DPi/180.D))

ra_refrac = lst - ha_refrac

end
