;+
; NAME:
;   radec_to_xyfocal
; PURPOSE:
;   Take RA and DEC values and return XFOCAL and YFOCAL for a plate
; CALLING SEQUENCE:
;   radec_to_xyfocal, ra, dec, xfocal, yfocal, racen=racen, $
;     deccen=deccen, airtemp=airtemp, mjd=mjd
; INPUTS:
;   ra,dec     - [N] arrays of locations (J2000 deg)
;   racen      - RA center for tile (J2000 deg)
;   deccen     - DEC center for tile (J2000 deg)
; OPTIONAL INPUTS:
;   lst        - LST of observation (defaults to racen)
;   airtemp    - Design temperature (in C, default to 5)
; OUTPUTS:
;   xfocal, yfocal - [N] arrays of focal plane positions
; COMMENTS:
;   Designed for the SDSS 2.5m at APO
;   PA=0 ALWAYS
; BUGS:
;   ALIGNMENT FIBERS HAVE OFFSET?????
; REVISION HISTORY:
;   26-Oct-2006  Written by MRB, NYU
;-
;------------------------------------------------------------------------------
pro radec_to_xyfocal, ra, dec, xfocal, yfocal, racen=racen, deccen=deccen, $
                      airtemp=airtemp, lst=lst, norefrac=norefrac, $
                      nodistort=nodistort

platescale = 217.7358D           ; mm/degree

;; from $PLATE_DIR/test/plParam.par
rcoeffs=[-0.000137627D, -0.00125238D, 1.5447D-09, 8.23673D-08, $
         -2.74584D-13, -1.53239D-12, 6.04194D-18, 1.38033D-17, $
         -2.97064D-23, -3.58767D-23] 

;; deal with atmospheric refraction
if(NOT keyword_set(norefrac)) then begin
    apo_refrac, ra, dec, racen, deccen, ra_refrac, dec_refrac, lst=lst, $
      airtemp=airtemp
endif else begin
    ra_refrac=ra
    dec_refrac=dec
endelse

;; convert to focal coordinates
;;radec_to_munu, ra_refrac, dec_refrac, mu, nu, node=racen-90.D,
;;incl=deccen
dradeg = 180.0d0/4.0/atan(1.0d0)
dtheta = djs_diff_angle(ra_refrac, dec_refrac, racen, deccen)
cosdtheta = cos(dtheta/dradeg)
fac = dtheta/tan(dtheta/dradeg)
r1 = ra_refrac/dradeg
d1 = dec_refrac/dradeg
d0 = deccen/dradeg
r0 = racen/dradeg
xsky = cos(d1) * sin(r1-r0)/cosdtheta * fac
ysky = (cos(d0)*sin(d1) - sin(d0)*cos(d1)*cos(r1-r0))/cosdtheta * fac
;recenter, ra_refrac, dec_refrac, racen, deccen, dlon, dlat
;w = where(dlon GE 180.0, nw)
;if (nw GT 0) then dlon[w] = dlon[w] - 360.0
xfocal = xsky*platescale
yfocal = ysky*platescale

;; apply radial distortions
if(NOT keyword_set(nodistort)) then begin
    rfocal= double(sqrt(xfocal^2+yfocal^2))
    correction=replicate(rcoeffs[0], n_elements(rfocal))
    for i=1L, n_elements(rcoeffs)-1L do begin
        correction=correction+rcoeffs[i]*((double(rfocal))^(double(i)))
    endfor
    rnew= rfocal+correction
    rscale=rnew/rfocal
    
    xfocal=xfocal*(correction/rfocal+1.D)
    yfocal=yfocal*(correction/rfocal+1.D)
endif

return
end
