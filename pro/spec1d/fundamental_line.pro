;+
; NAME:
;   fundamental_line
; PURPOSE:
;   Crush the fundamental plane to a measly line.
; CALLING SEQUENCE:
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS
; OUTPUTS:
; COMMENTS:
; BUGS:
;   Using RA and Dec for matches rather than OBJID because tsObj files
;     sometimes busted; search code for "HACK".
; EXAMPLES:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   25-Jun-2000  Written by Hogg (IAS)
;-
pro fundamental_line, spname=spname,tsname=tsname,omegam=omegam,omegal=omegal
;; set cosmographic parameters
  if NOT keyword_set(omegam) then omegam= 0.2
  if NOT keyword_set(omegal) then omegal= 0.0
;; set aperture radii for photo radial profile outputs, in arcsec
  photorad= [0.23,0.68,1.03,1.76,3.00,4.63,7.43,11.42,18.20,28.20,44.21, $
             69.00,107.81,168.20,263.00]
  maxnrad= 15
;; set aperture radii for measurements, in arcsec at one Hubble length!
  physrad= [1.0,2.0,3.0]
;; if no filenames set, use debugging defaults
  if NOT keyword_set(spname) then $
    spname= '/data/spectro/2d_3c/0306/2dnew/spPlate-0306-51690.fits'
  if NOT keyword_set(tsname) then $
    tsname= '/data/spectro/plates/tsObj-51637-306.fit'
;; read object IDs and the photo outputs (tsObj)
; HACK: USING RA AND DEC BECAUSE TSOBJ ARE SOMETIMES BUSTED
  plugmap= mrdfits(spname,4)
  ra= plugmap.ra
  dec= plugmap.dec
  tsobj= mrdfits(tsname,1)
;; load redshifts
;; begin loop over spectra
  for spectrum=0,n_elements(plugmap)-1 do begin
    z= 0.20
    help, z
;; compute aperture radii in arcsec, from redshifts and cosmography
    rad= physrad/angdidis(z,omegam,omegal)
;; identify relevant entry in the tsObj file
    obj= where((abs(tsobj.dec-dec[spectrum]) LT (1.0/3600.0)) AND $
               (abs(tsobj.ra-ra[spectrum]) LT $
                          (1.0/(3600.0*cos(dec[spectrum]*!DPI/180.0)))),nobj)
    obj= obj[0]
    help, nobj
;; check that the object exists, the object is unique, and that
;; there are enough entries in the i and r-band radial profiles
    if (nobj EQ 1) AND $
       (photorad[(tsobj.nprof)[3,obj]] GE max(rad)) AND $
       (photorad[(tsobj.nprof)[4,obj]] GE max(rad)) then begin
;; compute aperture luptitudes with interpolation
      lupr= interpol((tsobj.profmean)[*,3,obj],photorad,rad)
      lupi= interpol((tsobj.profmean)[*,4,obj],photorad,rad)
;; compute 2 (fairly) independent surface brightnesses
      sb= 10.0^(0.4*(24.0-lupr))
      print, sb
      area= !DPI*rad^2
      pseudosky= (sb[2]*area[2]-sb[1]*area[1])/(area[2]-area[1])
      help, pseudosky
      sbinner= sb[0]-pseudosky
      sbouter= (sb[1]*area[1]-sb[0]*area[0])/(area[1]-area[0])-pseudosky
help, sbinner,sbouter
    endif
stop
  endfor
;; done
  return
end
