;+
; NAME:
;   k-correct
; PURPOSE:
;   Take 5-band photometry, redshift, and desired rest-frame band and return
;     k-corrected magnitude.
; COMMENTS:
; CALLING SEQUENCE:
;   restphot= k-correct(obsphot,z)
; INPUTS:
;   obsphot   - 5-element vector containing SDSS photometry in observed frame
;   z         - redshift
; OPTIONAL INPUTS:
;   obserr    - uncertainty in obsphot
; KEYWORDS:
; OUTPUTS:
;   restphot  - 5-element vector containing SDSS photometry in rest frame
; OPTIONAL OUTPUTS:
;   resterr   - propagated uncertainty in restphot
; PROCEDURES CALLED:
; BUGS:
;   This uses an incredibly stupid interpolation scheme!
; REVISION HISTORY:
;   2000-Jun-28  Written by Hogg (IAS)
;-
