;+
; NAME:
;   correct_dlam
;
; PURPOSE:
;   Correct ADU/pixel to ADU/dlam
;
; CALLING SEQUENCE:
;   correct_dlam, flux, fluxivar, wset, dlam=dlam
;
; INPUTS:
;   flux       - object counts
;   fluxivar   - inverse variance
;   wset       - wavelength coefficient traceset
;
; OPTIONAL KEYWORDS:
;   dlam       - log lambda pixel size to convert to (1.0e-4 default)
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   flux and fluxivar are overwritten
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   divideflat
;   traceset2xy
;
; REVISION HISTORY:
;    4-Oct-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
pro correct_dlam, flux, fluxivar, wset, dlam=dlam

     if NOT keyword_set(dlam) then dlam=1.0e-4

     traceset2xy,wset, xx, central
     traceset2xy,wset, xx-0.5, lower
     traceset2xy,wset, xx+0.5, upper

     dlogimg = abs(upper - lower)
     divideflat, flux, fluxivar, (dlogimg/dlam), minval=0

return
end
   
