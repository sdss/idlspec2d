;+
; NAME:
;   correct_dlam
;
; PURPOSE:
;   Correct ADU/pixel to ADU/d(log-lambda)
;
; CALLING SEQUENCE:
;   correct_dlam, flux, fluxivar, wset, dlam=dlam
;
; INPUTS:
;   flux       - Flux
;   fluxivar   - Flux inverse variance
;   wset       - Wavelength coefficient trace-set
;
; OPTIONAL KEYWORDS:
;   dlam       - Log-lambda pixel size to convert to; default to 1.0d-4
;
; OUTPUTS:
;   flux       - (Modified.)
;   fluxivar   - (Modified.)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Make a map of the size of each pixel in delta-(log10-Angstroms),
;   and re-normalize the flux to ADU/d(log10-Angstroms).
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
;   04-Oct-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
pro correct_dlam, flux, fluxivar, wset, dlam=dlam, inverse=inverse

   if (NOT keyword_set(dlam)) then dlam = 1.0d-4

   traceset2xy, wset, xx, central
   traceset2xy, wset, xx-0.5, lower
   traceset2xy, wset, xx+0.5, upper

   dlogimg = abs(upper - lower)

   if keyword_set(inverse) then $
    divideflat, flux, invvar=fluxivar, (dlam/dlogimg), minval=0 $
   else $
    divideflat, flux, invvar=fluxivar, (dlogimg/dlam), minval=0

   return
end
;------------------------------------------------------------------------------
