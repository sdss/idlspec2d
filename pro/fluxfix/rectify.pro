;+
; NAME:
;   rectify
;
; PURPOSE:
;   Rectify spectra using a sliding median
;
; CALLING SEQUENCE:
;   normflux = rectify(flux, [ivar, nivar=, mask=, wave=])
;
; INPUTS:
;   flux -  flux vector [npix, nspec]
;
; OPTIONAL INPUTS:
;   ivar -  inverse variance -- needed if "nivar" is returned [npix, nspec]
;   mask - if set then mask strong optical absorption feature before doing
;          sliding median.  (To do this "wave" must also be supplied.)
;          On exit this keyword will hold a copy of the flux array with 
;          masked pixels set to NaN
;   wave - wavelength vector in Angstroms -- only needed if "mask" is set
;
; OUTPUT:
;   Rectified flux vector [npix, nspec]
; 
; OPTIONAL OUTPUT:
;   nivar - inverse variance of rectified flux [npix, nspec]
;   
; COMMENTS:
;   If mask is set 16 pixels around the folling lines are masked:      
;   H-delta, Ca_K, Ca_H, G-band, H-gamma H-beta, Mgb, H-alpha
;
; BUGS:
;   Medianing is done in pixel space not wavelength space -- so for two spectra
;   to be rectified in the same manner they must have identical wavelength
;   sampling.
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_median()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   12-Aug-2003  Created by C. Tremonti, Steward Observatory
;-
;------------------------------------------------------------------------------

function rectify, flux, ivar, nivar = nivar, mask = mask, wave = wave

ny = n_elements(flux[0,*])
mflux = flux 
nflux = flux * 0
if  keyword_set(ivar) then nivar = ivar * 0

if keyword_set(mask) and keyword_set(wave) then begin
  
  linectr = [3830.0, 3889.0, $
  ;         H-delta, Ca_k,   Ca_H,   G-band,     
             4101.7, 3933.7, 3968.5, 4300, 4310, $
  ;         H-gamma H-beta, Mgb, H-alpha
             4340.5, 4861.3, 5153.0, 6562.8]
   
  for iline = 0, n_elements(linectr) - 1 do begin 
    wtindx = where(wave gt linectr[iline] - 16 and wave lt linectr[iline] + 16)
    if wtindx[0] ne -1 then mflux[wtindx,*] = 'NaN'
  endfor
  mask = mflux
endif

for y = 0, ny - 1 do begin
;  smoothcont = smooth(djs_median(mflux[*,y], width = 250, $
  smoothcont = smooth(djs_median(mflux[*,y], width = 99, $
                                 boundary = 'reflect'), 25) 
  nflux[*,y] = flux[*,y] / smoothcont
  if keyword_set(ivar) then nivar[*,y] = ivar[*,y] * smoothcont^2
endfor

return, nflux

end

