
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

