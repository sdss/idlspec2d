;+
; NAME:
;   jeg_sphoto_coef
; 
; PURPOSE:
;   Compute coefficients which can be used to correct the spectrophotometry
;
; CALLING SEQUENCE:
;   coef = jeg_sphoto_coef(wave, flux, ivar, phot_mag, /plot, $
;                          x = x, xavg=xavg, x2avg=x2avg)
;
; INPUTS:
;   wave       - wavelength array [npix]
;   flux       - 2d flux array in 10-17 erg/s/cm^2/A [npix, nfiber]
;   ivar       - inverse variance [npix, nfiber]
;   phot_mag   - photo fiber mag [5, nfiber]
;
; OUTPUTS:
;   an array of coefficients for reconstructing the loss function [3, nfiber]
;     loss = 10.0^(-0.4 *(coef[0] + coef[1]*x + coef[2]x^2))
;     spec_obs = spec_true * loss 
;
; OPTIONAL OUTPUTS:
;   x      - the array of x values used to compute the loss [npix] (see paper
;            by Jim Gunn for full explanation)
;   xavg   - array of average x values for the g,r,i filters [3, npix] 
;   x2avg  - array of average x^2 values for the g,r,i filters [3, npix] 
;   doplot - if set a plot is generated showing the original spectrum in white,
;            the corrected spectrum in blue, the loss function in red, and
;            the 3 points used to compute the loss function in green
;
; COMMENTS:
;   See paper by Jim Gunn for a full explanation
;
; BUGS: 
;   Not sure what happens when the mags are bad
;
; PROCEDURES CALLED:
;   filter_thru
;   linterp
; 
; REVISION HISTORY:
;   2-Feb-2004  Written by C. Tremonti, U of Az
;-
;------------------------------------------------------------------------------

function jeg_sphoto_coef, wave, flux, ivar, phot_mag, doplot = doplot, $
         xavg = xavg, x2avg=x2avg

lam0 = 5018.0
d3900 = (lam0^2/3900.0^2 - 1) 
d9000 = (lam0^2/9000.0^2 - 1) 
x = (d3900 + d9000 + 2*(lam0^2/wave^2 - 1)) / (d3900 - d9000)

nfiber = n_elements(flux[0,*])
flambda2fnu = (wave*wave / 2.99792e18) # replicate(1,nfiber)
x2d = x # replicate(1,nfiber)

norm =  transpose(filter_thru(flux*flambda2fnu, waveimg=wave, $
                  mask=(ivar LE 0))) 

spec_mag = -2.5 * alog10(norm) + 2.5*17 -48.6

xavg =  transpose(filter_thru(flux*flambda2fnu*x2d, waveimg=wave, $
                    mask=(ivar LE 0))) / norm

x2avg =  transpose(filter_thru(flux*flambda2fnu*x2d^2, waveimg=wave, $
                     mask=(ivar LE 0))) / norm

coef = fltarr(3, nfiber)
 
for ii = 0, nfiber - 1 do begin

  a = [[1, xavg[1,ii], x2avg[1,ii]], $
       [1, xavg[2,ii], x2avg[2,ii]], $
       [1, xavg[3,ii], x2avg[3,ii]]]

  b = [spec_mag[1,ii] - phot_mag[1,ii], $
       spec_mag[2,ii] - phot_mag[2,ii], $
       spec_mag[3,ii] - phot_mag[3,ii]]

  svdc, a, w, u, v
  coef[*,ii] = svsol(u, w, v, b)

  if keyword_set(doplot) then begin
    loss = 10.0^(-0.4 *(coef[0,ii] + coef[1,ii]*x + coef[2,ii]*x^2))
    norm = median(flux[*,ii])
    plot, wave, smooth(flux[*,ii]/norm, 9), yr=[0, 3], /xs, nsum=9,$
          xtitle = 'Wavelength (A)', ytitle = 'Normalized Flux'
    oplot, wave, smooth(flux[*,ii]/loss/norm, 9), color=!blue, nsum=9
    oplot, wave, loss, color=!red, nsum=9

    linterp, x, wave, xavg[[1,2,3],ii], wlavg
    cfac = 10.0^((spec_mag[[1,2,3],ii] - phot_mag[[1,2,3],ii])/(-2.5))
    oplot, wlavg, cfac, psym=2, symsize = 3, thick=3, color=!green
  endif
endfor

return, coef

end
