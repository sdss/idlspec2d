;+
; NAME:
;     poly_fill
;
; PURPOSE:
;     Compare various interpolation methods.
;
; VERSION:
;     $Id: poly_fill.pro 50641 2014-11-17 20:39:17Z weaver $
;-
PRO poly_fill, lambda, flux, invvar, spectrum
    l = 250
    goodpts = WHERE( invvar[*,spectrum] GT 0, ngood )
    nflux = N_ELEMENTS(flux[*,spectrum])
    mingood = MIN(goodpts)
    maxgood = MAX(goodpts)
    smoothflux = SMOOTH(flux[*,spectrum],l,/EDGE_TRUNCATE)
    interpflux = djs_maskinterp(flux[*,spectrum],invvar[*,spectrum] EQ 0,/const)
    pixels = LINDGEN(nflux)
    newflux = interpflux
    IF mingood GT 0 THEN BEGIN
        damp1 = FLOAT(mingood < l)
        newflux *= 0.5*(1.0+ERF(FLOAT(pixels-mingood)/damp1))
    ENDIF
    IF maxgood LT nflux-1 THEN BEGIN
        damp2 = FLOAT(maxgood < l)
        newflux *= 0.5*(1.0+ERF(FLOAT(maxgood-pixels)/damp2))
    ENDIF
    djs_plot, lambda, flux[*,spectrum], color='white', xtitle='Wavelength [\AA]', $
        ytitle='Flux [10^{-17} erg cm^{-2} s^{-1} \AA^{-1}]'
    djs_oplot, lambda, interpflux, color='red'
    djs_oplot, lambda, newflux, color='green'
    djs_oplot, lambda, smoothflux, color='blue'
    RETURN
END
