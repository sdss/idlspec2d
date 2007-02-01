; Fit a low-order spline to a set of normalized standard/model flux
; vectors. We use semi-hand picked breakpoints, chosen to:
;  - avoid any stellar&telluric features: those are fit by other templates.
;  - be dense enough under the dichroic region to fit the strongly
;    variable features (most likely caused by humidity changes).
;  - be coarse enough to not fit to any small-scale effects. The
;    notion is that there are no small-scale fluxing variations due to
;    the optics & detectors, only due to telluric effects. And those
;    are supposed to be handled by the templates.
;  
; Derived from spflux_v5's spflux_bspline.
;
function spflux_basespline, color, loglam, flux, ivar, outmask=outmask, airmass=airmass, $
                            bkpts=bkpts

    isort = sort(loglam)

                                ; Ignore all non-telescope effects
    mask1 = 1 - spflux_masklines(loglam[isort], hwidth=12.e-4, /stellar, /telluric)

    if color EQ 'b' then begin
       maxwave = 10^max(loglam)-10 ; Find better endpoint? CPL
       rawbkpt = [3810.,3855,3910,3950,$
                  4030,4150,4250,4370,4450,4550,4650,4750,4840,4950,$
                  5050,5130,5230,5330,5430,5500,5550,5600,5650,5700,5750,5800,5850,5875]
       rawbkpt = [rawbkpt, fillspan(5970.,maxwave, spacing=20)]
    end else if color EQ 'r' then begin
                                ; The points red of 8930 are within a
                                ; telluric region, and so are masked
                                ; out. But we do need a point at the
                                ; end of the spectrum, and that region
                                ; is evidently fading right at the
                                ; end, so mask it back in for our
                                ; purposes. 
       tm = where(10^loglam[isort] GT 9190, cnt)
       if cnt gt 0 then mask1[tm] = 1

       tm = where(10^loglam[isort] GE 9040 AND 10^loglam[isort] LE 9060)
       mask1[tm] = 1
       if cnt gt 0 then mask1[tm] = 1
                                ; 
       tm = where(10^loglam[isort] GE 6935 AND 10^loglam[isort] LE 6940)
       mask1[tm] = 1
       if cnt gt 0 then mask1[tm] = 1
       
       minwave = 10^min(loglam)+10 ; Find better endpoint? CPL
       rawbkpt = [6510,6600,6700,6800,6855,6937,6975, $
                  7025,7075,7140, 7380,7475,7570, 7720,7820,7920, $
                  8020,8400,8475,8600,8700,8800,8900, $
                  9050,9200]
       rawbkpt = [fillspan(minwave,5870, cnt=4), $
                  fillspan(5960.,6260, spacing=25.), $
                  fillspan(6300.,6450, spacing=30.), $
                  rawbkpt]
    end else message, 'unknown color: ', color
    fullbkpt = alog10(rawbkpt)

                                ; Fix maxrej value, which is currently
                                ; 0, as it was in spflux_v5.
    outmask1 = 0
    nord = 3
    x2 = keyword_set(airmass) ? airmass[isort] : 0
    sset = bspline_iterfit(loglam[isort], flux[isort], $
                           invvar=ivar[isort]*mask1, lower=3, upper=3, fullbkpt=fullbkpt, $
                           maxrej=ceil(0.05*n_elements(XXXXXXX)), outmask=outmask1, nord=nord, $
                           x2=x2, npoly=2*keyword_set(airmass),requiren=1)
    if (max(sset.coeff) EQ 0) then $
       message, 'B-spline fit failed!!'
    if (arg_present(outmask)) then begin
       outmask = bytarr(size(loglam,/dimens))
       outmask[isort] = outmask1
    endif
    
    bkpts = fullbkpt
    return, sset
end
