
function atmdisp_cor, loglam, flux, plugtag, hdr, ngroff = ngroff, nrioff = nrioff, $
         title = title, mean_groff = groff_mean, mean_ngroff = ngroff_mean, $
         sig_groff = groff_sig, sig_ngroff = ngroff_sig

npix = n_elements(flux[*,0])
nfib = n_elements(flux[0,*])
std = where(strmatch(plugtag.objtype, '*R*STD*') eq 1)
wave = 10.0^loglam

;-------------------------------------------------------------------------------
; Read atmospheric dispersion vecotor matched in seeing & airmass
;-------------------------------------------------------------------------------

seeing = sxpar(hdr, 'SEEING50')
airmass = sxpar(hdr, 'AIRMASS')
airmass_fid = ['1.1', '1.2', '1.3', '1.4', '1.5']
junk = min(abs(airmass - airmass_fid), aindx)
seeing_fid = ['1.0', '1.5', '2.0', '2.5', '3.0']
junk = min(abs(seeing - seeing_fid), sindx)

atmdisp_file = filepath('atmdisp_vec_am' + airmass_fid[aindx] + 'see2.0.fit', $
               root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')

dispset = mrdfits(atmdisp_file, 1)
ndisp = n_elements(dispset.coeff[0,*])
traceset2xy, dispset, wave[*,0:ndisp-1], dispvec
mag = -2.5 * alog10(filter_thru(dispvec, wave = wave[*,0], /toair)) 
groffsynth = mag[*,1] - mag[*,2]  

;-------------------------------------------------------------------------------
; Compute (g-r) offsets of spectro and photo
;-------------------------------------------------------------------------------

flam2fnu = (wave*wave / 2.99792e18)
mag = -2.5 * alog10(filter_thru(flux * flam2fnu, wave = wave[*,0], /toair)) - 48.6 + 2.5*17.0
smag = transpose(mag[*,[1,2,3]]) 
pmag = plugtag.mag[[1,2,3]]

groff = (smag[0,*] - smag[1,*])  - (pmag[0,*] - pmag[1,*])
rioff = (smag[1,*] - smag[2,*])  - (pmag[1,*] - pmag[2,*])
xpos = plugtag.xfocal
ypos = plugtag.yfocal

;-------------------------------------------------------------------------------
; Fit 2d surface to (g-r) offsets
;-------------------------------------------------------------------------------

ok = where(abs(groff - 0.05) lt 0.5 and $
           strmatch(plugtag.objtype, '*SKY*') ne 1, nok)

; Before plot 
!P.MULTI = [0, 2, 2]
plot, xpos[ok], groff[ok], psym=3, xtitle = 'Plate X-position', title = title, $
      ytitle = '(g-r) [Spectro - Photo]', /xs, yr = [-0.2, 0.4], /nodata
djs_oplot, xpos[ok], groff[ok], psym=6, symsize=0.2, color='blue'
djs_oplot, xpos[std], groff[std], psym=6, symsize=0.2, color='red', thick=3
oplot, [-500, 500], [0, 0], thick=3

plot, ypos[ok], groff[ok], psym=3, xtitle = 'Plate Y-position', title = title, $
      ytitle = '(g-r) [Spectro - Photo]', /xs, yr = [-0.2, 0.4], /nodata
djs_oplot, ypos[ok], groff[ok], psym=6, symsize=0.2, color='blue'
djs_oplot, ypos[std], groff[std], psym=6, symsize=0.2, color='red', thick=3
oplot, [-500, 500], [0, 0], thick=3

; Do iterative fit
mask = groff * 0
mask[ok] = 1
acoeff = djs_sfit_iter(groff, xpos, ypos, 3, 3, yfit=zfit, mask=mask, maxrej = 25, $
         maxdev=0.4, lower=2.5, upper=3.0, freeiter=10, maxiter=19, outmask = outmask)

; After plot   
plot, xpos, groff, psym=3, xtitle = 'Plate X-position',  $
      ytitle = '(g-r) [Spectro - Photo]', /xs, yr = [-0.2, 0.4], /nodata
djs_oplot, xpos, groff - zfit, psym=6, symsize=0.2, color='blue'
djs_oplot, xpos[std], groff[std] - zfit[std], psym=6, symsize=0.3, color='red', thick=3
oplot, [-500, 500], [0, 0], thick=3

rej = where(outmask ne 1)
djs_oplot, xpos[rej], groff[rej] - zfit[rej], psym=2, symsize=0.5, color='green'

plot, ypos, groff, psym=3, xtitle = 'Plate Y-position', $
      ytitle = '(g-r) [Spectro - Photo]', /xs, yr = [-0.2, 0.4], /nodata
djs_oplot, ypos, groff - zfit, psym=6, symsize=0.2, color='blue'
djs_oplot, ypos[std], groff[std] - zfit[std], psym=6, symsize=0.3, color='red', thick=3
oplot, [-500, 500], [0, 0], thick=3
djs_oplot, ypos[rej], groff[rej] - zfit[rej], psym=2, symsize=0.5, color='green'

;----------------------------------------------------------------------------------
; Translate (g-r) offsets into correction vectors
;----------------------------------------------------------------------------------

zfit = zfit - median(zfit[std])
model = flux * 0.0

for ii = 0, nfib - 1 do begin
  dispseti = {func: 'legendre', xmin: 3500.0, xmax: 9400.0, coeff: fltarr(6)}
 
  for jj=0, 5 do begin
     linterp, groffsynth, reform(dispset.coeff[jj,*], ndisp), zfit[ii], coeffi
     dispseti.coeff[jj] = coeffi
  endfor
  traceset2xy, dispseti, wave[*,ii], modeli
  model[*,ii] = modeli
endfor

;----------------------------------------------------------------------------------
; Plot (g-r) histogram before/after correction
;----------------------------------------------------------------------------------

ymax = 50
!P.MULTI = [0, 2, 1]

nplothist, rioff[ok], bin=0.01, xr=[-0.4, 0.4], yr = [0, ymax], $
    xtitle = '(g-r) model - photo', ytitle = 'Number of Spectra', title = title
nplothist, rioff[ok], rixhist, riyhist, bin=0.01, color=djs_icolor('red'), /over, thick=3
nplothist, groff[ok], grxhist, gryhist, bin=0.01, color=djs_icolor('green'), /over, thick=3
oplot, [0, 0], [0, 1e4], thick=4
grfit = gaussfit(grxhist, gryhist, nterms=3, grcoef)
rifit = gaussfit(rixhist, riyhist, nterms=3, ricoef)
xyouts, -0.30, ymax*0.90, '(g-r) offset = ' + string(grcoef[1], format='(F6.3)')
xyouts, -0.30, ymax*0.80, '(g-r) sigma = ' + string(grcoef[2], format='(F6.3)')
xyouts, -0.30, ymax*0.70, '(r-i) offset = ' + string(ricoef[1], format='(F6.3)')
xyouts, -0.30, ymax*0.60, '(r-i) sigma = ' + string(ricoef[2], format='(F6.3)')

groff_mean = grcoef[1]
groff_sig = grcoef[2]

;----------------
; correct the flux

cflux = flux / model
flam2fnu = (wave*wave / 2.99792e18) 
mag = -2.5 * alog10(filter_thru(cflux * flam2fnu, wave = wave[*,0], /toair)) - 48.6 + 2.5*17.0
nsmag = transpose(mag[*,[1,2,3]]) 
      
ngroff = (nsmag[0,*] - nsmag[1,*])  - (pmag[0,*] - pmag[1,*])
nrioff = (nsmag[1,*] - nsmag[2,*])  - (pmag[1,*] - pmag[2,*])

nplothist, nrioff[ok], bin=0.01, xr=[-0.4, 0.4], yr = [0, ymax], $
    xtitle = '(g-r) model - photo', ytitle = 'Number of Spectra', title = title
nplothist, nrioff[ok], rixhist, riyhist, bin=0.01, color=djs_icolor('red'), /over, thick=3
nplothist, ngroff[ok], grxhist, gryhist, bin=0.01, color=djs_icolor('green'), /over, $
           thick=3
oplot, [0, 0], [0, 1e4], thick=4
grfit = gaussfit(grxhist, gryhist, nterms=3, grcoef)
rifit = gaussfit(rixhist, riyhist, nterms=3, ricoef)
xyouts, -0.30, ymax*0.90, '(g-r) offset = ' + string(grcoef[1], format='(F6.3)')
xyouts, -0.30, ymax*0.80, '(g-r) sigma = ' + string(grcoef[2], format='(F6.3)')
xyouts, -0.30, ymax*0.70, '(r-i) offset = ' + string(ricoef[1], format='(F6.3)')
xyouts, -0.30, ymax*0.60, '(r-i) sigma = ' + string(ricoef[2], format='(F6.3)')
xyouts, 0.0, ymax*0.95, 'Corrected', align=0.5, color=djs_icolor('red')

ngroff_mean = grcoef[1]
ngroff_sig = grcoef[2]

;----------------------------------------------------------------------------------
; Check correction as a function of plate xy

;!P.MULTI = [0, 1, 2]
;plot, xpos, zfit, psym=6, syms=0.2, /nodata, xtitle = 'Plate X-position', xr=[-300, 300], $
;      ytitle = '(g-r) Correction Residual', /xs, yr = [-0.2, 0.2]
;oplot, xpos, groff - ngroff - zfit, psym=6, symsize=0.2, color=!blue
;oplot, [-500, 500], [0, 0], thick=3

;plot, ypos, zfit, psym=6, syms=0.2, /nodata, xtitle = 'Plate Y-position', xr=[-300,300], $
;      ytitle = '(g-r) Correction Residual', /xs, yr = [-0.2, 0.2]
;oplot, ypos, groff - ngroff - zfit, psym=6, symsize=0.2, color=!blue
;oplot, [-500, 500], [0, 0], thick=3

!P.MULTI = 0

return, model

end
