; plot the coadded CCF

pro plot_coadd_ccf


coadd_ccf, lindgen(849), line='hb', result=result0
coadd_ccf, lindgen(849), line='he2', result=result1
coadd_ccf, lindgen(849), line='ha', result=result2


file = getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
target=mrdfits(file,1)
target=target[0:848]
ind=where(target.zfinal lt 0.96)
rmid=ind
coadd_ccf, rmid, line='mg2', result=result3

file='/data3/yshen/ftp/sdssrm/collab/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/qso_prop-0000-56837_lineshift.fits'
info=mrdfits(file,1)
logl5100=info.logl5100


fig='/data3/yshen/work/composite_lag/coadded_ccf.ps'
begplot, name=fig, /color, /landscape

plot_layout, 4, xypos=xypos, pmargin=[0.2,0.2], omargin=[0.1, 0.0, 1.05, 0.98]

ploterror, result0.tau, result0.ccf_coadd, result0.ccf_coadd_err, psym=5, $
  xtitle=textoidl('\tau'), ytitle='Coadded correlation', pos=xypos[*,0], yrange=[-0.1, 0.3]
oploterror, result1.tau, result1.ccf_coadd, result1.ccf_coadd_err, psym=5,color=cgcolor('green')
xyouts, -90, 0.25, textoidl('H\beta')
xyouts, -90, 0.20, textoidl('HeII'),color=cgcolor('green')

ploterror, result0.tau, result0.ccf_coadd, result0.ccf_coadd_err, psym=5, $
  xtitle=textoidl('\tau'), ytitle='Coadded correlation', pos=xypos[*,1], yrange=[-0.1, 0.3],/noerase
oploterror, result2.tau, result2.ccf_coadd, result2.ccf_coadd_err, psym=5,color=cgcolor('red')
;xyouts, [-90, 0.25], textoidl('H\beta')
xyouts, -90, 0.20, textoidl('H\alpha'),color=cgcolor('red')

ploterror, result0.tau, result0.ccf_coadd, result0.ccf_coadd_err, psym=5, $
  xtitle=textoidl('\tau'), ytitle='Coadded correlation', pos=xypos[*,2], yrange=[-0.1, 0.3],/noerase
oploterror, result3.tau, result3.ccf_coadd, result3.ccf_coadd_err, psym=5,color=cgcolor('cyan')
;xyouts, [-90, 0.25], textoidl('H\beta')
xyouts, -90, 0.20, textoidl('MgII'),color=cgcolor('cyan')


; now plot lum-divided sample
rmid=result0.rmid
lum_use=logL5100[rmid]
ind1=where(lum_use lt median(lum_use), complement=ind2)
print, 'Delta_logL=', median(lum_use[ind2]) - median(lum_use[ind1])
coadd_ccf, rmid[ind1], line='hb', result=result_lol
coadd_ccf, rmid[ind2], line='hb', result=result_hil
ploterror, result_lol.tau, result_lol.ccf_coadd, result_lol.ccf_coadd_err, psym=5, $
  xtitle=textoidl('\tau'), ytitle='Coadded correlation', pos=xypos[*,0], yrange=[-0.1, 0.3]
oploterror, result_hil.tau, result_hil.ccf_coadd, result_hil.ccf_coadd_err, psym=5,color=cgcolor('green')
xyouts, -90, 0.25, textoidl('H\beta, lo-L')
xyouts, -90, 0.20, textoidl('H\beta, hi-l'),color=cgcolor('green')

rmid=result1.rmid
lum_use=logL5100[rmid]
ind1=where(lum_use lt median(lum_use), complement=ind2)
print, 'Delta_logL=', median(lum_use[ind2]) - median(lum_use[ind1])
coadd_ccf, rmid[ind1], line='he2', result=result_lol
coadd_ccf, rmid[ind2], line='he2', result=result_hil
ploterror, result_lol.tau, result_lol.ccf_coadd, result_lol.ccf_coadd_err, psym=5, $
  xtitle=textoidl('\tau'), ytitle='Coadded correlation', pos=xypos[*,1], yrange=[-0.1, 0.3],/noerase
oploterror, result_hil.tau, result_hil.ccf_coadd, result_hil.ccf_coadd_err, psym=5,color=cgcolor('green')
xyouts, -90, 0.25, textoidl('HeII, lo-L')
xyouts, -90, 0.20, textoidl('HeII, hi-l'),color=cgcolor('green')

rmid=result2.rmid
lum_use=logL5100[rmid]
ind1=where(lum_use lt median(lum_use), complement=ind2)
print, 'Delta_logL=', median(lum_use[ind2]) - median(lum_use[ind1])
coadd_ccf, rmid[ind1], line='ha', result=result_lol
coadd_ccf, rmid[ind2], line='ha', result=result_hil
ploterror, result_lol.tau, result_lol.ccf_coadd, result_lol.ccf_coadd_err, psym=5, $
  xtitle=textoidl('\tau'), ytitle='Coadded correlation', pos=xypos[*,2], yrange=[-0.1, 0.3],/noerase
oploterror, result_hil.tau, result_hil.ccf_coadd, result_hil.ccf_coadd_err, psym=5,color=cgcolor('green')
xyouts, -90, 0.25, textoidl('H\alpha, lo-L')
xyouts, -90, 0.20, textoidl('H\alpha, hi-l'),color=cgcolor('green')

rmid=result3.rmid
lum_use=logL5100[rmid]
ind1=where(lum_use lt median(lum_use), complement=ind2)
print, 'Delta_logL=', median(lum_use[ind2]) - median(lum_use[ind1])
coadd_ccf, rmid[ind1], line='mg2', result=result_lol
coadd_ccf, rmid[ind2], line='mg2', result=result_hil
ploterror, result_lol.tau, result_lol.ccf_coadd, result_lol.ccf_coadd_err, psym=5, $
  xtitle=textoidl('\tau'), ytitle='Coadded correlation', pos=xypos[*,3], yrange=[-0.1, 0.3],/noerase
oploterror, result_hil.tau, result_hil.ccf_coadd, result_hil.ccf_coadd_err, psym=5,color=cgcolor('green')
xyouts, -90, 0.25, textoidl('MgII, lo-L')
xyouts, -90, 0.20, textoidl('MgII, hi-l'),color=cgcolor('green')

endplot

cgfixps, fig


end
