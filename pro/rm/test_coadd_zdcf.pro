; perform some tests of the coadded zdcf

pro test_coadd_zdcf

; peak significance of individual ccf
file='/data3/yshen/ftp/sdssrm/collab/prepspec/ACBFJ/ccf_output'
readcol,file,format='l,a,x,x,x,x,d', rmid, line, peak_sig

; now load the sample prop
file='/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/qso_prop-0000-56837.fits'
fits=mrdfits(file,1)
logl5100=fits[rmid].logL5100
logl5100_orig=logl5100
zall=fits[rmid].z

; correct for host L5100
file='/data3/yshen/work/agn_host/decomp_final.fits'
decomp=mrdfits(file,1)
for i=0L, n_elements(logl5100) - 1 do begin
   ind=where(decomp.fiber eq rmid[i]+1)
   if ind[0] ne -1 then logl5100[i]=decomp[ind].LOGL5100_QSO
endfor

figfile='/data3/yshen/work/composite_lag/L_z_dist_hi_ccf_peak_sig.eps'
begplot,name=figfile, /color,/landscape

rmid_uniq = rmid[uniq(rmid)]
plot,[0],[0], xrange=[0.,1.1], yrange=[43, 46], xtitle='Redshift', ytitle=textoidl('log L_{5100} [erg s^{-1}]')
oplot, fits[rmid_uniq].z, fits[rmid_uniq].logL5100, psym=1, color=cgcolor('gray')

sz=3
min=0.999 ; threshold for peak signficance in the ccf
ind=where(line eq 'ha' and peak_sig gt min and logl5100 gt 0,nn)
nobj=nn & medz=median(zall[ind]) & medl=median(logl5100[ind])
oplot, zall[ind], logL5100[ind], psym=1,color=cgcolor('green')
oplot, [medz], [medl], psym=1, symsize=sz, color=cgcolor('green')
ind=where(line eq 'hb' and peak_sig gt min and logl5100 gt 0,nn)
nobj=[nobj,nn] & medz=median(zall[ind]) & medl=median(logl5100[ind])
oplot, zall[ind], logL5100[ind], psym=symcat(6),color=cgcolor('cyan')
oplot, [medz], [medl], psym=symcat(6), symsize=sz, color=cgcolor('cyan')
ind=where(line eq 'he2' and peak_sig gt min and logl5100 gt 0,nn)
nobj=[nobj,nn] & medz=median(zall[ind]) & medl=median(logl5100[ind])
oplot, zall[ind], logL5100[ind], psym=symcat(5),color=cgcolor('blue')
oplot, [medz], [medl], psym=symcat(5), symsize=sz, color=cgcolor('blue')
ind=where(line eq 'mg2' and peak_sig gt min and logl5100 gt 0,nn)
nobj=[nobj,nn] & medz=median(zall[ind]) & medl=median(logl5100[ind])
oplot, zall[ind], logL5100[ind], psym=symcat(9),color=cgcolor('red')
oplot, [medz], [medl], psym=symcat(9), symsize=sz, color=cgcolor('red')

items=textoidl(['H\alpha','H\beta', 'HeII', 'MgII'])+' ('+string(nobj, format='(i0)')+')'
colors=cgcolor(['green','cyan','blue','red'])
sym=[1,6,5,9]
xyouts, 0.16, 0.89, 'CCF peak sig>'+string(min, format='(f0.3)'),/norm
for i=0,3 do begin
  legend, items[i], box=0, pos=[0.15, 0.87-0.05*i], /norm, color=colors[i],psym=symcat(sym[i]),textcolor=colors[i]
endfor

endplot
cgfixps, figfile

figfile='/data3/yshen/work/composite_lag/coadd_zdcf_hi_ccf_peak_sig.eps'
begplot,name=figfile, /color,/landscape
plot_layout, 4, xypos=xypos
linename=['ha','hb','he2','mg2']
lineuse=['ha','hb','he2_4686','mg2']
for i=0, 3 do begin
   ind=where(line eq linename[i] and peak_sig gt min and logl5100 gt 0,nn)
   coadd_zdcf, rmid[ind], line=lineuse[i], coadd=coadd
   ploterror, coadd.tau, coadd.wmea_dcf3,coadd.edcf_wmea3, xtitle=textoidl('\tau [Observed Days]'), ytitle='Coadded ZDCF', psym=-5, title=linename[i]

endfor
endplot
cgfixps, figfile


figfile='/data3/yshen/work/composite_lag/coadd_zdcf_15.eps'
begplot,name=figfile, /color,/landscape
linename=['ha','hb','he2','mg2']
lineuse=['ha','hb','he2_4686','mg2']
rmid_15=[101, 191, 229, 267, 272, 320, 457, 589, 645, 694, 767, 769, 775, 789, 840]
logl_15=dblarr(15) & z_15=dblarr(15)
for i=0l, 14 do begin
   ind=where(rmid eq rmid_15[i])
   if ind[0] eq -1 then print, 'RMID ', rmid_15[i], ' not found'
   logl_15[i]= logL5100[ind[0]] & z_15[i]= zall[ind[0]]
endfor
for i=0, 3 do begin
   coadd_zdcf, rmid_15, line=lineuse[i], coadd=coadd
   title=linename[i]+', med-z='+string(median(z_15),format='(f0.2)')+', med-logL='+string(median(logL_15),format='(f0.2)')
   ploterror, coadd.tau, coadd.wmea_dcf3,coadd.edcf_wmea3, xtitle=textoidl('\tau [Observed Days]'), ytitle='Coadded ZDCF', psym=-5, title=title

endfor
endplot
cgfixps, figfile



end
