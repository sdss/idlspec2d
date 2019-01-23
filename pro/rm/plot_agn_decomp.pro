; make plots for the AGN decomposition paper

pro plot_agn_decomp, f14=f14,overwrite=overwrite, old_data_only=old_data_only, $
 no_hist=no_hist

font='helvetica'

file='/data3/quasar/yshen/work/agn_host/decomp_final.fits'
;file='/data3/quasar/yshen/work/agn_host/decomp_final_vp06_logL5100_qso.fits'

result=mrdfits(file,1)

logL5100_q=result.logL5100_qso
logL5100=result.logL5100_QSO
f_H_5100=result.f_H_5100
ind=where(f_H_5100 gt 0.1 and logL5100 gt 0)
logL5100[ind]= alog10(10.D^logL5100[ind]/(1. - f_H_5100[ind]))
FWHM=result.fwhm_hb

; find those that have been successfully decomposed and measured a sigma
;ind=where(result.sigma_shen gt 70. and result.sigma_shen lt 400 and result.sigma_ok_final eq 1 and result.logL5100_qso gt 0, ngood)
ind = where(result.sigma_ok_shen eq 1,ngood)
indd=where(result.sigma_ok_shen eq 1 and result.z gt 0.6, ngood_highz)

fig1='/data3/quasar/yshen/work/agn_host/figs/L_z_dist.eps'
if ~keyword_set(no_hist) then ysize=6 else ysize=4
begplot,name=fig1,xsize=6,ysize=ysize,/encap,/color, font=font

logL5100_tot=result.logL5100_tot
if ~keyword_Set(no_hist) then begin
  pos=[0.16, 0.4, 0.96, 0.98] 
  xtickname=replicate(' ',6L)
endif else begin
  pos=[0.16, 0.16, 0.96, 0.98]
  xtitle='Redshift'
endelse
plot, result.z, logL5100_tot, xtickname=xtickname, $
  ytitle=textoidl('logL_{5100,tot} [erg s^{-1}]'),$
  psym=symcat(9), xrange=[0, 1.1], yrange=[42.7,45.8],pos=pos,/ysty,/xsty,/nodata,xtitle=xtitle
; overplot the ShenJ08 sample
file='/home/yshen/Research/Projects/EV1/data/aj268062_mrt1_ShenJ.txt'
readcol,file, format='x,d,x,x,x,x,x,x,d,d,d,d,x,d',zz,L5100,L5100_2, sigma, sigma_c, fwhm
L5100 = L5100_2 ; use the Host not subtracted AGN L
; only keep objects with measurements
ind = where(sigma gt 0 and fwhm gt 0)
L5100 = L5100[ind] & sigma = sigma[ind] & fwhm = fwhm[ind]
fwhm = fwhm*1d3
fwhm = fwhm*10.D^0.05 ; upscale by 0.05 dex to match the FWHM in Shen++(2011) for the 111 common obj
; note that L5100 in ShenJ08 needs to multiply by 5100A
logL5100 = alog10(10.D^L5100*5100.D)
oplot, zz,logL5100,psym=symcat(9,thick=3),symsize=0.3,color=cgcolor('gray')

ind = where(result.sigma_ok_shen eq 1,ngood)
indd=where(result.sigma_ok_shen eq 1 and result.z gt 0.6, ngood_highz)
if ~keyword_set(old_data_only) then begin
  oplot, result.z, logL5100_tot,  psym=symcat(9)
  oplot, result[ind].z,logL5100_tot[ind],psym=symcat(16),color=cgcolor('red')
endif
legend, 'Shen08', /norm, pos=[0.4, 0.95],box=0,psym=symcat(9,thick=3),symsize=0.3,color=cgcolor('gray'),textcolor=cgcolor('dark gray')

; overplot the Woo et al. sample
file='/data3/quasar/yshen/work/agn_host/woo_data'
readcol,file,format='d,d,d,d',zz,l5100,sigma,logmbh
oplot, zz, alog10(l5100*1d44), psym=symcat(5,thick=5),symsize=1,color=cgcolor('cyan')
legend, 'Woo et al.', /norm, pos=[0.6, 0.95],box=0,psym=symcat(5,thick=5), $
  symsize=1,color=cgcolor('cyan'),textcolor=cgcolor('cyan')
if ~keyword_set(old_data_only) then xyouts, 0.98, 43., 'A'

ntot=n_elements(result)
items=[textoidl(' N=') + string(ntot,format='(i0)'), $
   textoidl(' N_{good}=') + string(ngood,format='(i0)')]
if ~keyword_set(old_data_only) then begin
  legend, items[0],/norm,pos=[0.21,0.95],box=0,psym=symcat(9)
  legend, items[1],/norm,pos=[0.21,0.90],box=0,psym=symcat(16), color=cgcolor('red'),textcolor=cgcolor('red')
endif

if ~keyword_set(no_hist) then begin
  plothist, result.z, bin=0.1, xrange=[0,1.1],xtitle='Redshift',pos=[0.16,0.11,0.96,0.4], $
  /xsty, /noerase, ytitle=textoidl('N_{qso}'), yrange=[0,48],/ysty,xticklen=0.05
  plothist, result[ind].z,bin=0.1,/over,color=cgcolor('red')
  xyouts, 0.08, 38, 'B'
endif

endplot

if not keyword_set(F14) then $
 fig2='/data3/quasar/yshen/work/agn_host/figs/Mse_sigma_vp06.eps' else $
 fig2='/data3/quasar/yshen/work/agn_host/figs/Mse_sigma_f14.eps'

begplot,name=fig2,xsize=6,ysize=6,/encap,/color,/cmyk, font=font
sigma=result.sigma_SHEN & sigma_err=result.sigma_SHEN_err
if ~keyword_set(f14) then begin
  mass=result.logmbh_vp06 & mass_err=result.logmbh_vp06_err
  ytitle=textoidl('M_{BH,vir} (H\beta, VP06) [M')+sunsymbol()+']' 
endif else begin
  mass=result.logmbh_f14 & mass_err=result.logmbh_f14_err
  ytitle=textoidl('M_{BH,vir} (H\beta, F14) [M')+sunsymbol()+']'
endelse
thick=2
plot, [0], [0], xrange=[20, 700], yrange=[1d6,1d10], /xlog,/ylog, $
  xtitle=textoidl('\sigma_* [km s^{-1}]'),ytitle=ytitle,xsty=5, $
  ysty=5, pos=[0.16, 0.12, 0.98, 0.98],/nodata
  

; now fit to the data, taking into account errors in both measurables
; using Brandon's MCMC linear regression routine of logMbh on logsigma
outfile='/data3/quasar/yshen/work/agn_host/mcmc_reg/post_all.fits'
if file_test(outfile) eq 0 or keyword_set(overwrite) then begin
  xarr=alog10(sigma[ind]/200.D) & xerr=sigma_err[ind]/alog(10.D)/sigma[ind]
  yarr=mass[ind] & yerr=mass_err[ind]
  linmix_err, xarr,yarr, post, xsig=xerr, ysig=yerr
  mwrfits, post, outfile, /create
endif else post=mrdfits(outfile, 1)
outfile='/data3/quasar/yshen/work/agn_host/mcmc_reg/post_z_gt_0d6.fits'
if file_test(outfile) eq 0 or keyword_set(overwrite) then begin
  xarr=alog10(sigma[indd]/200.) & xerr=sigma_err[indd]/alog(10.D)/sigma[indd]
  yarr=mass[indd] & yerr=mass_err[indd]
  linmix_err, xarr,yarr, post2, xsig=xerr, ysig=yerr
  mwrfits, post2, outfile, /create
endif else post2=mrdfits(outfile,1)

alpha=median(post.alpha) & beta1=median(post.beta) & int_scat1=sqrt(median(post.sigsqr))
alpha2=median(post2.alpha) & beta2=median(post2.beta) & int_scat2=sqrt(median(post2.sigsqr))
print, alpha, beta1, int_scat1
print, alpha2, beta2, int_scat2
;oplot, [10., 1000.], 10.D^alpha*[10.,1000.]^beta1, thick=6
; now get 2sigma confidence
logsigma_arr=1. + findgen(21)*0.1
ymin=dblarr(21) & ymax=dblarr(21)
ymin2=dblarr(21) & ymax2=dblarr(21)
for i=0, 20 do begin
  logm=post.alpha + post.beta*(logsigma_arr[i] - alog10(200.D))
  ymin[i]=10.D^quantile_1d(0.025, logm)
  ymax[i]=10.D^quantile_1d(0.975, logm)
  logm=post2.alpha + post2.beta*(logsigma_arr[i] - alog10(200.D))
  ymin2[i]=10.D^quantile_1d(0.025, logm)
  ymax2[i]=10.D^quantile_1d(0.975, logm)
endfor
;polyfill, 10.D^[logsigma_arr, reverse(logsigma_arr)], [ymin2, reverse(ymax2)],noclip=0,/fill,color=cgcolor('red'), /line_fill, orientation=45
polyfill, 10.D^[logsigma_arr, reverse(logsigma_arr)], [ymin, reverse(ymax)],noclip=0,/fill,color=cgcolor('gray')
axis, xaxis=0,xrange=[20, 700],xtitle=textoidl('\sigma_* [km s^{-1}]'),/xsty,/xlog
axis, xaxis=1,xrange=[20, 700],/xsty, /xlog, xtickname=replicate(' ',6L)
axis, yaxis=0,yrange=[1d6,1d10],/ylog,/ysty, ytitle=ytitle
axis, yaxis=1,yrange=[1d6,1d10],/ylog,/ysty, ytickname=replicate(' ', 6L)

oplot, [10., 1000.], 10.D^alpha*([10.,1000.]/200.D)^beta1, thick=6
;oplot, [10., 1000.], 10.D^alpha2*[10.,1000.]^beta2, thick=6,color=cgcolor('red')
;oplot, 10.^logsigma_arr, ymin2, thick=6,color=cgcolor('red'),line=1
;oplot, 10.^logsigma_arr, ymax2, thick=6,color=cgcolor('red'),line=1

oploterror, sigma[ind], 10.D^mass[ind], sigma_err[ind], 10.D^mass[ind]*alog(10.D)*mass_err[ind], psym=symcat(9),thick=thick
; oploterror, sigma[indd], 10.D^mass[indd], sigma_err[indd], 10.D^mass[indd]*alog(10.D)*mass_err[indd], psym=symcat(9),color=cgcolor('red'),errcolor=cgcolor('red'),thick=thick
oplot, sigma[indd], 10.D^mass[indd],psym=symcat(16),color=cgcolor('red'),symsize=0.6

; Spearman's test
spear=r_correlate(sigma[ind], 10.D^mass[ind])
spear1=r_correlate(sigma[indd],10.D^mass[indd])
print,string(spear)

; Eqn 7 of KH13
oplot, [10, 1000], 0.309*([10,1000]/200.)^4.38*1d9, line=2
xyouts, 80, 2d6, 'KH13'

legend, ' all (88): r='+string(spear[0],format='(f0.2)') + $
  ', p='+string(spear[1],format='(e0.1)'),psym=symcat(9),box=0,pos=[0.18,0.95],/norm
legend, ' z>0.6 ('+string(ngood_highz,format='(i0)')+'): r='+string(spear1[0],format='(f0.2)') $
 +', p='+string(spear1[1],format='(e0.1)'),psym=symcat(16),box=0,pos=[0.18,0.9], $
  /norm,color=cgcolor('red'),textcolor=cgcolor('red')

endplot

end


pro plot_zevo_mse_sigma, overwrite=overwrite, sfont=sfont

if keyword_set(sfont) then font='helvetica'

file='/data3/quasar/yshen/work/agn_host/decomp_final.fits'
result=mrdfits(file,1)

ind=where(result.sigma_ok_shen eq 1)
result=result[ind]
ind=sort(result.z)
result=result[ind]
zz=result.z
mass=result.logmbh_vp06 & mass_err=result.logmbh_vp06_err
sigma=result.sigma_shen & sigma_err=result.sigma_shen_err
xarr=alog10(sigma/200.D) & xerr=sigma_err/alog(10.D)/sigma
yarr=mass & yerr=mass_err

file='/data3/quasar/yshen/work/agn_host/mcmc_reg/post_all.fits'
post_all=mrdfits(file,1)
alpha=median(post_all.alpha) & beta1=median(post_all.beta)

figfile='/data3/quasar/yshen/work/agn_host/figs/zevo_msigma.eps'
begplot, name=figfile, /color,/cmyk,/encap,xsize=19,ysize=5,font=font
charsize=2

; divide into 4 redshift bins
pos0=[0.05, 0.12, 0.2825, 0.98]
xtitle=textoidl('\sigma_* [km s^{-1}]')
ytitle=textoidl('M_{BH,vir} (H\beta, VP06) [M')+sunsymbol()+']'

file='/data3/quasar/yshen/work/agn_host/woo_data'
readcol,file,format='d,d,d,d',zz_woo,l5100_woo,sigma_woo,logmbh_woo

tag=['A','B','C','D']
for i=0L, 3L do begin

  dx=pos0[2]-pos0[0] & dy=0.    ; pos0[3]-pos0[1]
  ;pos=[ pos0[0] + dx*(i mod 2), pos0[1] - dy*(i/2), $
  ;      pos0[2] + dx*(i mod 2), pos0[3] - dy*(i/2) ]
  pos=[ pos0[0] + dx*(i), pos0[1] - dy*(i/2), $
        pos0[2] + dx*(i), pos0[3] - dy*(i/2) ]

  imin=i*22L & imax=(i+1)*22L - 1L 
  zmin=min(zz[imin:imax], max=zmax)
  zmed=median(zz[imin:imax])
  ztag='z='+string(zmin,format='(f0.2)')+'-'+string(zmax,format='(f0.2)')+ $
    ' (<z>=' + string(zmed,format='(f0.2)')+')'

  print, zmin, zmax, zmed

  xarr1=xarr[imin:imax] & yarr1=yarr[imin:imax]
  xerr1=xerr[imin:imax] & yerr1=yerr[imin:imax]

  print, r_correlate(xarr1, yarr1)

  outfile='/data3/quasar/yshen/work/agn_host/mcmc_reg/post_z'+string(zmin,format='(f0.2)')+'-'+string(zmax,format='(f0.2)')+'.fits'
  if file_test(outfile) eq 0 or keyword_set(overwrite) then begin
    linmix_err, xarr1,yarr1, post, xsig=xerr1, ysig=yerr1
    mwrfits, post, outfile, /create 
  endif else post=mrdfits(outfile,1)

  plot, [0], [0], xrange=[20, 700], yrange=[2d6,9d9], /xlog,/ylog, $
    xsty=5, ysty=5, pos=pos,/nodata,/noerase, charsize=charsize
  ; now get 2sigma confidence
  logsigma_arr=1. + findgen(21)*0.1
  ymin=dblarr(21) & ymax=dblarr(21)
  ymin1=dblarr(21) & ymax1=dblarr(21)
  for j=0, 20 do begin
    logm=post.alpha + post.beta*(logsigma_arr[j] - alog10(200.D))
    ; 2sigma
    ymin[j]=10.D^quantile_1d(0.025, logm)
    ymax[j]=10.D^quantile_1d(0.975, logm)
    ; 1sigma
    ymin1[j]=10.D^quantile_1d(0.16, logm)
    ymax1[j]=10.D^quantile_1d(0.84, logm)
  endfor
  polyfill, 10.D^[logsigma_arr, reverse(logsigma_arr)], [ymin, reverse(ymax)],noclip=0,/fill,color=cgcolor('sky blue')
  polyfill, 10.D^[logsigma_arr, reverse(logsigma_arr)], [ymin1, reverse(ymax1)],noclip=0,/fill,color=cgcolor('cyan')
  ;if i eq 0 or i eq 1 then xname=replicate(' ',6L) else tmp=temporary(xname)
  if i ne 0 then yname=replicate(' ',6L) else tmp=temporary(yname)
  axis, xaxis=0,xrange=[20, 700],/xsty,/xlog, xtickname=xname
  axis, xaxis=1,xrange=[20, 700],/xsty, /xlog, xtickname=replicate(' ',6L)
  axis, yaxis=0,yrange=[2d6,9d9],/ylog,/ysty, ytickname=yname
  axis, yaxis=1,yrange=[2d6,9d9],/ylog,/ysty, ytickname=replicate(' ', 6L)

  alpha2=median(post.alpha) & beta2=median(post.beta)

  oplot, [10., 1000.], 10.D^alpha*([10.,1000.]/200.D)^beta1, thick=8
  oplot, [10., 1000.], 10.D^alpha2*([10.,1000.]/200.D)^beta2, thick=8,color=cgcolor('red')
  ; Eqn 7 of KH13
  oplot, [10, 1000], 0.309*([10,1000]/200.)^4.38*1d9, line=2
  if i eq 0 then xyouts, 80, 3d6, 'KH13', charsize=charsize

  thick = 2
  oploterror, 10.D^xarr1*200.D, 10.D^yarr1, xerr1*alog(10.D)*10.D^xarr1*200.D, 10.D^yarr1*alog(10.D)*yerr1, psym=symcat(9),thick=thick, hatlength=!D.X_VSIZE / 200.
  xyouts, 25, 4d9, ztag, charsize=charsize

  if i eq 0 then begin ; overplot z=0.36 Woo data
     ind_woo=where(zz_woo le 0.4)
     oplot, sigma_woo[ind_woo], 10.D^logmbh_woo[ind_woo],psym=symcat(5,thick=5) $
       ,color=cgcolor('red')
    legend, pos=[23, 3.5d9], ' Woo et al. (z=0.36)', $
      color=cgcolor('red'),box=0,psym=symcat(5,thick=5), $
      textcolor=cgcolor('red'),charsize=charsize
  endif
  if i eq 1 then begin ; overplot z=0.57 Woo data
     ind_woo=where(zz_woo gt 0.4)
     oplot, sigma_woo[ind_woo], 10.D^logmbh_woo[ind_woo],psym=symcat(5,thick=5) $
       ,color=cgcolor('red')
     legend, pos=[23, 3.5d9], ' Woo et al. (z=0.57)', $
        color=cgcolor('red'),box=0,psym=symcat(5,thick=5), $
      textcolor=cgcolor('red'), charsize=charsize
  endif
  xyouts, 400., 3d6, tag[i], charsize=charsize

endfor

xyouts, 0.5, 0.02, align=0.5, xtitle, /norm, charsize=charsize
xyouts, 0.02, 0.5, align=0.5, ytitle, /norm, charsize=charsize, orien=90

endplot

end

; compare the results in ShenJ08
pro plot_comp_shen08

font='helvetica'

file='/data3/quasar/yshen/work/agn_host/decomp_final.fits'
result=mrdfits(file,1)


file='/home/yshen/Research/Projects/EV1/data/aj268062_mrt1_ShenJ.txt'
readcol,file, format='x,d,x,x,x,x,x,x,d,d,d,d,x,d',zz,L5100,L5100_2, sigma, sigma_c, fwhm
L5100 = L5100_2 ; use the Host not subtracted AGN L
; only keep objects with measurements
ind = where(sigma gt 0 and fwhm gt 0)
L5100 = L5100[ind] & sigma = sigma[ind] & fwhm = fwhm[ind]
fwhm = fwhm*1d3
fwhm = fwhm*10.D^0.05 ; upscale by 0.05 dex to match the FWHM in Shen++(2011) for the 111 common obj

; note that L5100 in ShenJ08 needs to multiply by 5100A
logL5100 = alog10(10.D^L5100*5100.D)

logmbh_vp06=0.91 + 0.5*alog10(10.D^logL5100/1d44) + 2.*alog10(fwhm)

figfile='/data3/quasar/yshen/work/agn_host/figs/comp_shenJ08.eps'
begplot,name=figfile,xsize=6,ysize=6.5,/encap,/color,font=font

ytitle=textoidl('M_{BH,vir} (H\beta, VP06) [M')+sunsymbol()+']'
xrange=[40,500] & yrange=[3d6, 6d9]
plot, [0], [0], xrange=xrange, yrange=yrange, /xlog,/ylog, $
  xtitle=textoidl('\sigma_* [km s^{-1}]'),ytitle=ytitle,xsty=5, $
  ysty=5, pos=[0.15, 0.12, 0.98, 0.88],/nodata

fitsfile='/data3/quasar/yshen/work/agn_host/mcmc_reg/post_all.fits'
post=mrdfits(fitsfile,1)

alpha=median(post.alpha) & beta1=median(post.beta) & int_scat1=sqrt(median(post.sigsqr))
print, alpha, beta1, int_scat1
; now get 2sigma confidence
logsigma_arr=1. + findgen(21)*0.1
ymin=dblarr(21) & ymax=dblarr(21)
ymin2=dblarr(21) & ymax2=dblarr(21)
for i=0, 20 do begin
  logm=post.alpha + post.beta*(logsigma_arr[i] - alog10(200.D))
  ymin[i]=10.D^quantile_1d(0.025, logm)
  ymax[i]=10.D^quantile_1d(0.975, logm)
endfor
polyfill, 10.D^[logsigma_arr, reverse(logsigma_arr)], [ymin, reverse(ymax)],noclip=0,/fill,color=cgcolor('gray')
axis, xaxis=0,xrange=xrange,xtitle=textoidl('\sigma_* [km s^{-1}]'),/xsty,/xlog
axis, xaxis=1,xrange=xrange,/xsty, /xlog, xtickname=replicate(' ',6L)
axis, yaxis=0,yrange=yrange,/ylog,/ysty, ytitle=ytitle
axis, yaxis=1,yrange=yrange,/ylog,/ysty, ytickname=replicate(' ', 6L)
oplot, [10., 1000.], 10.D^alpha*([10.,1000.]/200.D)^beta1, thick=6

;specify color range based on logL5100
min=42.7 & max=44.7

cgloadct, 13, /silent
thick=2
ind=where(result.sigma_ok_shen eq 1, nnn)
colors=floor( (result[ind].logL5100_tot - min)/(max-min)*255L)
colors=colors > 0 
colors=colors < 255
for i=0L, nnn-1 do begin
;  oploterror, result[ind[i]].sigma_shen, 10.D^result[ind[i]].logmbh_vp06, $
;    result[ind[i]].sigma_shen_err, 10.D^result[ind[i]].logmbh_vp06*alog(10.D)*result[ind[i]].logmbh_vp06_err, $
;    psym=symcat(16,thick=3),symsize=1, thick=thick,color=colors[i],errcolor=colors[i]
  oplot, [result[ind[i]].sigma_shen], [10.D^result[ind[i]].logmbh_vp06], $
    psym=symcat(9,thick=8),symsize=1,color=colors[i]
endfor

colors=floor( (logL5100 - min)/(max-min)*255L)
colors=colors > 0
colors=colors < 255
nnn=n_elements(sigma)
for i=0, nnn-1 do oplot,[sigma[i]], [10.D^logmbh_vp06[i]], psym=symcat(16,thick=3), $
   symsize=0.5,color=colors[i]

; now plot the Woo et al. data
file='/data3/quasar/yshen/work/agn_host/woo_data'
readcol,file,format='d,d,d,d',zz,l5100,sigma,logmbh
logL5100 = alog10(L5100*1d44)
colors=floor( (logL5100 - min)/(max - min)*255L  )
colors=colors > 0
colors=colors < 255
nnn=n_elements(sigma)
for i=0, nnn-1 do oplot,[sigma[i]], [10.D^logmbh[i]], psym=symcat(5,thick=8), $
   symsize=1., color=colors[i]

legend, box=0, ' this work', psym=symcat(9,thick=6),symsize=1,pos=[0.16, 0.86],/norm
legend, box=0, ' Shen08', psym=symcat(16,thick=3), $
   symsize=0.5,pos=[0.16,0.82],/norm
legend, box=0, ' Woo et al.', psym=symcat(5,thick=5),pos=[0.16,0.78],/norm,symsize=1.0

; Eqn 7 of KH13
oplot, [10, 1000], 0.309*([10,1000]/200.)^4.38*1d9, line=2
xyouts, 52, 4d6, 'KH13'

;pos=[0.15, 0.12, 0.98, 0.9]
pos1 = [0.15, 0.94, 0.98, 0.99]
cgcolorbar,/norm, pos = pos1, color=cgcolor('firebrick',255), range=[min, max], minor=5,ncolors=255
xyouts, 0.02, 0.91, textoidl('logL_{5100,tot}'),color=fsc_color('firebrick',255), /norm
;xyouts, 0.95, 0.91, textoidl('[erg s^{-1}]'),color=fsc_color('firebrick',255), /norm

endplot

end


pro plot_msigma_para_evo

; plot the evolution in the best-fit msigma relation

zz=[0.6,0.76, 0.26, 0.53, 0.70, 0.84]
alpha=[8.377,8.395,8.324,8.372,8.388,8.364]
alpha1=[-0.066,-0.088,-0.180,-0.138,-0.235,-0.091]
alpha2=[0.067,0.088,0.170,0.137,0.233,0.097]

beta=[1.535,1.081,1.695,1.592,1.199,0.960]
beta1=[-0.303,-0.492,-0.604,-0.660,-1.216,-0.588]
beta2=[0.304,0.511,0.614,0.677,1.231,0.598]

sca=[0.406,0.440,0.410,0.414,0.536,0.381]
sca1=[-0.031,-0.048,-0.063,-0.062,-0.083,-0.068]
sca2=[0.035,0.059,0.082,0.082,0.109,0.091]

syms=[16,9,6,6,6,6]

figfile='/data3/quasar/yshen/work/agn_host/figs/para_evo.eps'
begplot,name=figfile,xsize=6,ysize=8,/encap,/color,font=font
color=cgcolor(['black','red','cyan', 'cyan','cyan','cyan'])

pos0=[0.15, 0.68, 0.98, 0.95]
xtitle=textoidl('Redshift')
ytitle=textoidl(['\alpha', '\beta', 'Intrinsic Scatter'])
thick = 6
ticklen=0.04
xrange=[0.15,0.95]

pos=pos0
dx=pos0[2]-pos0[0] & dy=pos0[3]-pos0[1]+0.02

title=textoidl('log(M_{BH,vir}/M')+sunsymbol()+textoidl(')=\alpha+\betalog(\sigma_*/200 km s^{-1})')
plot, zz, alpha,/nodata, xrange=xrange,ytitle=ytitle[0],pos=pos,/xsty,$
  xticklen=ticklen,xtickname=replicate(' ', 6L), yrange=[8.12,8.6],/ysty,title=title
for jj=0, 5L do begin
   oploterror, zz[jj],alpha[jj],alpha1[jj],/lobar,thick=thick,psym=symcat(syms[jj]),$
    color=color[jj],errcolor=color[jj], errthick=thick
   oploterror, zz[jj],alpha[jj],alpha2[jj],/hibar,thick=thick,psym=symcat(syms[jj]),$
    color=color[jj],errcolor=color[jj], errthick=thick
endfor

pos=[ pos[0], pos[1] - dy, $
      pos[2], pos[3] - dy ]
plot, zz, beta, /nodata, xrange=xrange,ytitle=ytitle[1],pos=pos,/noerase,/xsty,$
 xticklen=ticklen,yrange=[0., 2.6],xtickname=replicate(' ', 6L),/ysty
for jj=0, 5L do begin
   oploterror, zz[jj],beta[jj],beta1[jj],/lobar,thick=thick,psym=symcat(syms[jj]),$
    color=color[jj],errcolor=color[jj],errthick=thick
   oploterror, zz[jj],beta[jj],beta2[jj],/hibar,thick=thick,psym=symcat(syms[jj]),$
    color=color[jj],errcolor=color[jj],errthick=thick
endfor

pos=[ pos[0], pos[1] - dy, $
      pos[2], pos[3] - dy ]
plot, zz, sca, /nodata, xrange=xrange,ytitle=ytitle[2],pos=pos,/noerase,/xsty,$
 xticklen=ticklen,yrange=[0.3, 0.65],xtitle='Redshift'
for jj=0, 5L do begin
   oploterror, zz[jj],sca[jj],sca1[jj],/lobar,thick=thick,psym=symcat(syms[jj]),$
    color=color[jj],errcolor=color[jj],errthick=thick
   oploterror, zz[jj],sca[jj],sca2[jj],/hibar,thick=thick,psym=symcat(syms[jj]),$
    color=color[jj],errcolor=color[jj],errthick=thick
endfor



endplot

end
