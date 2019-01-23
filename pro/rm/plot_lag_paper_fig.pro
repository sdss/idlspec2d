; plot an R_L relation
; loglum and loglum_err are the log and log err of L5100 in units of erg/s
; rblr, rerr_hi and rerr_lo are the BLR in units of light days

pro plot_R_L_relation, loglum,loglum_err,rblr,rerr_hi, rerr_lo,type=type,figfile=figfile

; first read the R-L data from Bentz et al.
file='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/tables/Bentz13_data'
readcol,file, format='a,d,d,d,x,x,d,d',obj,tau,err_plus,err_minus,L5100,L5100_err

L5100=10.D^L5100/1d40  ;& L5100_err=L5100_err*alog(10.D)*L5100
L5100_err_hi=L5100*10.D^L5100_err - L5100
L5100_err_lo=L5100 - L5100*10.D^(-L5100_err)

if ~keyword_set(figfile) then figfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/R_L_relation.eps'
begplot, name=figfile,/color,/encap, /cmyk,ysize=5
ang = string(197B)
thick=2. & symsize=0.5
pos=[0.12,0.13, 0.95,0.98]
plot, L5100, tau, /xlog,/ylog, xrange=[2d41,1d46]/1d40,yrange=[1,400],/xsty,/ysty, $
 /nodata, xtitle=textoidl('\lambdaL_{\lambda}(5100 ')+ang+')'+textoidl(' [erg s^{-1}]'), $
 ytitle=textoidl('BLR Size [light days]'),pos=pos,xtickname=textoidl(['10^{42}','10^{43}','10^{44}','10^{45}', '10^{46}']),xticklen=0.032

;print, L5100_err
;oploterror, L5100, tau, L5100_err, replicate(0., n_elements(tau)),psym=symcat(9,thick=thick),symsize=symsize,thick=thick,/nohat
oploterror, L5100, tau, L5100_err_hi,err_plus, psym=symcat(9,thick=thick),/nohat, /hibar,thick=thick,symsize=symsize,color=cgcolor('gray'),errcolor=cgcolor('dark gray')
oploterror, L5100, tau, L5100_err_lo,err_minus, psym=symcat(9,thick=thick),/nohat, /lobar,thick=thick,symsize=symsize,color=cgcolor('gray'),errcolor=cgcolor('dark gray')
xyouts, 0.18, 0.9, /norm, 'Bentz et al. (2013)',color=cgcolor('dark gray')

; overplot the best fit in Bentz++13
k=1.527 & alpha=0.533
xx=[1d41,1d47]
yy=10.D^(k + alpha*alog10(xx/1d44))
oplot, xx/1d40, yy, color=cgcolor('dark gray')


; now overplot new data
symsize=1.2
if n_elements(loglum) gt 0 then begin
  
  lum=10.D^loglum/1d40
  lum_err_hi=lum*10.D^loglum_err - lum
  lum_err_lo=lum - lum*10.D^(-loglum_err)

  ind_uniq=uniq(type,sort(type))
  ngroup=n_elements(ind_uniq)
  color=cgcolor(['black', 'red', 'green', 'blue'])
  sym=[16,6,15,5]
  nnn=n_elements(loglum)
  flag=lonarr(nnn)
  pos=[0.17,0.89]
  for i=0L, ngroup-1 do begin
    indd=where(type eq type[ind_uniq[i]] )
    flag[indd] = i
    oploterror,[lum[indd]],[rblr[indd]],[lum_err_hi[indd]], [rerr_hi[indd]],psym=symcat(sym[i]),$
     color=color[i], errcolor=color[i],/hibar,symsize=symsize
    oploterror,[lum[indd]],[rblr[indd]],[lum_err_lo[indd]], [rerr_lo[indd]],psym=symcat(sym[i]),$
     color=color[i], errcolor=color[i],/lobar,symsize=symsize
    legend, type[ind_uniq[i]], /norm,pos=pos,box=0, psym=symcat(sym[i]),color=color[i],textcolor=color[i]
    pos[1]=pos[1]-0.05
  endfor

endif


endplot

end

pro plot_M_RM_SE, logvp, logvp_err, logMse, logMse_err, sigma_rms, sig_err,fwhm_mean,fwhm_err, type=type,fvir=fvir

if ~keyword_set(fvir) then fvir=5.5

if ~keyword_set(figfile) then figfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/Mrm_Mse.eps'
begplot, name=figfile,/color,/encap, xsize=4, ysize=4

pos=[0.15,0.17, 0.98,0.98]
range=[6,9]
plot, [0],[0], /nodata,xrange=range, yrange=range,xtitle=textoidl('log M_{RM} (f=5.5) [M')+sunsymbol()+']', $
 ytitle=textoidl('log M_{SE,H\beta} [M')+sunsymbol()+']',pos=pos

oplot,range,range,line=2
ind_uniq=uniq(type,sort(type))
ngroup=n_elements(ind_uniq)
color=cgcolor(['black', 'red', 'green', 'blue'])
sym=[16,6,15,5]

pos=[0.17,0.92]
for i=0L, ngroup-1 do begin
  indd=where(type eq type[ind_uniq[i]] )
  oploterror, [logvp[indd]+alog10(fvir)], logMse[indd], logvp_err[indd], logMse_err[indd],psym=symcat(sym[i]),$
  color=color[i], errcolor=color[i]
  legend, type[ind_uniq[i]], /norm,pos=pos,box=0, psym=symcat(sym[i]),color=color[i],textcolor=color[i]
    pos[1]=pos[1]-0.05
endfor
endplot

figfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/sig_rms_fwhm_mean.eps'
begplot, name=figfile,/color,/encap, xsize=4, ysize=4
pos=[0.2,0.17, 0.98,0.98]
xrange=[2.8,3.35] & yrange=[3,4]
plot, [0],[0], /nodata,xrange=xrange, yrange=yrange,/xsty,/ysty, $
  xtitle=textoidl('log\sigma_{rms} [km s^{-1}]'), $
 ytitle=textoidl('log FWHM_{mean,H\beta} [km s^{-1}]'), pos=pos ;, ytickname=textoidl(['1000','10^4'])
pos=[0.22,0.92]
color=cgcolor(['black', 'red', 'green', 'blue'])
for i=0L, ngroup-1 do begin
  indd=where(type eq type[ind_uniq[i]] )
  oploterror, [alog10(sigma_rms[indd])], alog10(fwhm_mean[indd]),(sig_err[indd]/sigma_rms[indd])/alog(10.D), fwhm_err[indd]/fwhm_mean[indd]/alog(10.D),psym=symcat(sym[i]),$
  color=color[i], errcolor=color[i]
  ;oploterror, [(sigma_rms[indd])],(fwhm_mean[indd]), (sig_err[indd]), fwhm_err[indd],psym=symcat(sym[i]),$
  ;color=color[i], errcolor=color[i]
 
  legend, type[ind_uniq[i]], /norm,pos=pos,box=0, psym=symcat(sym[i]),color=color[i],textcolor=color[i]
    pos[1]=pos[1]-0.05
endfor
endplot


end
