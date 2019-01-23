; plot an R_L relation
; loglum and loglum_err are the log and log err of L5100 in units of erg/s
; rblr, rerr_hi and rerr_lo are the BLR in units of light days

pro plot_R_L_relation, loglum,loglum_err,rblr,rerr_hi, rerr_lo,type=type,figfile=figfile,ctype=ctype

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
xyouts, 0.18, 0.9, /norm, textoidl('Bentz et al. (2013, H\beta)'),color=cgcolor('dark gray')

; overplot the best fit in Bentz++13
k=1.527 & alpha=0.533
xx=[1d41,1d47]
yy=10.D^(k + alpha*alog10(xx/1d44))
oplot, xx/1d40, yy, color=cgcolor('dark gray')


; now overplot new data
if n_elements(loglum) gt 0 then begin
  
  lum=10.D^loglum/1d40
  lum_err_hi=lum*10.D^loglum_err - lum
  lum_err_lo=lum - lum*10.D^(-loglum_err)

  ind_uniq=uniq(type,sort(type))
  ngroup=n_elements(ind_uniq)
  color=cgcolor(['black', 'red', 'green', 'blue'])
  sym=[9,6,15,5]
  nnn=n_elements(loglum)
  flag=lonarr(nnn)
  pos=[0.17,0.89]
  for i=0L, ngroup-1 do begin
    indd=where(type eq type[ind_uniq[i]] )
    flag[indd] = i
    oploterror,[lum[indd]],[rblr[indd]],[lum_err_hi[indd]], [rerr_hi[indd]],psym=symcat(sym[i]),$
     color=color[i], errcolor=color[i],/hibar
    oploterror,[lum[indd]],[rblr[indd]],[lum_err_lo[indd]], [rerr_lo[indd]],psym=symcat(sym[i]),$
     color=color[i], errcolor=color[i],/lobar
    if ~keyword_set(ctype) then legend, type[ind_uniq[i]], /norm,pos=pos,box=0, psym=symcat(sym[i]),color=color[i],textcolor=color[i] else $
        legend, ' '+textoidl(ctype[ind_uniq[i]]), /norm,pos=pos,box=0, psym=symcat(sym[i]),color=color[i],textcolor=color[i]
    pos[1]=pos[1]-0.05
  endfor

endif


endplot

end
