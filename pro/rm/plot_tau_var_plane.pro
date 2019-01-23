; follow Keith's suggestion to plot tau against line chi^2
; highlight those with peak sig>0.99, and the 14 "detections"

pro plot_tau_var_plane

file='/data3/quasar/yshen/work/lags/prepspec/output/first_100_obj'
readcol,file,format='l,a,d,d,x,x,d',rmid,line,chi2,lag,peak_sig


figfile='/data3/quasar/yshen/work/lags/prepspec/output/tau_var.eps'
begplot,name=figfile,/encap,/color,xsize=6,ysize=4

syms=[5,9,16]
syms2=[17,16]
ind=where(line eq 'hb')
indd=where(line eq 'hb' and peak_sig gt 0.99)
pos=[0.15, 0.17,0.95,0.98]
plot,lag[ind], chi2[ind], /ylog, xtitle=textoidl('\tau_{obs}'), ytitle=textoidl('line \chi^2') $
 , psym=symcat(syms[0]),pos=pos,ytickname=textoidl(['10','100','10^3','10^4'])
oplot, lag[indd],chi2[indd],psym=symcat(syms[0]),color=cgcolor('red')

ind=where(line eq 'mg2')
indd=where(line eq 'mg2' and peak_sig gt 0.99)
oplot,lag[ind], chi2[ind],psym=symcat(syms[1])
oplot, lag[indd], chi2[indd],psym=symcat(syms[1]),color=cgcolor('red')

items=['Hbeta','MgII']
legend, pos=[0.2, 0.94], items[0],/norm, box=0,psym=symcat(syms[0])
legend, pos=[0.2, 0.89], items[1],/norm, box=0,psym=symcat(syms[1])

; now overplot the 14 "detections"
id=[775,797,272,320,252,377,160,101,229,767,694,519,457,140]
line_d=['hb','hb','hb','hb','hb','hb','hb','mg2','mg2','mg2','hb','hb','mg2','mg2']

nnn=n_elements(id)
for i=0,nnn-1 do begin
  ind=where(rmid eq id[i] and line eq line_d[i])
  if line_d[i] eq 'hb' then symind=0 
  if line_d[i] eq 'mg2' then symind=1
  oplot,lag[ind], chi2[ind],psym=symcat(syms2[symind])

endfor


endplot

end
