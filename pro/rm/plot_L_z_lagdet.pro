; plot the distribution in the L5100-z plane for 
; different samples: local RM AGN, SDSS-RM lag detections, etc


pro plot_L_z_lagdet, bennert=bennert,hst=hst,plot_mass=plot_mass

figfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/L_z_dist.eps'
begplot,name=figfile,xsize=7,ysize=5,/encap,/color

red, omegalambda=0.7,omega0=0.3,h100=0.7
zlist = indgen(101)*0.1
t_lb_list = getage(0,/gyr) - getage(zlist, /gyr)
t_arr = [0,1,2,3,4,5,6,getage(0,/gyr) - getage(0.8, /gyr)]
xtickname = [string(t_arr[0:6], format='(i0)'),string(t_arr[7],format='(f0.1)') ]
xtickv = interpol(zlist,t_lb_list,t_arr)


plot,[0],[0],/nodata,xrange=[0,0.8],yrange=[42, 46],xtitle='Redshift', $
  ytitle=textoidl('logL_{5100,AGN} [erg s^{-1}]'),pos=[0.12,0.12,0.95,0.9],xsty=5

axis, xaxis=0, xrange=[0,0.8],/xsty, xtitle = 'Redshift'
axis, xaxis=1, xrange=[0,0.8], xtitle = textoidl('Lookback Time [Gyr]'), xticks=7, $
   xtickv=xtickv, xtickname = xtickname

; now plot the local RM objects
file = '/home/yshen/Research/IDL/lib/Projects/reverberation_mapping/lc_data/sample_shen'
readcol,file,format='x,x,d,d,d', lag_rm,logL_rm, z_rm, m_rm
oplot, z_rm, logL_rm, psym = symcat(9,thick=3),color=cgcolor('red')
logmbh_rm=alog10(m_rm*1d6)
print, 'RM sample (N, zmed):', n_elements(z_rm), median(z_rm)

; plot the Bennert10 points
if keyword_set(bennert) then begin
  file='/data3/quasar/yshen/work/lags/prepspec/output/paper/Bennert10'
  readcol,file,format='x,x,d,x,x,x,x,x,x,x,x,d,x,d', zz, LL,logmbh
  logL=alog10(1d44*LL)
  oplot, zz, logL, psym=2,color=cgcolor('dark gray')
endif

; plot the HST cyc23 targets
if keyword_set(hst) then begin
  file='/data3/quasar/yshen/work/lags/prepspec/output/paper/hst23'
  readcol,file,format='x,x,d,x,x,x,x,x,x,d,d', z, logmbh_shen,logl_shen
  oplot, z, logL_shen, psym=symcat(16),color=cgcolor('blue'),symsize=1.5
  logmbh_shen=logmbh_shen+alog10(5.5)
endif else begin ; plot the reported detections in Shen++2015
  file='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/tables/lag_results'
  readcol,file,format='x,x,d,x,x,x,x,x,x,x,x,x,x,d',z,logL_shen
  oplot, z, logL_shen, psym=symcat(16),color=cgcolor('blue'),symsize=1.5
  print, 'SDSS-RM (N, zmed):', n_elements(z), median(z)
endelse

if keyword_Set(hst) then begin
  legend, 'this proposal', psym=symcat(16), pos=[0.68,0.85], box=0,/norm,color=cgcolor('blue'),$
  textcolor=cgcolor('blue')
endif else begin
  legend, 'this work', psym=symcat(16), pos=[0.68,0.88], box=0,/norm,color=cgcolor('blue'),$
  textcolor=cgcolor('blue')
endelse
legend, 'local RM AGN', box=0, psym = symcat(9,thick=3),color=cgcolor('red'),/norm, $
  pos=[0.68, 0.83],textcolor=cgcolor('red')
if keyword_set(bennert) then legend, 'Bennert10', box=0, psym = symcat(2),color=cgcolor('dark gray'),/norm, $
  pos=[0.68, 0.90],textcolor=cgcolor('dark gray')

endplot

if keyword_set(plot_mass) then begin
figfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/Mbh_z_dist.eps'
begplot,name=figfile,xsize=7,ysize=5,/encap,/color
plot,[0],[0],/nodata,xrange=[0,0.8],yrange=[6, 10],xtitle='Redshift', $
  ytitle=textoidl('logM_{BH} [M')+sunsymbol()+']',pos=[0.12,0.12,0.95,0.98],/xsty
oplot, z_rm,logmbh_rm,psym = symcat(9,thick=3),color=cgcolor('red')
oplot,zz,logmbh,psym=2,color=cgcolor('dark gray')
oplot,z,logmbh_shen,psym=symcat(16),color=cgcolor('blue'),symsize=1.5
legend, 'this proposal (RM BH mass)', psym=symcat(16), pos=[0.4,0.85], $
   box=0,/norm,color=cgcolor('blue'), textcolor=cgcolor('blue')
legend, 'local RM AGN (RM BH mass)', box=0, psym = symcat(9,thick=3),color=cgcolor('red'),/norm, $
  pos=[0.4, 0.95],textcolor=cgcolor('red')
legend, 'Bennert10 (SE BH mass)', box=0, psym = symcat(2),color=cgcolor('dark gray'),/norm, $
  pos=[0.4, 0.90],textcolor=cgcolor('dark gray')
endplot
endif

end
