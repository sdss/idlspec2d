; compare RM masses and SE masses for the 16 spectro-only detections

pro plot_Mrm_Mse

file='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/tables/lag_results'

readcol,file,format='x,x,x,x,a,x,x,x,x,x,d,d,d,x,x,x,x,x,d,d',line,Mrm, Err_lo, Err_hi, Mse, Err
Mrm=Mrm + alog10(5.5)

figfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/Mrm_Mse.eps'
pos=[0.15,0.15, 0.96,0.98]
range=[6.5,9]
begplot, name=figfile,/color,/encap, xsize=5, ysize=5,/cmyk

plot,[0],[0],/nodata,xrange=range, yrange=range,xtitle=textoidl('log M_{RM} (f=5.5) [M')+sunsymbol()+']', $
 ytitle=textoidl('log M_{SE,H\beta} [M')+sunsymbol()+']',pos=pos
oplot,range,range,line=2

ind=where(line eq 'H\beta')
oploterror, Mrm[ind], Mse[ind], Err_lo[ind], Err[ind], /lobar, psym=symcat(16)
oploterror, Mrm[ind], Mse[ind], Err_hi[ind], Err[ind], /hibar, psym=symcat(16)
ind=where(line eq 'MgII')
oploterror, Mrm[ind], Mse[ind], Err_lo[ind], Err[ind], /lobar, psym=symcat(9),color=cgcolor('red'),errcolor=cgcolor('red')
oploterror, Mrm[ind], Mse[ind], Err_hi[ind], Err[ind], /hibar, psym=symcat(9),color=cgcolor('red'),errcolor=cgcolor('red')

pos=[0.2,0.95]
legend, textoidl('H\beta'),/norm,pos=pos,box=0, psym=symcat(16)
pos[1]=pos[1]-0.05
legend, textoidl('MgII'),/norm,pos=pos,box=0, psym=symcat(9),color=cgcolor('red'),textcolor=cgcolor('red')
endplot

figfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/Mrm_Mse_diff.eps'
pos=[0.15,0.15, 0.96,0.98]
range=[-1.5,1.5]
begplot, name=figfile,/color,/encap, xsize=5, ysize=5,/cmyk
yrange=[6.5,9.2]
plot,[0],[0],/nodata,xrange=range, yrange=yrange,xtitle=textoidl('log (M_{SE,H\beta}/M_{RM})'), $
 ytitle=textoidl('log M_{RM} [M')+sunsymbol()+']',pos=pos,/ysty,/xsty
oplot, [0,0],yrange
oplot, [-0.5,-0.5],yrange,line=1
oplot, [0.5,0.5],yrange,line=1

ind=where(line eq 'H\beta')
oploterror, Mse[ind]-Mrm[ind], Mrm[ind], Err_hi[ind], Err_lo[ind], /lobar, psym=symcat(16)
oploterror, Mse[ind]-Mrm[ind], Mrm[ind], Err_lo[ind], Err_hi[ind], /hibar, psym=symcat(16)
ind=where(line eq 'MgII')
oploterror, Mse[ind]-Mrm[ind], Mrm[ind], Err_hi[ind], Err_lo[ind], /lobar, psym=symcat(9),color=cgcolor('red'),errcolor=cgcolor('red')
oploterror, Mse[ind]-Mrm[ind], Mrm[ind], Err_lo[ind], Err_hi[ind], /hibar, psym=symcat(9),color=cgcolor('red'),errcolor=cgcolor('red')
pos=[0.2,0.95]
legend, textoidl('H\beta'),/norm,pos=pos,box=0, psym=symcat(16)
pos[1]=pos[1]-0.05
legend, textoidl('MgII'),/norm,pos=pos,box=0, psym=symcat(9),color=cgcolor('red'),textcolor=cgcolor('red')

endplot

end
