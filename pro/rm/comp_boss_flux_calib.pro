;compare the flux calibration using 70 SDSS-RM std stars
; during 2014 (32 epochs) for v5_7_1, v5_10_10 and my custom version

pro comp_boss_flux_calib

rmid = [850:919]

epoch = [1:32]

; my custom flux calib
rm_calc_rms_spec,rmid, epoch=epoch,result=result_all, run2d='v5_7_1',run1d='v5_7_1',calibdir='wh_skysub/'
wave=result_all[0].wave
fracrms = result_all.rmsflux/result_all.medflux
; BOSS pipeline v5_7_1
rm_calc_rms_spec,rmid, epoch=epoch,result=result_all1, run2d='v5_7_1',run1d='v5_7_1',calibdir='/'
wave1=result_all1[0].wave
fracrms1 = result_all1.rmsflux/result_all1.medflux
; eBOSS pipeline v5_10_10
rm_calc_rms_spec,rmid, epoch=epoch,result=result_all2, run2d='v5_10_10',run1d='v5_10_10',calibdir='/'
wave2=result_all2[0].wave
fracrms2 = result_all2.rmsflux/result_all2.medflux

title='70 SDSS-RM flux standard stars in 2014 [56660-56837]; 32 epochs'
plot, wave,median(fracrms,dim=2), xrange=[3d3,1.2d4],/xsty,yrange=[0,0.2], xtitle='Obs Wavelength [A]', ytitle='RMS/AVG (median over 70 stars)', title=title
oplot, wave1, median(fracrms1,dim=2), color=cgcolor('cyan')
oplot, wave2, median(fracrms2,dim=2), color=cgcolor('red')

items=['v5_10_10 (eboss)', 'v5_7_1 (boss)', 'sdss-rm custom']
colors=cgcolor(['red', 'cyan', 'white'])

legend, items, color=colors, box=0,pos=[0.6,0.9],/norm,textcolor=colors


end
