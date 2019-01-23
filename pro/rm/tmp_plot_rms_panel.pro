; Plot quick-check plots for Jon

pro tmp_plot_rms_panel, result=result, epoch_ID=epoch_ID

  if n_elements(result) eq 0 then $
  rm_get_rms_prop,result=result, epoch_ID=epoch_ID


  figfile=getenv('IDLRM_DIR') + '/misc/rms_var_check.eps'
  begplot, name=figfile,/cmyk,/color,xsize=10,ysize=8,/encap

  ind=where(result.hbeta_rms gt 0)
  plot, [0],[0],xrange=[0.01,10],yrange=[0.01,10],/xlog,/ylog,$
   xtitle='Median Line Error', ytitle='Line RMS Variation', $
   pos=[0.12, 0.58, 0.48, 0.98],xtickname=['0.01','0.1','1.0','10'],ytickname=['0.01','0.1','1.0','10']
  oplot, result[ind].hbeta_mederr,result[ind].hbeta_rms,psym=symcat(9),color=cgcolor('red')
  ind=where(result.mgii_rms gt 0)
  oplot, result[ind].mgii_mederr,result[ind].mgii_rms,psym=symcat(9),color=cgcolor('blue')
  ind=where(result.OIII_rms gt 0)
  oplot, result[ind].oiii_mederr,result[ind].oiii_rms,psym=1
  oplot, [0.01,10],[0.01,10],line=2
  xyouts,0.015,5,textoidl('H\beta'),color=cgcolor('red')
  xyouts,0.015,3,'MgII',color=fsc_color('blue')
  xyouts,1.0,0.017,'[OIII]5007'
  if keyword_set(epoch_id) then begin
     epstr=string(epoch_id[0]+1,format='(i0)')
     for i=1L,n_elements(epoch_id)-1 do epstr=epstr+','+string(epoch_id[i]+1,format='(i0)')
    xyouts, 0.05,5,'ep:'+epstr
  endif
  plot, [0],[0],xrange=[0.01,10],yrange=[0.01,10],/xlog,/ylog,$
   xtitle='[OIII] RMS Variation', ytitle='Line RMS Variation', $
   pos=[0.58, 0.58, 0.94, 0.98],/noerase,xtickname=['0.01','0.1','1.0','10'],ytickname=['0.01','0.1','1.0','10']
  ind=where(result.hbeta_rms gt 0 and result.oiii_rms gt 0)
  oplot, result[ind].oiii_rms,result[ind].hbeta_rms,psym=symcat(9),color=cgcolor('red')
  ind=where(result.mgii_rms gt 0 and result.oiii_rms gt 0)
  oplot, result[ind].oiii_rms,result[ind].mgii_rms,psym=symcat(9),color=cgcolor('blue')
  oplot, [0.01,10],[0.01,10],line=2
  xyouts,0.015,5,textoidl('H\beta'),color=cgcolor('red')
  xyouts,0.015,3,'MgII',color=fsc_color('blue')

  plot, [0],[0],xrange=[0.01,10],yrange=[0.01,10],/xlog,/ylog,$
   xtitle='L3000 RMS Variation', ytitle='Line RMS Variation', $
   pos=[0.12, 0.09, 0.48, 0.49],/noerase,xtickname=['0.01','0.1','1.0','10'],ytickname=['0.01','0.1','1.0','10']
  ind=where(result.hbeta_rms gt 0 and result.L3000_rms gt 0)
  oplot, result[ind].L3000_rms,result[ind].hbeta_rms,psym=symcat(9),color=cgcolor('red')
  ind=where(result.mgii_rms gt 0 and result.L3000_rms gt 0)
  oplot, result[ind].L3000_rms,result[ind].mgii_rms,psym=symcat(9),color=cgcolor('blue')
  oplot, [0.01,10],[0.01,10],line=2
  xyouts,0.015,5,textoidl('H\beta'),color=cgcolor('red')
  xyouts,0.015,3,'MgII',color=fsc_color('blue')

  plot, [0],[0],xrange=[0.002,10],yrange=[0.01,10],/xlog,/ylog,/xsty, $
   xtitle='Median L3000 Error', ytitle='L3000 RMS Variation', $
   pos=[0.58, 0.09, 0.94, 0.49],/noerase,xtickname=['0.01','0.1','1.0','10'],ytickname=['0.01','0.1','1.0','10']
  ind=where(result.L3000_rms gt 0)
  oplot, result[ind].L3000_mederr,result[ind].L3000_rms,psym=symcat(9)
  oplot, [0.01,10],[0.01,10],line=2
  oplot, [0.01,10],[0.04,0.04],line=1,color=cgcolor('cyan'),thick=8
  xyouts, 0.04,0.025,'4% spectrophotometry',color=cgcolor('cyan')


  endplot

end
