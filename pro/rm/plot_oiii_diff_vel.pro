; plot the oiii velocity difference from different
; methods of measuring the oiii line

pro plot_oiii_diff_vel

   file='/data3/yshen/work/lineshifts/lineshift.fits'

   result=mrdfits(file,1)
   Rfe=dblarr(n_elements(result)) - 1.
   ind=where( (result.Hbeta_br)[3,*] gt 0)
   Rfe[ind] = (result.REW_FE_4434_4684)[ind] / (result.Hbeta_br)[3,ind]  

   tags=tag_names(result)
   cspeed=2.9979246d5
   wave=5008.24D

   figfile='/data3/yshen/work/lineshifts/oiii_diff_vel.ps'
   begplot,name=figfile,/color,/landscape
   npanel = 4
   charsize=1.2 & len=0.04
   plot_layout, npanel, xypos=xypos, omargin=[0.1, -0.02, 0.98, 0.96], pmargin=[0.15,0.25]

   ; now try v_full_peak
   i1 = where(tags eq 'OIII5007C_ERR') ; err in the line
   i2 = where(tags eq 'OIII5007_ERR') ; err in the line
   j1 = where(tags eq 'OIII5007C' )
   j2 = where(tags eq 'OIII5007' )

   ind = where( (result.(i1))[0,*] gt 0 and (result.(i2))[0,*] gt 0  $
          and (result.(i1))[2,*] gt 0 and (result.(i1))[2,*] le 1.D/3./alog(10.D) $ ; the line flux is 3sig detection
          and (result.(i2))[2,*] gt 0 and (result.(i2))[2,*] le 1.D/3./alog(10.D) )
   vel1 = ( (result[ind].(j1))[0,*] - wave)/wave*cspeed
   vel1_err = (result[ind].(i1))[0,*]/wave*cspeed
   vel2 = ( (result[ind].(j2))[0,*] - wave)/wave*cspeed
   vel2_err = (result[ind].(i2))[0,*]/wave*cspeed
   logL=result[ind].logL5100
   yrange=[-60.,60.]

   indd = where(vel1_err gt 0 and vel1_err lt 500. and vel2_err gt 0 and vel2_err lt 500.)
   med_err=median(sqrt(vel1_err[indd]^2 + vel2_err[indd]^2))
   xdata=logL[indd] & ydata=(vel2 - vel1)[indd]
   plot, xdata, ydata, psym = 2,xtitle=textoidl('logL5100 [erg/s]'), $
      title=textoidl('[OIII]_a-[OIII]_c' +' [km/s]'), pos=xypos[*,0], yrange=yrange,symsize=0.2, $
      xrange=xrange,/xsty, charsize=charsize, noerase=0, xticklen=len, yticklen=len
   vmed=median(ydata)
   print, vmed
   oplot, [41.,46.], replicate(vmed,2), color=cgcolor('red')

   indd = where(vel1_err gt 0 and vel1_err lt 500. and vel2_err gt 0 and vel2_err lt 500. and Rfe[ind] gt -0.1)
   med_err=median(sqrt(vel1_err[indd]^2 + vel2_err[indd]^2))
   xdata=Rfe[ind[indd]] & ydata=(vel2 - vel1)[indd]
   plot, xdata, ydata, psym = 2,xtitle=textoidl('R_{FeII}'), $
      title=textoidl('[OIII]_a-[OIII]_c' +' [km/s]'), pos=xypos[*,1], yrange=yrange,symsize=0.2, $
      xrange=[-0.1,3],/xsty, charsize=charsize, noerase=1, xticklen=len, yticklen=len
   vmed=median(ydata)
   print, vmed
   oplot, [-1.,5], replicate(vmed,2), color=cgcolor('red')

   ; now try v_cen_50per
   i1 = where(tags eq 'OIII5007C_ERR') ; err in the line
   i2 = where(tags eq 'OIII5007_ERR') ; err in the line
   j1 = where(tags eq 'OIII5007C' )
   j2 = where(tags eq 'OIII5007' )

   ind = where( (result.(i1))[0,*] gt 0 and (result.(i2))[4,*] gt 0  $
          and (result.(i1))[2,*] gt 0 and (result.(i1))[2,*] le 1.D/3./alog(10.D) $ ; the line flux is 3sig detection
          and (result.(i2))[2,*] gt 0 and (result.(i2))[2,*] le 1.D/3./alog(10.D) )
   vel1 = ( (result[ind].(j1))[0,*] - wave)/wave*cspeed
   vel1_err = (result[ind].(i1))[0,*]/wave*cspeed
   vel2 = ( (result[ind].(j2))[4,*] - wave)/wave*cspeed
   vel2_err = (result[ind].(i2))[4,*]/wave*cspeed
   logL=result[ind].logL5100
   yrange=[-60.,60.]

   indd = where(vel1_err gt 0 and vel1_err lt 500. and vel2_err gt 0 and vel2_err lt 500.)
   med_err=median(sqrt(vel1_err[indd]^2 + vel2_err[indd]^2))
   xdata=logL[indd] & ydata=(vel2 - vel1)[indd]
   plot, xdata, ydata, psym = 2,xtitle=textoidl('logL5100 [erg/s]'), $
      title=textoidl('[OIII]_{cent}-[OIII]_c' +' [km/s]'), pos=xypos[*,2], yrange=yrange,symsize=0.2, $
      xrange=xrange,/xsty, charsize=charsize, noerase=1, xticklen=len, yticklen=len
   vmed=median(ydata)
   print, vmed
   oplot, [41.,46.], replicate(vmed,2), color=cgcolor('red')

   indd = where(vel1_err gt 0 and vel1_err lt 500. and vel2_err gt 0 and vel2_err lt 500. and Rfe[ind] gt -0.1)
   med_err=median(sqrt(vel1_err[indd]^2 + vel2_err[indd]^2))
   xdata=Rfe[ind[indd]] & ydata=(vel2 - vel1)[indd]
   plot, xdata, ydata, psym = 2,xtitle=textoidl('R_{FeII}'), $
      title=textoidl('[OIII]_{cent}-[OIII]_c' +' [km/s]'), pos=xypos[*,3], yrange=yrange,symsize=0.2, $
      xrange=[-0.1,3],/xsty, charsize=charsize, noerase=1, xticklen=len, yticklen=len
   vmed=median(ydata)
   print, vmed
   oplot, [-1.,5], replicate(vmed,2), color=cgcolor('red')

   endplot
   cgfixps, figfile

end

