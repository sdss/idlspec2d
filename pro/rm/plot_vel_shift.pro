; plot the velocity shifts among different lines
 
pro plot_vel_shift, hist=hist, fitfile=fitfile, figfile=figfile

if ~keyword_set(fitfile) then begin
 ;file='/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/qso_prop-0000-56837_lineshift.fits'
 ;file='/data3/yshen/work/lineshifts/lineshift.fits' 
  file='/data3/yshen/work/lineshifts/lineshift_outlierfix.fits' ; this fixes the one [oii] vel outlier
endif else file=fitfile

result=mrdfits(file,1)
cspeed=2.9979246d5
tags=tag_names(result)

; make a structure to store all the peak vel measurements
nobj=n_elements(result)
output={RMID:0L, z:0.D, logL1350:[0.D, -1.D], logL1700:[0.D, -1.D],logL3000:[0.D, -1.D],logL5100:[0.D, -1.D], $
  Hbeta_br:[0.D, -1.D], OIII5008c:[0.D, -1.D], OIII5008a:[0.D, -1.D], CaII3934:[0.D, -1.D], OII3728:[0.D, -1.D], $
  NeV3426:[0.D, -1.D], MgII:[0.D, -1.D], CIII:[0.D, -1.D], CIIIa:[0.D, -1.D], HeII1640:[0.D, -1.D], CIV:[0.D, -1.D], $
  SiIV:[0.D, -1.D]}
output=replicate(output, nobj)
line_out=['Hbeta_br', 'OIII5007c', 'OIII5007', 'CaII3934', 'OII3728', 'NeV3426', 'MgII', 'CIII_br', 'CIII', 'HeII1640', 'CIV_br', 'SIIV_OIV']
line_out=strupcase(line_out)
wave_out=[4862.68, 5008.24, 5008.24, 3934.78, 3728.48, 3426.84, 2798.75, 1908.73, 1908.73, 1640.42, 1549.06, (1396.76 + 1402.06)*0.5]
; now compile the peak vel shifts
for ii=0L, n_elements(line_out) - 1 do begin
   i1 = where(tags eq line_out[ii] + '_ERR') ; err in the line
   j1 = where(tags eq line_out[ii] )

   ind = where( (result.(i1))[0,*] gt 0 and (result.(i1))[2,*] gt 0 and (result.(i1))[2,*] le 1.D/3./alog(10.D)) 
   ;the line flux is 3sig detection

   wave1=wave_out[ii]
   ;vel1=0.*(result.(j1))[0,*] & vel1_err=0.*(result.(j1))[0,*] - 1.D
   vel1 = ( (result[ind].(j1))[0,*] - wave1)/wave1*cspeed
   vel1_err = (result[ind].(i1))[0,*]/wave1*cspeed

   print, (tag_names(output))[ii+6], " ",  line_out[ii], " ", wave_out[ii]
   output[ind].(ii+6) = [vel1, vel1_err]

   ;message, 'stop'

endfor
output.RMID=result.RM_ID & output.z=result.z
for ii=0, nobj - 1 do begin
  output[ii].logL1350=[result[ii].logL1350, result[ii].logL1350_err]
  output[ii].logL1700=[result[ii].logL1700, result[ii].logL1700_err]
  output[ii].logL3000=[result[ii].logL3000, result[ii].logL3000_err]
  output[ii].logL5100=[result[ii].logL5100, result[ii].logL5100_err]
endfor
outfile='/data3/yshen/work/lineshifts/peak_vel_cat.fits'
mwrfits, output, outfile, /create

if ~keyword_set(figfile) then figfile='/data3/yshen/work/lineshifts/vel_shift.eps'
;'/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/vel_shift.ps'
begplot,name=figfile,/color,/landscape

charsize=0.8

; if n_elements(hist) eq 0 then hist=1
histpeak=0.3 & histbin=0.1 & len=0.04

line=['Hbeta_br', 'OIII5007c', 'CaII3934', 'OII3728', 'NeV3426', 'MgII', 'CIII_br', 'CIII', 'HeII1640', 'CIV_br', 'OIII5007', 'SIIV_OIV', 'LYA_BR','SIIII1892','ALIII1857']

line=strupcase(line)
linetag=textoidl(['H\beta_{br}', '[OIII]_c', 'CaII', '[OII]', '[NeV]', 'MgII', 'CIII]', 'CIII]_a', 'HeII', 'CIV', '[OIII]_a', 'SiIV', 'Ly\alpha','SIIII1892','ALIII1857'])
wave=[4862.68, 5008.24, 3934.78, 3728.48, 3426.84, 2798.75, 1908.73, 1908.73, 1640.42, 1549.06, 5008.24, (1396.76 + 1402.06)*0.5, 1215.67]
;(1396.76 + 1402.06)*0.5

; print how many objects were fit for each line
file1='/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/qso_prop-0000-56837_lineshift.fits'
result_tt=mrdfits(file1,1)
tags_tt=tag_names(result_tt)
line_tt=['SIIII1892', 'ALIII1857']
for ii=0L, n_elements(line_tt) - 1 do begin
  i1=where(tags_tt eq line_tt[ii] + '_ERR')
  j1=where(tags_tt eq line_tt[ii])
  ind_tmp=where( (result_tt.(i1))[0,*] gt 0 and (result_tt.(j1))[0,*] gt 0, ntmp)
  print, line_tt[ii], ntmp
endfor

pair = [ [2, 0] , $ ; broad Hbeta-CaII
         [2, 1] , $ ; OIII-CaII
         [2,10] , $ ; OIII_a - CaII
         [2, 3] , $ ; OII-CaII
         [2, 4] , $ ; NeV-OII
         [3, 5] , $ ; MgII-OII
         [3, 6] , $ ; CIII-OII
         [3, 7] , $ ; CIII_all-OII
         [3, 8] , $ ; HeII-OII
         [8, 9] , $ ; CIV-HeII
         [5, 9] , $ ; CIV-MgII
         [5 ,11] ] ; $ ; SiIV/OIV-MgII
;         [5, 12] ]  ; Lya-MgII
lum=['logL5100', 'logL5100', 'logL5100', 'logL3000','logL3000', 'logL3000', 'logL1700', 'logL1700', 'logL1700', 'logL1700', 'logL1700', 'logL1700','logL1700']
lum_ref=[44.,44.,44.,44.5,44.5,44.5,45.,45.,45.,45.,45.,45.]
; lum=strupcase(lum)

; set up plot layout
npanel = (size(pair))[2]
plot_layout, npanel, xypos=xypos, omargin=[0.06, 0.01, 0.98, 0.96], pmargin=[0.05,0.12]
for i=0L, npanel - 1 do begin

   if i eq 0 then noerase=0 else noerase=1
   pos = xypos[*, i]
   i1 = where(tags eq line[pair[0,i]] + '_ERR') ; err in the line
   i2 = where(tags eq line[pair[1,i]] + '_ERR') ; err in the line
   j1 = where(tags eq line[pair[0,i]] )
   j2 = where(tags eq line[pair[1,i]] )

   ind = where( (result.(i1))[0,*] gt 0 and (result.(i2))[0,*] gt 0  $
          and (result.(i1))[2,*] gt 0 and (result.(i1))[2,*] le 1.D/3./alog(10.D) $ ; the line flux is 3sig detection
          and (result.(i2))[2,*] gt 0 and (result.(i2))[2,*] le 1.D/3./alog(10.D) )
   jj=where( tags eq strupcase(lum[i]))
   if lum[i] eq 'logL5100' then xrange=[43, 45]
   if lum[i] eq 'logL3000' then xrange=[43.5,46]
   if lum[i] eq 'logL1700' then xrange=[44,46.5]

   logL=result[ind].(jj) 
   wave1=wave[pair[0,i]] & wave2=wave[pair[1,i]]
   vel1 = ( (result[ind].(j1))[0,*] - wave1)/wave1*cspeed
   vel1_err = (result[ind].(i1))[0,*]/wave1*cspeed
   vel2 = ( (result[ind].(j2))[0,*] - wave2)/wave2*cspeed
   vel2_err = (result[ind].(i2))[0,*]/wave2*cspeed
   ; default vel error tolerance is 500 km/s 
   indd = where(vel1_err gt 0 and vel1_err lt 500. and vel2_err gt 0 and vel2_err lt 500.)
   ;indd = where(vel1_err gt 0 and vel1_err lt 50. and vel2_err gt 0 and vel2_err lt 50.)
   tmp_arr = (vel2 - vel1)[where(vel1_err gt 0 and vel1_err lt 500. and vel2_err gt 0 and vel2_err lt 500.)]
   tmp_med = median(tmp_arr)
   tmp_std = stddev(tmp_arr)
   tmp_std2 = 0.5*(quantile_1d(0.84, tmp_arr) - quantile_1d(0.16, tmp_arr))
   ; print, tmp_med, tmp_std
   ;indd = where(vel1_err gt 0 and vel1_err lt 500. and vel2_err gt 0 and vel2_err lt 500. and abs(vel2 - vel1 - tmp_med) le 3.*tmp_std2 )
   med_err=median(sqrt(vel1_err[indd]^2 + vel2_err[indd]^2))
   err16 = quantile_1d(0.16, sqrt(vel1_err[indd]^2 + vel2_err[indd]^2))
   err84 = quantile_1d(0.84, sqrt(vel1_err[indd]^2 + vel2_err[indd]^2))
   xdata=logL[indd] & ydata=(vel2 - vel1)[indd]/1000.
   per16=quantile_1d(0.16, ydata) & per84=quantile_1d(0.84, ydata)
   per1_x=quantile_1d(0.01, xdata) & per99_x=quantile_1d(0.99, xdata)
   print, 'per1 logl=', per1_x

   vmed=median(ydata) ; vscat=sqrt( (stddev(ydata))^2 - (med_err/1000.)^2)
   ; calc vscat in a crude way
   vscat=0.5*(per84 - per16)
   ;print, vscat
   vscat=sqrt(vscat^2 - (med_err/1000.D)^2 )

   nobj=n_elements(ydata)

   ; here we remove a few objects with logl=0
   ; note that Fig. 1 in the paper was incorrectly produced w/o this cut, but the results are nevertheless consistent
   ind_nonzero = where(xdata gt 0)
   xdata = xdata[ind_nonzero] & ydata = ydata[ind_nonzero]

   nobj2=n_elements(ydata)

   ind_tmp=where(ydata ge per16 and ydata le per84)
   ;yrange=round(max( abs(ydata[ind_tmp]) ) * [-1,1]*10)/10.
   yrange=[-1., 1.]

   title=textoidl(linetag[pair[1,i]] + '-' +linetag[pair[0,i]])
   plot, xdata, ydata, psym = 2,xtitle=textoidl(lum[i] + ' [erg/s]'), $
      title=title, pos=pos, yrange=yrange,symsize=0.2, $
      xrange=xrange,/xsty, charsize=charsize, noerase=noerase, xticklen=len, yticklen=len
   oplot, xrange,[0,0], line=2

   ;plot a moving average
   ; in the paper, Fig. 1 was produced with nmin=2
   moving_average, xdata, ydata, xbin=0.35, xmin=per1_x, xmax=per99_x, nmin=5 $
      , xarr, yarr, yerr ;, /boots
   ind_good=where(yerr gt 0)
   oploterror, xarr[ind_good], yarr[ind_good], yerr[ind_good],psym=symcat(3),color=cgcolor('green'),errcolor=cgcolor('green')
   imax = n_elements(xdata)

   ; do a linear regression with all the data points
   sixlin, xdata - lum_ref[i], ydata, aa, siga, bb, sigb
   print, title + ' N=', imax
   print, 'med err(v12)=', med_err
   ; Spearman's test
   spearman=r_correlate(xdata, ydata)
   print, 'Spearman r and p:', spearman
   ;if title eq '[NeV]-CaII' then begin
      nsample = 1000L
      boots_ind = boot_indices(imax, nsample=nsample)
      r_arr = dblarr(nsample) & p_arr = dblarr(nsample)
      for iboot=0, nsample-1 do begin
         tmp = r_correlate(xdata[boots_ind[iboot,*]], ydata[boots_ind[iboot,*]])
         r_arr[iboot] = tmp[0] & p_arr[iboot] = tmp[1]
      endfor
      print, 'Bootstrap Spearman r and p'
      print, median(r_arr), median(p_arr)
      print, quantile_1d(0.16, r_arr), quantile_1d(0.16,p_arr)
      print, quantile_1d(0.84, r_arr), quantile_1d(0.84,p_arr)
   ;endif


   ;oplot, xrange, replicate(vmed,2), color=cgcolor('red')
   ; plot a regression fit
   regfit=mpfitexpr('P(0) + P(1)*x', xarr[ind_good]-lum_ref[i], yarr[ind_good], yerr[ind_good], [0., 0.], perror=perror,/quiet)
   ;print, title, lum_ref[i]
   print, 'linreg (a):', aa[0]*1000., '+-', siga[0]*1000.
   print, '(b):', bb[0]*1000., '+-', sigb[0]*1000.
   print, 'binned fit (a):', regfit[0]*1000., '+-',perror[0]*1000.
   print, '(b):', regfit[1]*1000., '+-',perror[1]*1000.
   print, '-------------------------------'

   oplot, xrange, regfit[0] + regfit[1]*(xrange - lum_ref[i]), color=cgcolor('red')

   yoff0 = med_err/1000.
   yoff1 = (err16 - med_err)/1000. & yoff2 = (err84 - med_err)/1000.
   oploterror, xrange[1]-0.5, yrange[1]*0.85 - yoff0, med_err/1000., psym=3, color=cgcolor('red'), /nohat
   oploterror, xrange[1]-0.6, yrange[1]*0.85 -yoff0 - yoff1, err16/1000., psym=3, color=cgcolor('red'), /nohat
   oploterror, xrange[1]-0.4, yrange[1]*0.85 -yoff0 - yoff2, err84/1000., psym=3, color=cgcolor('red'), /nohat
   xyouts, xrange[0]+0.2, yrange[0]*0.85, 'N='+string(nobj,format='(i0)'), charsize=charsize
   xyouts, xrange[0]+0.2, yrange[1]*0.75, string(vmed*1d3,format='(i0)')+ '/' $
         + string(vscat*1d3,format='(i0)') + ' km/s',charsize=charsize
   if keyword_set(hist) then begin
      plothist_old, ydata, bin=histbin, xhist,yhist, /rotate,color=cgcolor('cyan'), /noerase, $
      xsty=5,ysty=5,pos=pos,yrange=yrange, xrange=[0,2], peak=histpeak, $
      xtickname=replicate(' ',10L), ytickname=replicate(' ',10L), thick=5
   endif

endfor
xyouts, 0.5, 0.01, textoidl('log (Continuum Luminosity) [erg s^{-1}]'), /norm, align=0.5
xyouts, 0.025, 0.5, textoidl('Velocity [1000 km s^{-1}]'), /norm, align=0.5, orient=90
endplot
cgfixps, figfile

figfile='/data3/yshen/work/lineshifts/vel_shift_hist.eps'
begplot,name=figfile,/color,/landscape
plot_layout, npanel, xypos=xypos, omargin=[0.05, 0.0, 0.99, 0.96], pmargin=[0.065,0.1]
for i=0L, npanel - 1 do begin

   if i eq 0 then noerase=0 else noerase=1
   pos = xypos[*, i]
   i1 = where(tags eq line[pair[0,i]] + '_ERR') ; err in the line
   i2 = where(tags eq line[pair[1,i]] + '_ERR') ; err in the line
   j1 = where(tags eq line[pair[0,i]] )
   j2 = where(tags eq line[pair[1,i]] )

   ind = where( (result.(i1))[0,*] gt 0 and (result.(i2))[0,*] gt 0  $
          and (result.(i1))[2,*] gt 0 and (result.(i1))[2,*] le 1.D/3./alog(10.D) $ ; the line flux is 3sig detection
          and (result.(i2))[2,*] gt 0 and (result.(i2))[2,*] le 1.D/3./alog(10.D) )
   jj=where( tags eq strupcase(lum[i]))

   logL=result[ind].(jj)
   wave1=wave[pair[0,i]] & wave2=wave[pair[1,i]]
   vel1 = ( (result[ind].(j1))[0,*] - wave1)/wave1*cspeed
   vel1_err = (result[ind].(i1))[0,*]/wave1*cspeed
   vel2 = ( (result[ind].(j2))[0,*] - wave2)/wave2*cspeed
   vel2_err = (result[ind].(i2))[0,*]/wave2*cspeed
   indd = where(vel1_err gt 0 and vel1_err lt 500. and vel2_err gt 0 and vel2_err lt 500.)
   med_err=median(sqrt(vel1_err[indd]^2 + vel2_err[indd]^2))
   xdata=logL[indd] & ydata=(vel2 - vel1)[indd]/1000.

   nobj=n_elements(ydata)

   per16=quantile_1d(0.16, ydata) & per84=quantile_1d(0.84, ydata)
   per5_y=quantile_1d(0.05, ydata) & per95_y=quantile_1d(0.95, ydata)
   ind_tmp=where(ydata ge per16 and ydata le per84)

   binsz=stddev(ydata)/5.
   xrange=binsz*5*3.*[-1,1]

   plothist, ydata, bin=binsz, ytitle=textoidl('N_{obj}'), xrange=xrange, $
      xtitle=textoidl(linetag[pair[1,i]] + '-' +linetag[pair[0,i]] +' [1000 km/s]'), pos=pos, $
       charsize=charsize, noerase=noerase, xticklen=len, yticklen=len, xhist,yhist
   ;oplot, [0,0], line=2
   fit1=mpfitfun('gauss1', xhist,yhist,yhist,weights=1.D, [median(ydata), binsz*5., 10.],/quiet)
   xgrid=0.02*findgen(301L) - 3.
   oplot, xgrid, gauss1(xgrid, fit1), color=cgcolor('red')

   ; plot corrected L-dependence
   if i eq 8 or i eq 10 or i eq 11 then begin
      if i eq 8 then begin
         aa=-231. & bb=-282.
      endif
      if i eq 10 then begin
         aa=-242. & bb=-438.
      endif
      if i eq 11 then begin
         aa=-123. & bb=-345.
      endif
      aa=aa/1000.D & bb=bb/1000.D
      ydata_corr = ydata - (aa + bb*(xdata - lum_ref[i]))
      plothist, ydata_corr, bin=binsz, /over, line=2,color=cgcolor('cyan'), xhist2,yhist2
      fit2=mpfitfun('gauss1', xhist2,yhist2,yhist2,weights=1.D, [median(ydata_corr), binsz*5., 10.],/quiet)
      oplot, xgrid, gauss1(xgrid, fit2), color=cgcolor('blue'),line=2
      int_scat2=sqrt( (fit2[1]*1d3)^2 - med_err^2 )
      xyouts, pos[0]+0.015, pos[3]-0.06, string(fit2[0]*1d3,format='(i0)')+ '/' $
         + string(int_scat2, format='(i0)') + ' km/s',charsize=charsize, /norm,color=cgcolor('blue')
   endif

   int_scat=sqrt( (fit1[1]*1d3)^2 - med_err^2 )
   xyouts, pos[0]+0.015, pos[3]-0.03, string(fit1[0]*1d3,format='(i0)')+ '/' $
         + string(int_scat, format='(i0)') + ' km/s',charsize=charsize, /norm
   xyouts, pos[2]-0.06, pos[3]-0.03, 'N='+string(nobj,format='(i0)'), charsize=charsize,/norm

   vmed=median(ydata) ; vscat=sqrt( (stddev(ydata))^2 - (med_err/1000.)^2)
   ; calc vscat in a crude way
   vscat=0.5*(per84 - per16)
   vscat=sqrt(vscat^2 - (med_err/1000.D)^2 )
endfor
endplot
cgfixps, figfile


end
