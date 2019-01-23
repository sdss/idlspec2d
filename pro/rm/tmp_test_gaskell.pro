; perform some tests to dispute Gaskell's comments

pro tmp_test_gaskell


file='/data3/yshen/work/lineshifts/peak_vel_cat.fits'
result = mrdfits(file,1)


logl3000=(result.logl3000)[0,*]
logl3000_err = (result.logl3000)[1,*]

voii = (result.OII3728)[0,*]
voii_err = (result.oii3728)[1,*]
vcaii = (result.CAII3934)[0,*]
vcaii_err = (result.CAII3934)[1,*]

goto, test2

maxerr = 1d6
ind = where(logl3000_err gt 0 and voii_err gt 0 and voii_err lt maxerr and vcaii_err gt 0 and vcaii_err lt maxerr)
toterr = sqrt( voii_err[ind]^2 + vcaii_err[ind]^2  )
mederr_oii = median(voii_err[ind])
print, mederr_oii
mederr_caii = median(vcaii_err[ind])
print, mederr_caii
mederr = median(sqrt(voii_err[ind]^2 + vcaii_err[ind]^2))
mederr_oii_all = median(voii_err)
mederr_caii_all = median(vcaii_err)
mederr_all = sqrt(mederr_oii_all^2 + mederr_caii_all^2)


figfile='/data3/yshen/work/lineshifts/oii_caii.ps'
begplot, name=figfile, /color,/landscape


plot, logl3000[ind], voii[ind] - vcaii[ind], psym=symcat(2),symsize=0.5, $
  xrange=[43.0, 46.5], yrange=[-1500, 2500], /xsty, /ysty, $
  xtitle = textoidl('Log \lambdaL_{\lambda,3000} '), ytitle = textoidl('V_{[OII]}-V_{CaII} (km s^{-1})')
; oplot, [46.4,46.4], [1500, 1500.+mederr], color=cgcolor('red'), thick=5


nnn = n_elements(ind)
toterr2 = sqrt(toterr^2 + 46.^2)
mock = randomn(seed, nnn)*toterr + 1800.
mock2 = randomn(seed2, nnn)*toterr2 + 1800.
;oplot, logl3000[ind], mock, psym=symcat(2),symsize=0.5,color=cgcolor('red')
oplot, logl3000[ind], mock2, psym=symcat(2),symsize=0.5,color=cgcolor('green')


endplot
cgfixps, figfile


test2: 
; OIII - OII
vo3 = (result.OIII5008C)[0,*]
vo3_err = (result.OIII5008C)[1,*]
logl5100 = (result.logl5100)[0,*]
logl5100_err = (result.logl5100)[1,*]


maxerr = 500.  
ind = where(logl5100_err gt 0 and voii_err gt 0 and voii_err lt maxerr and vo3_err gt 0 and vo3_err lt maxerr and $
  abs(vo3 - voii) lt 500  )
imax = n_elements(ind)

figfile = '/data3/yshen/work/lineshifts/oiii_oii.ps'
begplot, name=figfile, /color,/landscape

xrange=[43,46.5]
plot, logl5100[ind], vo3[ind] - voii[ind], xrange=xrange, yrange=[-500,500],/xsty,/ysty,$
  xtitle= textoidl('Log \lambdaL_{\lambda,5100} '), ytitle = textoidl('V_{[OIII]c}-V_{[OII]} (km s^{-1})'), $
  psym=symcat(2),symsize=0.5
oplot, [43,46.5], [0,0]

xdata = logl5100[ind] & ydata = vo3[ind] - voii[ind]  & ydata_err = sqrt( vo3_err[ind]^2 + voii_err[ind]^2 )
lum_ref = 44.

per1_x=quantile_1d(0.01, xdata) & per99_x=quantile_1d(0.99, xdata)
moving_average, xdata, ydata, xbin=0.35, xmin=per1_x, xmax=per99_x, nmin=2 $
      , xarr, yarr, yerr ;, /boots
oploterror, xarr, yarr, yerr,psym=symcat(3),color=cgcolor('green'),errcolor=cgcolor('green')
ind_good=where(yerr gt 0)
regfit=mpfitexpr('P(0) + P(1)*x', xarr[ind_good]-lum_ref, yarr[ind_good], yerr[ind_good], [0., 0.], perror=perror,/quiet)
xrange1=[43., 45.7]
oplot, xrange1, regfit[0] + regfit[1]*(xrange1 - lum_ref), color=cgcolor('green')


sixlin, xdata - lum_ref, ydata, aa, siga, bb, sigb
spearman=r_correlate(xdata, ydata)
print, 'Spearman r and p:', spearman
nsample = 5000L
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

print, 'linreg (a):', aa[0], '+-', siga[0]
print, '(b):', bb[0], '+-', sigb[0]
oplot, xrange1, (xrange1 -lum_ref)*bb[0] + aa[0], color=cgcolor('red')

print, 'binned fit (a):', regfit[0], '+-',perror[0]
print, '(b):', regfit[1], '+-',perror[1]
print, '-------------------------------'
xyouts, 44, 450, '    16%     84%'
xyouts, 43.2, 400, 'Spearman r=' + string(spearman[0], format='(f0.2)') + '  ' + string(quantile_1d(0.16, r_arr), format='(f0.2)') + '  ' + string(quantile_1d(0.84, r_arr), format='(f0.2)')
xyouts, 43.2, 350, 'Spearman p=' + string(spearman[1], format='(e7.1)') + '  ' + string(quantile_1d(0.16, p_arr), format='(e7.1)') + '  ' + string(quantile_1d(0.84, p_arr), format='(e7.1)')


xyouts, 45.5, 400, textoidl('\sigma_{V1}, \sigma_{V_2}<') + string(maxerr,format='(i0)') + textoidl('km s^{-1}')
xyouts, 43.2, -400, 'linreg (Y|X) b='+string(bb[0],format='(i0)')+textoidl('\pm')+string(sigb[0],format='(i0)')+textoidl(' km s^{-1}'),color=cgcolor('red')
xyouts, 43.2, -350, 'Binned fit b=' + string(regfit[1],format='(i0)')+textoidl('\pm')+string(perror[1],format='(i0)')+textoidl(' km s^{-1}'),color=cgcolor('green')

linmix_err, xdata - lum_ref, ydata, post, ysig=ydata_err
oplot, xrange1, median(post.alpha) + median(post.beta)*(xrange1 - lum_ref), color=cgcolor('cyan')
err = 0.5*(quantile_1d(0.84, post.beta) - quantile_1d(0.16,post.beta))
xyouts, 45.3, -400, 'Bayesian slope:' + string(median(post.beta), format='(i0)') + textoidl('\pm') + string(err, format='(i0)'),color=cgcolor('cyan')


endplot
cgfixps, figfile


end
