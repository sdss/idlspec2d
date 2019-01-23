; get the regression results for the MCMC runs
; on the M-sigma regression

pro get_reg_results

outdir='/data3/quasar/yshen/work/agn_host/mcmc_reg/'
tag=['post_all.fits', 'post_z_gt_0d6.fits', 'post_z0.12-0.42.fits', 'post_z0.42-0.60.fits', $
     'post_z0.61-0.76.fits', 'post_z0.76-1.00.fits']

for i=0L, 5L do begin
  result=mrdfits(outdir+tag[i],1,/silent)

  alpha=(result.alpha)
  beta=(result.beta)
  scat=sqrt(result.sigsqr)

  a0=median(alpha) & a1=quantile_1d(0.16, alpha) & a2=quantile_1d(0.84, alpha)
  str1='$' + string(a0,format='(f0.3)')+'_{-'+string(a0-a1,format='(f0.3)') + $
    '}^{+' + string(a2-a0,format='(f0.3)')+'}$ & '

  a0=median(beta) & a1=quantile_1d(0.16, beta) & a2=quantile_1d(0.84, beta)
  str2='$' + string(a0,format='(f0.3)')+'_{-'+string(a0-a1,format='(f0.3)') + $
    '}^{+' + string(a2-a0,format='(f0.3)')+'}$ & '

  a0=median(scat) & a1=quantile_1d(0.16, scat) & a2=quantile_1d(0.84, scat)
  str3='$' + string(a0,format='(f0.3)')+'_{-'+string(a0-a1,format='(f0.3)') + $
    '}^{+' + string(a2-a0,format='(f0.3)')+'}$'

  print, str1+str2+str3

endfor

end
