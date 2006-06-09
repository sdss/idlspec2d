   hist_file = 'sky_res_v5.fits' 
   skyall = mrdfits(hist_file,0)
   skyrej = mrdfits(hist_file,1)
   lrgall = mrdfits(hist_file,2)
   lrgrej = mrdfits(hist_file,3)

   hist_file4 = 'sky_res_v4.fits' 
   skyall4 = mrdfits(hist_file4,0)
   skyrej4 = mrdfits(hist_file4,1)
   lrgall4 = mrdfits(hist_file4,2)
   lrgrej4 = mrdfits(hist_file4,3)
   sig = findgen(241)*0.1 - 12.

   s5 = total(skyrej,2) / total(s5) * 10.
   e5 = sqrt(s5 > 1) 
   a5 = [1.0, 0.02, 0.96, 1.0]
   low = where(sig LT 5.0)
   res5 = mpfitfun('gauss_l2', sig, alog10(s5), 1./sqrt(e5), a5)
   res5 = mpfitfun('gauss_l2', sig[low], alog10(s5[low]), 1./sqrt(e5[low]),res5)

   fit5 = 10^gauss_l2(sig, res5, lvls=lvls5)
   sky_lvl5 = [reverse(res5[1] + lvls5), res5[1], res5[1] - lvls5]
   print, sky_lvl5, format='(11(f8.3))'
   g5 = exp(-(sig-res5[1])^2/2./lvls5[0]^2)/sqrt(2.0 * !Pi * lvls5[0]^2)

   a4 = [1.0, 0.02, 0.96, 1.0]
   s4 = total(skyrej4,2) / total(skyrej4)*10.
   e4 = sqrt(s4 > 1)
   res4 = mpfitfun('gauss_l2', sig, alog10(s4), 1./sqrt(e4), a4)
   res4 = mpfitfun('gauss_l2', sig[low], alog10(s4[low]), 1./sqrt(e4[low]),res4)
   fit4 = 10^gauss_l2(sig, res4, lvls=lvls4)
   g4 = exp(-(sig-res4[1])^2/2./lvls4[0]^2)/sqrt(2.0 * !Pi * lvls4[0]^2)
   sky_lvl4 = [reverse(res4[1] + lvls4), res4[1], res4[1] - lvls4]
   print, sky_lvl4, format='(11(f8.3))'

   dfpsplot, 'plot_sky_reshist.ps', /color, /square
   !p.thick = 6
   djs_plot, sig, s4, xr=[-12, 12], yr=[1.0e-6,0.5], /xs, ps=10,/ylog,/ys, $
     /nodata, charthick=4., chars=1.5, thick=6, $
     xtitle='Normalized Sky Flux Residuals', $
     ytitle='Probability Density', xthick=6, ythick=6
   polyfill, [reverse(sig),sig], [reverse(s5), g5]>10^!y.crange[0] 
   djs_oplot, sig, fit5, color='red'
   oplot, sig, s4, ps=10

   dfpsclose
 

   l5 = total(lrgrej,2) / total(lrgrej)*10.
   e5 = sqrt(l5 > 1)
   a5 = [1.0, 0.02, 0.96, 1.0]
   res5 = mpfitfun('gauss_l2', sig, alog10(l5), 1./sqrt(e5), a5)
   res5 = mpfitfun('gauss_l2', sig[low], alog10(l5[low]), 1./sqrt(e5[low]),res5)
   fit5 = 10^gauss_l2(sig, res5, lvls=lvls5)
   lrg_lvl5 = [reverse(res5[1] + lvls5), res5[1], res5[1] - lvls5]
   print, lrg_lvl5, format='(11(f8.3))'
   g5 = exp(-(sig-res5[1])^2/2./lvls5[0]^2)/sqrt(2.0 * !Pi * lvls5[0]^2)

   l4 = total(lrgrej4,2) / total(lrgrej4)*10.
   a4 = [1.0, 0.02, 0.96, 1.0]
   e4 = sqrt(l4 > 1)
   res4 = mpfitfun('gauss_l2', sig, alog10(l4), 1./sqrt(e4), a4)
   res4 = mpfitfun('gauss_l2', sig[low], alog10(l4[low]), 1./sqrt(e4[low]),res4)
   fit4 = 10^gauss_l2(sig, res4, lvls=lvls4)
   g4 = exp(-(sig-res4[1])^2/2./lvls4[0]^2)/sqrt(2.0 * !Pi * lvls4[0]^2)
   lrg_lvl4 = [reverse(res4[1] + lvls4), res4[1], res4[1] - lvls4]
   print, lrg_lvl4, format='(11(f8.3))'

   dfpsplot, 'plot_lrg_reshist.ps', /color, /square
   !p.thick = 6
   djs_plot, sig, l4, xr=[-12, 12], yr=[1.0e-6,0.5], /xs, ps=10,/ylog,/ys, $
     /nodata, charthick=4., chars=1.5, thick=6, $
     xtitle='Normalized LRG Flux Residuals', $
     ytitle='Probability Density', xthick=6, ythick=6
   polyfill, [reverse(sig),sig], [reverse(l5), g5]>10^!y.crange[0] 
   djs_oplot, sig, fit5, color='red'
   oplot, sig, l4, ps=10
   dfpsclose

   print, sky_lvl5, format='(11(f8.3))'
   print, sky_lvl4, format='(11(f8.3))'
   print, lrg_lvl5, format='(11(f8.3))'
   print, lrg_lvl4, format='(11(f8.3))'
 end
