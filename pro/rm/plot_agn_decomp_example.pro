; plot a 3-panel figure for each of the lowz objects in the AGN-host decomposing sample

pro plot_agn_decomp_example

file='/data3/quasar/yshen/work/agn_host/decomp.fits'
result=mrdfits(file,1)

figfile='/data3/quasar/yshen/work/agn_host/showcase.ps'
begplot,name=figfile,/color

outdir='/data3/quasar/yshen/work/agn_host/'
fitsdir = outdir + 'fits/'
nnn=n_elements(result)

charsize = 1.2 & thick = 4. & xticks = 2L & xminor = 5L
linethick = 0.1 & symsize = 3.
ang = string(197B) & len=0.04
pos1 = [0.1, 0.70, 0.97, 0.965]
pos2 = [0.1, 0.38, 0.97, 0.65]
pos3 = [0.1, 0.06, 0.97, 0.33]
if ~keyword_set(qso_npca) then qso_npca=10
if ~keyword_set(gal_npca) then gal_npca=5

linelist = [3725.94, 3727.24, 3970.072, 4101.73, 4340.46, 4363.21, 4685.71, $
      5158.89, 5199.08, 5302.86, $  ; 4861.3632, 4958.911, 5006.843
      6300.32, 6548.05, 6562.801, $
      6583.45, 6716.44, 6730.82]
nline=n_elements(linelist)
vaclist = linelist
airtovac, vaclist
vaclist = alog10(vaclist)
mwidth = 6.e-4 ; Mask out any pixels within +/- 420 km/s

ytitle = textoidl('Flux Density f_\lambda (10^{-17} erg s^{-1} cm^{-2} ') $
      + ang + textoidl('^{-1})')
for i=0L, nnn-1 do begin
  
  wave=result[i].wave
  flux=result[i].flux
  err=result[i].err
  rm_pca_decomp,wave,flux,err,z=result[i].z,ra=result[i].ra,dec=result[i].dec,/deredden,result=decomp_result, $
      qso_npca=qso_npca,gal_npca=gal_npca
  ind_fit=decomp_result.fit_ind

  restwave=wave/(1.+result[i].z)  
  title='SDSS J'+result[i].sdss_name 
  xrange=[3400,max(restwave)]
  ind_s = where(restwave ge 3400. and wave le 9900.)
  yrange = [-0.5, max(smooth(flux[ind_s],3))*1.2]
  if ~keyword_set(nsmooth) then nsmooth=1L
  plot, restwave, smooth(flux,nsmooth), xrange=xrange, yrange=yrange, /ystyle $
      , charsize=charsize, xthick=thick, ythick=thick, charthick=thick, _extra=extra $
      , /xsty,xticklen=len, yticklen = len/2, pos=pos1, thick=linethick,title=title
  oplot, restwave, result.err, thick=linethick,color=cgcolor('gray')
  gal_recon=result[i].flux_recon - (result[i].flux - result[i].flux_gal)
  qso_recon=result[i].flux_recon - gal_recon
  oplot, restwave[ind_fit], (result[i].flux_recon)[ind_fit], color=cgcolor('red')
  oplot, restwave[ind_fit], (qso_recon)[ind_fit],color=cgcolor('blue')
  oplot, restwave[ind_fit], (gal_recon)[ind_fit],color=cgcolor('dark green')
  xyouts, 0.85,  pos1[3]-0.03, 'z='+string(result[i].z,format='(f5.3)'), /norm $
        , charsize=charsize, charthick=thick, color=cgcolor('red')
  xyouts, 0.7,  pos1[3]-0.03, textoidl('f_H=') + string(result[i].f_H,format='(f0.2)'),/norm,charsize=charsize,$
      charthick=thick,color=cgcolor('red')
  xyouts, 0.15,  pos1[3]-0.03, 'decomposition', /norm, charsize=charsize, charthick=thick
 
  ; plot sigma measurements
  galflux_fit=result[i].flux_gal
  galivar_fit=result[i].ivar
  galivar_fit_use=galivar_fit
  ind = where( (restwave ge 4760. and restwave le 5020.) or $ ; Hbeta
                 (restwave ge 6400. and restwave le 6765.) or $ ; Halpha
                 ;(restwave ge 4320. and restwave le 4400.) or $ ; Hgamma
                 ;(restwave ge 4080. and restwave le 4120.) or $ ; Hdelta
                 ;(restwave ge 4668. and restwave le 4696.) or $ ; HeII
                 (restwave lt 4000.) or (restwave gt 5350) )
  if ind[0] ne -1 then galivar_fit_use[ind]=0.
  vdans = rm_vdispfit(galflux_fit, galivar_fit_use, alog10(wave), zobj=result[i].z, yfit=yfit, $
         eigenfile='spEigenElodie.fits', columns=lindgen(5), npoly=5, dzpix=7, /return_chisq)
  xrange2=[4125., min([5350., max(restwave)])]
  ind_s = where(restwave ge xrange2[0] and restwave le xrange2[1] and wave le 9900. and galivar_fit_use gt 0)
  s_flux=smooth(galflux_fit[ind_s], 3)
  yrange2=[min(s_flux),max(s_flux)*1.2]
  ;ploterror, restwave, galflux_fit, err, xrange=xrange2,charsize=charsize, xthick=thick, ythick=thick, $
  ;    charthick=thick,/xsty,xticklen=len, yticklen = len/2, pos=pos2, thick=linethick, $
  ;    /noerase, yrange=yrange2,/ysty
  plot, restwave, galflux_fit, xrange=xrange2,charsize=charsize, xthick=thick, ythick=thick, $
      charthick=thick,/xsty,xticklen=len, yticklen = len/2, pos=pos2, thick=linethick, $
      /noerase, yrange=yrange2,/ysty
  maskind=where(galivar_fit_use eq 0,nbad)
  galflux_fit_plot=galflux_fit
  galflux_fit_plot[maskind] = 0. ; this sets the flux within the masks to be zero
  ;oplot, restwave[maskind],replicate(yrange2[0],nbad), psym=4,color=cgcolor('cyan'),symsize=0.2
  polyfill, [4760.,5020.,5020.,4760.],[yrange2[0], yrange2[0],yrange2[1],yrange2[1]],/line_fill, $
      orientation=45,noclip=0,color=cgcolor('gray')

  for iline=0,nline-1 do begin
       xx1=linelist[iline]*(1. - mwidth*alog(10.))
       xx2=linelist[iline]*(1. + mwidth*alog(10.))
       yy1=yrange2[0] & yy2=yrange2[1]
       polyfill, [xx1,xx2, xx2, xx1], [yy1, yy1, yy2, yy2], /line_fill, $
          orientation=45, noclip=0,color=cgcolor('cyan')
  endfor
  ind_plot=where(galivar_fit_use gt 0)
  if ind_plot[0] ne -1 then oplot, restwave[ind_plot], yfit[ind_plot], $
     color=cgcolor('red'), psym=symcat(16),symsize=0.2 ; 4125.7221-6796.7293
  xyouts, 0.3, pos2[3]-0.03, textoidl('\sigma_{*}=') $
      + string(round(result[i].sigma_shen),format='(i0)')+textoidl('\pm') $
      + string(round(result[i].sigma_shen_err),format='(i0)') + textoidl(' km s^{-1}'),/norm, $
       charsize=charsize,charthick=thick,color=cgcolor('red')
  xyouts, 0.15, pos2[3]-0.03, 'galaxy', /norm, charsize=charsize, charthick=thick
  xyouts, 0.3, pos2[3]-0.06, textoidl('S/N=')+string(result[i].medsn_sigma,format='(f0.1)'),/norm, $
        charsize=charsize,charthick=thick,color=cgcolor('red')

  ; plot quasar fit
  gal_recon=result[i].flux_recon - (result[i].flux - result[i].flux_gal)
  if result[i].f_H gt 0.1 then qsoflux_fit=result[i].flux-gal_recon else $
      qsoflux_fit=result[i].flux
  objtag=string(result[i].plate,format='(i4.4)')+'-'+string(result[i].mjd,format='(i5.5)') $
         +'-'+string(result[i].fiber,format='(i4.4)')
  fitsfile=fitsdir+objtag+'.fits'
  para=mrdfits(fitsfile, 1, /silent)
  conti_fit=para.conti_fit & line_fit=para.line_fit
  linename=strtrim(para.linename)
  xrange3=[4500, min([5400., max(restwave)])]
  ind_s = where(restwave ge xrange3[0] and restwave le xrange3[1] and wave le 9900.)
  s_flux=smooth(qsoflux_fit[ind_s], 3)
  yrange3=[-max(s_flux)*0.1,max(s_flux)*1.2]
  ; ploterror,restwave, qsoflux_fit, err, xrange=xrange3,charsize=charsize, xthick=thick, ythick=thick, $
  ;     charthick=thick,/xsty,xticklen=len, yticklen = len/2, pos=pos3, thick=linethick, $
  ;    xtitle='Rest Wavelength ['+ang+']',/noerase, yrange=yrange3,/ysty
  plot,restwave, qsoflux_fit, xrange=xrange3,charsize=charsize, xthick=thick, ythick=thick, $
      charthick=thick,/xsty,xticklen=len, yticklen = len/2, pos=pos3, thick=linethick, $
      xtitle='Rest Wavelength ['+ang+']',/noerase, yrange=yrange3,/ysty

  xyouts,0.03,0.5,ytitle,/norm,charsize=charsize,charthick=thick,align=0.5,orien=90
  xyouts, 0.15, pos3[3]-0.03, 'quasar', /norm, charsize=charsize, charthick=thick

  ; now overplot the model
  conti_fit=para.conti_fit & line_fit=para.line_fit
  f_fe_balmer_model = fe_flux_balmer(restwave, conti_fit[3:5])
  f_pl_model = conti_fit[6]*(restwave/3000.0)^conti_fit[7]
  f_poly_model = f_poly_conti(restwave, conti_fit[11:*])
  f_conti_model = f_pl_model + f_fe_balmer_model + f_poly_model
  ind=where(linename eq 'Hbeta_br')
  gauss_broad_hb=line_fit[ind]
  ; determine the empirical FWHM
  gau_arr= get_multi_gaussian_prop(gauss_broad_hb)
  ind=where(linename eq 'Hbeta_na' or strmatch(linename, 'OIII4959*') eq 1 or strmatch(linename, 'OIII5007*') eq 1) 
  gauss_narrow=line_fit[ind]
  ngauss=n_elements(ind)/3L
  f_line_model = manygauss(alog(restwave),gauss_broad_hb)+manygauss(alog(restwave),gauss_narrow)
  f_all_model=f_line_model+f_pl_model+f_poly_model+f_fe_balmer_model

  oplot, restwave, f_pl_model+f_poly_model, color=cgcolor('brown')
  oplot, restwave, f_fe_balmer_model, color=cgcolor('blue')
  oplot, restwave, f_all_model, color=cgcolor('red')
  oplot, restwave, manygauss(alog(restwave), gauss_broad_hb),color=cgcolor('green')
  oplot, exp([gau_arr[3], gau_arr[4]]), 0.5*[gau_arr[5],gau_arr[5]]
  for igauss=0, ngauss-1 do oplot, restwave, onegauss(alog(restwave), gauss_narrow[igauss*3:igauss*3+2]),color=cgcolor('cyan')

  ;message, 'diagnose'

  splog, 'Finished: ', i+1, '/', nnn

endfor



endplot
end
