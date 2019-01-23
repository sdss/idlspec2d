; plot the PS1 light curve
pro ir_plot_ps1

file='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/spitzer/spitzer_c11-12_raw_lc.fits'

result=mrdfits(file,1)
nnn=n_elements(result)

figfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/spitzer/ps1_lc.ps'
begplot,name=figfile,/color,/landscape

for i=0L, nnn-1 do begin

  mjd=result[i].lc_mjd_g & lc=result[i].lc_mag_g & lc_err=result[i].lc_err_g
  ind=where(mjd gt 0)
  mjd=mjd[ind] & lc=lc[ind] & lc_err=lc_err[ind]
  minmjd=min(mjd, max=maxmjd)
  
  ;yrange=[min(lc), max(lc)]
  g_med=median(lc)
  yrange=[g_med+1.5, g_med-1.5]

  pos=[0.08, 0.12, 0.9, 0.95]
  xrange=[55200, 56600]
  xtickname=['55200','55400','55600','55800','56000', $
'56200', '56400', '56600']
  ra=result[i].ra_2 & dec=result[i].dec_2
  radec2string, ra,dec, outstring,rahr=rahr, ramin=ramin, rasec=rasec, $
     decdeg=decdeg, decmin=decmin, decsec=decsec
  title='SDSSJ '+ rahr+ramin+rasec+decdeg+decmin+decsec+', '+'z='+string(result[i].zfinal,format='(f0.4)')
  ploterror, mjd, lc, lc_err,psym=symcat(9),xtitle='MJD',ytitle='PS1 g mag', $
    xrange=xrange, yrange=yrange, title=title, pos=pos, /xsty, xtickname=xtickname
  ; overplot a bunch of information
  xyouts, 0.5, 0.89, /norm, textoidl('logL_{bol}=')+string(result[i].logLbol_2,format='(f0.3)')
  xyouts, 0.7, 0.89, /norm, textoidl('\tau_{IR,obs}=')+string(result[i].TAU_IR_DAY_OBS,format='(i0)')+' days'
  xyouts, 0.15, 0.89, 'W1='+string(result[i].W1MPRO_2,format='(f0.2)')+' ['+string(result[i].W1_MJY, format='(f0.3)')+' mJy]', /norm
  xyouts, 0.5, 0.85, /norm, 'RMS g='+string((result[i].ps1_RMS_mag)[0],format='(f0.2)')

  mjd2datelist, xrange[0], datelist=date0
  mjd2datelist, xrange[1], datelist=date1  
  xyouts, pos[0], pos[3]+0.01, align=0.5, date0,/norm
  xyouts, pos[2], pos[3]+0.01, align=0.5, date1,/norm

endfor
endplot

end

; prepare for IR monitoring programs (ground K, and Spitzer 3.6,4.5um)
pro make_raw_sample_ir

file=getenv('IDLRM_DIR')+'/etc/target_ir_info.fits'
result=mrdfits(file,1)
m_thre = 0.
ps1_rms_g =  result.ps1_rms_mag[0,*] 
ind = where(result.w1mpro gt 0 and result.w1sigmpro gt 0 and ps1_rms_g gt m_thre $
    and result.w1_mjy gt 0.035 and result.tau_ir_day_obs lt 365.25*4., nn)

outfile=getenv('IDLRM_DIR')+'/etc/spitzer_c11-12_raw.fits'
mwrfits, result[ind], outfile, /create

print, nn

end

pro rm_ir_prep, overwrite=overwrite

  ; ugriz: 3551, 4686, 6166, 7480, 8932
  ; K: 2.2um; IRAC-Warm: 3.6, 4.5um

  outfile=getenv('IDLRM_DIR')+'/etc/target_ir_info.fits'
  if file_test(outfile) eq 0 or keyword_set(overwrite) then begin

  ; first compile a fits file including IR properties (mag, flux, lag)
  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  target=target[0:848]  ; keep only RM targets
  nobj=n_elements(target)

  ; in units of A/s
  c = 2.99792458d18
  d10pc = 3.08d19

  ; read in SED file from Richards++(2006)
  ; units are log[Hz] and log[erg/s]
  SED_file = '/home/yshen/Research/IDL/lib/Projects/quasar/quasar_bolometric_corr.txt'
  readcol, sed_file, format = 'f,f', log_nu, log_nuLnu, skipline = 29L, /silent

  result={K_ab:0.D, K_mJy:0.D, IRAC36_ab:0.D, IRAC36_mJy:0.D, $
          IRAC45_ab:0.D, IRAC45_mJy:0.D,Mi_z2:0.D,logLbol:0.D, $
          tau_IR_pc:0.D, tau_IR_day_obs:0.D, $
          w1mpro:0.D, w1sigmpro:-1.D, w2mpro:0.D, w2sigmpro:-1.d, $
          w3mpro:0.D, w3sigmpro:-1.D, w4mpro:0.D, w4sigmpro:-1.D, $
          w1_ab:0.D, w1_mJy:0.D, w2_ab:0.D, w2_mJy:0.D, w3_ab:0.D, w3_mJy:0.D, $
          w4_ab:0.D, w4_mJy:0.D, $
          ps1_rms_mag:replicate(-1.D,5)}
  result=replicate(result, nobj)
  result=struct_addtags(target, result)
  ugriz=result.psfmag

  ; now compute K and IRAC-Warm using the average SED
  for i=0L, nobj-1 do begin
  
    ; i, K, IRAC3.6, 4.5
    nu_arr=c/( [7480., 2.2d4, 3.6d4, 4.5d4]/(1. + result[i].zfinal) )
    log_nu_arr = alog10(nu_arr)
    log_Lnu_arr = interpol(log_nuLnu, log_nu, log_nu_arr) - log_nu_arr

    ; AB mag, assuming the sinh psfmag is approximately ab
    lups=ugriz[*,i]
    maggies=k_lups2maggies(lups)
    k_abfix, maggies, replicate(0.02,5)
    maggie2mag=22.5-2.5*alog10(maggies*1d9)
    imag=maggie2mag[3] ; this is now ab

    m_ab = -2.5*(log_Lnu_arr - log_Lnu_arr[0]) + imag

    ; mJy
    mJy = 10.^( m_ab*(-0.4) )*3631d3

    result[i].K_ab=m_ab[1] & result[i].K_mJy=mJy[1]
    result[i].IRAC36_ab=m_ab[2] & result[i].IRAC36_mJy=mJy[2]
    result[i].IRAC45_ab=m_ab[3] & result[i].IRAC45_mJy=mJy[3]

    Mi_z2=get_abs_mag(imag, result[i].zfinal)
    logL2500 = -0.4*(mi_z2 + 48.6 + 2.5*alog10(3.)) + alog10(4d*!PI*d10pc^2) + alog10(3d18/2500d)
    logL5100 = logL2500 - 0.5*alog10(5100./2500.)
    Lbol = 10.D^logL5100 * 9.26
    result[i].logLbol = alog10(Lbol)
    result[i].Mi_z2 = Mi_z2

    ; inner torus dust radius from Nenkova et al. (2008, paperII, eqn 1)
    Tsub=1500. ; K
    Rd_pc = 0.4*(Lbol/1d45)^0.5*(1500./Tsub)^2.6
    result[i].tau_IR_pc = Rd_pc
    result[i].tau_IR_day_obs = Rd_pc*3.08568d18/2.99792458d10/(3600.D*24.)*(1. + result[i].zfinal)

  endfor

  ; now add PS1 variability properties
  file = '/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/PS1/targets_PS1MD07b.fits'  
  ps1=mrdfits(file,1)
  spherematch, result.ra, result.dec, ps1.ra_dr10, ps1.dec_dr10, 0.2/3600.D, match0, match1, dist1
  ; NB, ps1 filters are grizy
  result[match0].ps1_rms_mag=ps1[match1].stdev
  result[match0].w1mpro=ps1[match1].w1mpro & result[match0].w1sigmpro=ps1[match1].w1sigmpro
  result[match0].w2mpro=ps1[match1].w2mpro & result[match0].w2sigmpro=ps1[match1].w2sigmpro
  result[match0].w3mpro=ps1[match1].w3mpro & result[match0].w3sigmpro=ps1[match1].w3sigmpro
  result[match0].w4mpro=ps1[match1].w4mpro & result[match0].w4sigmpro=ps1[match1].w4sigmpro
   
  result[match0].w1_ab = result[match0].w1mpro + 2.683
  result[match0].w1_mJy = 10.^( result[match0].w1_ab*(-0.4) )*3631d3
  result[match0].w2_ab = result[match0].w2mpro + 3.319
  result[match0].w2_mJy = 10.^( result[match0].w2_ab*(-0.4) )*3631d3
  result[match0].w3_ab = result[match0].w3mpro + 5.242
  result[match0].w3_mJy = 10.^( result[match0].w3_ab*(-0.4) )*3631d3
  result[match0].w4_ab = result[match0].w4mpro + 6.604
  result[match0].w4_mJy = 10.^( result[match0].w4_ab*(-0.4) )*3631d3

  ; output the file
  mwrfits, result, outfile, /create

  endif else result = mrdfits(outfile, 1)

  ; now make diagnostic plots

  figfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/spitzer/flux_tau.ps'
  begplot, name=figfile,/color, encap=0
  !p.multi = [0,1,2]
  ; first plot for 3.6um 
  ind = where(result.w1mpro gt 0 and result.w1sigmpro gt 0, n0)
  ps1_rms_g = result.ps1_rms_mag[0,*]
  m_thre=0.25 ; > this m_thre will be flagged
  ind1 = where(result.w1mpro gt 0 and result.w1sigmpro gt 0 and ps1_rms_g gt m_thre, n1)
  plot, result[ind].zfinal, result[ind].w1_Mjy, yrange=[0.01, 3],/ysty, /ylog, $
    xtitle='Redshift', ytitle = textoidl('W1(3.4 \mum) Flux Density [mJy]'), psym=symcat(7)
  oplot, result[ind1].zfinal, result[ind1].w1_mjy, psym=symcat(9,thick=5),color=cgcolor('red')
  xyouts, 0.6, 0.91,'W1 detection: '+ string(n0, format='(i0)'),/norm
  xyouts, 0.6, 0.88, textoidl('RMS g_{PS1}>')+string(m_thre,format='(f0.2)')+': '+ string(n1, format='(i0)'),/norm,color=cgcolor('red')

  ;plot, result[ind].zfinal, result[ind].TAU_IR_DAY_OBS, yrange=[100., 5d4],/ysty, /ylog, $
  ;  xtitle='Redshift', ytitle = textoidl('\tau_{hot dust} [observed days]'), psym=symcat(7)
  ;oplot, result[ind1].zfinal, result[ind1].TAU_IR_DAY_OBS, psym=symcat(9,thick=5),color=cgcolor('red')
  ;oplot, [0,5], [365.25,365.25]*2, line=2,thick=5
  ;xyouts, 3,480, '2yr Spitzer Cy11-12'
  ;oplot, [0,5], [365.25,365.25]*4, line=2, thick=5,color=cgcolor('cyan')
  ;xyouts, 3, 1000, '4yr early opt LC', color=cgcolor('cyan')
  plot, result[ind].w1_mjy, result[ind].TAU_IR_DAY_OBS, xrange=[0.01,3],/xlog, yrange=[100., 5d4], $
    /ysty, /ylog, /xsty,  $
    xtitle=textoidl('W1(3.4 \mum) Flux Density [mJy]'), ytitle = textoidl('\tau_{hot dust} [observed days]'), psym=symcat(7)
  oplot, result[ind1].w1_mjy, result[ind1].TAU_IR_DAY_OBS, psym=symcat(9,thick=5),color=cgcolor('red')
  oplot, [0.01,3], [365.25,365.25]*2, line=2,thick=5
  ;xyouts, 0.1,480, '2yr Spitzer Cy11-12'
  oplot, [0.01,3], [365.25,365.25]*4, line=2, thick=5,color=cgcolor('cyan')
  ;xyouts, 0.1, 1000, '4yr early opt LC', color=cgcolor('cyan')
  item=['2yr Spitzer Cy11-12', '4yr early opt LC']
  colors=cgcolor(['black','cyan'])
  legend, item,box=0,pos=[0.5, 0.15],/norm, color=cgcolor(['black','cyan']),line=[2,2],textcolor=colors

   ; then plot for 4.5um 
  ind = where(result.w2mpro gt 0 and result.w2sigmpro gt 0, n0)
  ps1_rms_g = result.ps1_rms_mag[0,*]
  m_thre=0.25 ; > this m_thre will be flagged
  ind1 = where(result.w2mpro gt 0 and result.w2sigmpro gt 0 and ps1_rms_g gt m_thre, n1)
  plot, result[ind].zfinal, result[ind].w2_Mjy, yrange=[0.01, 3],/ysty, /ylog, $
    xtitle='Redshift', ytitle = textoidl('W2(4.6 \mum) Flux Density [mJy]'), psym=symcat(7)
  oplot, result[ind1].zfinal, result[ind1].w2_mjy, psym=symcat(9,thick=5),color=cgcolor('red')
  xyouts, 0.6, 0.91,'W2 detection: '+ string(n0, format='(i0)'),/norm
  xyouts, 0.6, 0.88, textoidl('RMS g_{PS1}>')+string(m_thre,format='(f0.2)')+': '+ string(n1, format='(i0)'),/norm,color=cgcolor('red')

  plot, result[ind].zfinal, result[ind].TAU_IR_DAY_OBS, yrange=[100., 5d4],/ysty, /ylog, $
    xtitle='Redshift', ytitle = textoidl('\tau_{hot dust} [observed days]'), psym=symcat(7)
  oplot, result[ind1].zfinal, result[ind1].TAU_IR_DAY_OBS, psym=symcat(9,thick=5),color=cgcolor('red')
  oplot, [0,5], [365.25,365.25]*2, line=2,thick=5
  xyouts, 3.,480, '2yr Spitzer Cy11-12'
  oplot, [0,5], [365.25,365.25]*4, line=2, thick=5,color=cgcolor('cyan')
  xyouts, 3., 1000, '4yr early opt LC', color=cgcolor('cyan')


  !p.multi = 0
  endplot



end
