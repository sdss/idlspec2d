;+
; NAME:
;
; PURPOSE:
;
; INPUTS:
;
; OPTIONAL INPUTS;
;
; OUTPUTS:
;
;
;--------------------------

pro rm_agn_host_decomp,plate,fiber,mjd,tag=tag,figfile=figfile,calibdir=calibdir,$
      ra=ra,dec=dec,z=z,nsmooth=nsmooth, qso_npca=qso_npca, gal_npca=gal_npca, $
      pop_qsofit_err=pop_qsofit_err,outdir=outdir, nofits=nofits, update_sigma=update_sigma, $
      refit=refit

  if n_elements(plate) eq 0 then begin ; default to fit all z<1.1 RM quasars
    splog, 'Load default target list: '
    file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
    target=mrdfits(file,1)
    ; keep only quasars (i.e., std and sky will have zfinal=0)
    ind = where(target.zfinal le 1.1 and target.zfinal gt 0)
    ra=target[ind].ra & dec=target[ind].dec & z=target[ind].zfinal
    plate=0 & fiber=ind+1 & mjd=56837L
  endif

  if ~keyword_set(outdir) then outdir='/data3/quasar/yshen/work/agn_host/'
  if file_test(outdir) eq 0 then spawn, 'mkdir ' + outdir
  fitsdir = outdir + 'fits/'
  if file_test(fitsdir) eq 0 then spawn, 'mkdir ' + fitsdir
  errdir = outdir + 'err/'
  if file_test(errdir) eq 0 then spawn, 'mkdir ' + errdir

  if ~keyword_set(tag) then tag='decomp'
  outfile=outdir + tag + '.fits'
  if file_test(outfile) eq 1 and keyword_set(update_sigma) then begin
     output0=mrdfits(outfile,1)
     nofits = 1L
  endif


  ; setup cosmology
  red, omegalambda=0.7,omega0=0.3,h100=0.7

  nplate=n_elements(plate) & nfiber=n_elements(fiber) & nmjd=n_elements(mjd)
  if nplate lt nfiber then plate=replicate(plate[0],nfiber)
  if nmjd lt nfiber then mjd=replicate(mjd[0],nfiber)

  if ~keyword_set(calibdir) then calibdir='wh_skysub/'  
  if ~keyword_set(nsmooth) then nsmooth=1L

  ; setup the output Figure
  if ~keyword_set(figfile) then figfile=outdir + tag+'.ps'
  begplot,name=figfile,/landscape,/color
  charsize = 1. & thick = 4. & xticks = 2L & xminor = 5L
  linethick = 0.1 & symsize = 3.
  ang = string(197B) & len=0.04     
  pos = [0.08, 0.5, 0.97, 0.96]  ;else pos = [0.12, 0.14, 0.96, 0.98]
  pos2 = [0.08, 0.06, 0.97, 0.42]
  ytitle = textoidl('Flux Density f_\lambda (10^{-17} erg s^{-1} cm^{-2} ') $
      + ang + textoidl('^{-1})')
  xtitle = textoidl('Rest Wavelength (')+ang+')'

  nobj=n_elements(plate)
  rm_readspec,plate[0],fiber[0],mjd=mjd[0],calibdir=calibdir,wave=wave
  npix=n_elements(wave)
  output={sdss_name:'', ra:0.D, dec:0.D, plate:0L, fiber:0L, mjd:0L, z:0.D, $
       wave:wave, flux:dblarr(npix), $
       err:dblarr(npix), ivar:dblarr(npix), flux_recon:dblarr(npix), $
       flux_gal:dblarr(npix), f_H:0., f_H_5100:0.,decomp_redchi2:0.,npca_qso:0L, npca_gal:0L, $
       medsn_sigma:0., npix_sigma:0L, $
       sigma_shen:0.D, sigma_shen_err:-1.D, sigma_redchi2:0.D, $
       sigma2:0.D, sigma2_err:-1.D, sigma2_redchi2:0.D, $
       sigma_greene:0.D, sigma_greene_err:-1.D, $
       logL5100_qso:0.D, logL5100_qso_err:-1.D, qso_c_redchi2:0.D, $
       FWHM_hb:0.D, FWHM_hb_err:-1.D, FWHM_OIII5007:0.D, FWHM_OIII5007_err:-1.D, $
       FWHM_OIII5007c:0.D, FWHM_OIII5007c_err:-1.D, $
       FWHM_OIII5007w:0.D, FWHM_OIII5007w_err:-1.D, $
       qso_l_redchi2:0.D, logMbh_hb:0.D, logMbh_hb_err:-1.D }
  output=replicate(output, nobj)
  output.ra=ra & output.dec=dec & output.plate=plate & output.fiber=fiber
  output.mjd=mjd & output.z=z
 
  for i=0L, nobj - 1 do begin

    rm_readspec,plate[i],fiber[i],mjd=mjd[i],calibdir=calibdir, $
      wave=wave,flux=flux,flerr=err,invvar=ivar
    restwave=wave/(1.+z[i])

    ;-------------------------- 
    ; decomposing the spectrum
    if ~keyword_set(qso_npca) then qso_npca=10 
    if ~keyword_set(gal_npca) then gal_npca=5
    rm_pca_decomp,wave,flux,err,z=z[i],ra=ra[i],dec=dec[i],/deredden,result=result, $
      qso_npca=qso_npca,gal_npca=gal_npca
    ; get residual in the PCA fitting range
    residual = flux[result.fit_ind] - (result.flux_recon)[result.fit_ind]
    norm_residual = residual * sqrt(ivar[result.fit_ind])
    ind_fit=result.fit_ind

    xrange=[3400,max(restwave)]
    ind_s = where(restwave ge 3400. and wave le 9900.)
    yrange = [-0.5, max(smooth(flux[ind_s],20))*1.2]
    plot, restwave, smooth(result.flux,nsmooth), xrange=xrange, yrange=yrange, /ystyle $
      , charsize=charsize, xthick=thick, ythick=thick, charthick=thick, _extra=extra $
      , /xsty,xticklen=len, yticklen = len/2, pos=pos, thick=linethick,title='AGN/host decomposition'
    oplot, restwave, result.err, thick=linethick,color=cgcolor('gray')
    oplot, restwave[ind_fit], (result.flux_recon)[ind_fit], color=cgcolor('red')
    oplot, restwave[ind_fit], (result.qso_recon)[ind_fit],color=cgcolor('blue')
    oplot, restwave[ind_fit], (result.gal_recon)[ind_fit],color=cgcolor('dark green')

    objtag=string(plate[i],format='(i4.4)')+'-'+string(mjd[i],format='(i5.5)') $
         +'-'+string(fiber[i],format='(i4.4)')
    radec2string,ra[i],dec[i],dummy,rahr=rahr,ramin=ramin,rasec=rasec $
        , decdeg=decdeg,decmin=decmin,decsec=decsec
    sdss_name=rahr+ramin+rasec+decdeg+decmin+decsec
    xyouts, 0.15, 0.922, objtag, /norm, charsize = charsize, charthick=thick
    xyouts, 0.15, 0.892, '68% norm residual=' + $
      string( quantile_1d(0.68, abs(norm_residual)), format='(f4.2)'), /norm,charsize=charsize,charthick=thick

    xyouts, 0.35, 0.922, /norm, SDSS_name, charsize = charsize, charthick = thick
    xyouts, 0.85, 0.922, 'z='+string(z[i],format='(f5.3)'), /norm $
        , charsize=charsize, charthick=thick, color=cgcolor('red')
    xyouts,0.5,0.01,xtitle,/norm,charsize=charsize,charthick=thick,align=0.5
    xyouts,0.02,0.5,ytitle,/norm,charsize=charsize,charthick=thick,align=0.5,orien=90
    f_H=result.f_H
    xyouts, 0.65, 0.922, textoidl('f_H=') + string(f_H,format='(f0.2)'),/norm,charsize=charsize,$
      charthick=thick
    
    output[i].sdss_name=sdss_name
    output[i].wave=wave & output[i].flux=result.flux & output[i].err=result.err
    output[i].ivar=result.ivar & output[i].flux_recon=result.flux_recon
    output[i].flux_gal=result.flux-result.qso_recon & output[i].npca_qso=result.npca_qso
    output[i].npca_gal=result.npca_gal & output[i].decomp_redchi2=result.chi2/result.dof
    output[i].f_H=result.f_H & output[i].f_H_5100=result.f_H_5100
    ; set the pixels outside ind_fit to be zero
    ind_min=min(ind_fit) & ind_max=max(ind_fit)
    if ind_min gt 0 then output[i].flux_gal[0:ind_min-1]=0 
    if ind_max lt npix-1 then output[i].flux_gal[ind_max+1:*]=0

    ;-------------------------
    ; fit host stellar velocity dispersion
    ; the wavelength range in the vdisp fit is always smaller than the PCA galaxy/qso template range
    galflux_fit=result.flux-result.qso_recon  
    galivar_fit=result.ivar
    galivar_fit_use=galivar_fit
    wave_max_use=4600. ; only fit dispersion to the short wavelength
    ; mask the regions around broad Hbeta and broad Halpha
    ind = where( (restwave ge 4760. and restwave le 5020.) or $ ; Hbeta
                 (restwave ge 6400. and restwave le 6765.) or $ ; Halpha
                 ;(restwave ge 4320. and restwave le 4400.) or $ ; Hgamma
                 ;(restwave ge 4080. and restwave le 4120.) or $ ; Hdelta
                 ;(restwave ge 4668. and restwave le 4696.) or $ ; HeII
                 (restwave lt 4000.) or (restwave gt 5350) )
    if ind[0] ne -1 then galivar_fit_use[ind]=0.
    ind_good=where(restwave ge 4125. and restwave le 5350 and galivar_fit_use gt 0, npix_sigma)
    output[i].npix_sigma=npix_sigma
    output[i].medsn_sigma=median(galflux_fit[ind_good]*sqrt(galivar_fit_use[ind_good]) )

    vdans = rm_vdispfit(galflux_fit, galivar_fit_use, alog10(wave), zobj=z[i], yfit=yfit, $
         eigenfile='spEigenElodie.fits', columns=lindgen(5), npoly=5, dzpix=7, /return_chisq)
    output[i].sigma_SHEN = vdans.vdisp & output[i].sigma_SHEN_err = vdans.vdisp_err
    output[i].sigma_redchi2 = vdans.vdispchi2/vdans.vdispdof
    if file_test(outfile) eq 1 and keyword_set(update_sigma) then begin
        output0[i].sigma_SHEN = vdans.vdisp & output0[i].sigma_SHEN_err = vdans.vdisp_err
        output0[i].sigma_redchi2 = vdans.vdispchi2/vdans.vdispdof
    endif

    ; this sigma fit includes the Ca H/K bands
    ind2 = where( (restwave ge 4760. and restwave le 5020.) or $
      (restwave ge 6400. and restwave le 6765.) ) ; or restwave gt wave_max_use )
    galivar_fit_use2=galivar_fit
    if ind2[0] ne -1 then galivar_fit_use2[ind2]=0.
    vdans2 = rm_vdispfit(galflux_fit, galivar_fit_use2, alog10(wave), zobj=z[i], yfit=yfit2, $
         eigenfile='spEigenElodie.fits', columns=lindgen(5), npoly=5, dzpix=7, /return_chisq )
    output[i].sigma2=vdans2.vdisp & output[i].sigma2_err=vdans2.vdisp_err
    output[i].sigma2_redchi2=vdans2.vdispchi2/vdans2.vdispdof
    if file_test(outfile) eq 1 and keyword_set(update_sigma) then begin
        output0[i].sigma2=vdans2.vdisp & output0[i].sigma2_err=vdans2.vdisp_err
        output0[i].sigma2_redchi2=vdans2.vdispchi2/vdans2.vdispdof
    endif

    xrange2=[3900., min([6100., max(restwave)]) ]
    s_flux=smooth(galflux_fit[ind_s], 10)
    yrange2=[min(s_flux),max(s_flux)]
    plot, restwave, galflux_fit, xrange=xrange2,  charsize=charsize, xthick=thick, ythick=thick, $
      charthick=thick,/xsty,xticklen=len, yticklen = len/2, pos=pos2, thick=linethick, $
       title='stellar velocity dispersion',/noerase, yrange=yrange2,/ysty
    maskind=where(galivar_fit_use eq 0,nbad)
    galflux_fit_plot=galflux_fit
    galflux_fit_plot[maskind] = 0. ; this sets the flux within the masks to be zero
    ; oplot, restwave, galflux_fit_plot, thick=linethick
    oplot, restwave[maskind],replicate(yrange2[0],nbad), psym=4,color=cgcolor('cyan'),symsize=0.2
    ; now plot the narrow masks from rm_vdispfit 
    linelist = [3725.94, 3727.24, 3970.072, 4101.73, 4340.46, 4363.21, 4685.71, $
       4861.3632, 4958.911, 5006.843, 5158.89, 5199.08, 5302.86, $
       6300.32, 6548.05, 6562.801, $
       6583.45, 6716.44, 6730.82]
    nline=n_elements(linelist)
    vaclist = linelist
    airtovac, vaclist
    vaclist = alog10(vaclist)
    mwidth = 6.e-4 ; Mask out any pixels within +/- 420 km/s
    for iline=0,nline-1 do begin
       xx1=linelist[iline]*(1. - mwidth*alog(10.))
       xx2=linelist[iline]*(1. + mwidth*alog(10.))
       yy1=yrange2[0] & yy2=yrange2[1]
       polyfill, [xx1,xx2, xx2, xx1], [yy1, yy1, yy2, yy2], /line_fill, $
          orientation=45, /noclip,color=cgcolor('cyan')
    endfor

    ind2=where(restwave le 7000.)
    if ind2[0] ne -1 then oplot, restwave[ind2], yfit2[ind2], color=cgcolor('blue') ; all

    ind_plot=where(galivar_fit_use gt 0)
    if ind_plot[0] ne -1 then $
      oplot, restwave[ind_plot], yfit[ind_plot],color=cgcolor('red'), psym=symcat(16),symsize=0.1 ; 4000-5350, same as what Greene uses; but the Elodie eigenvectors run from 4125.7221 -- 6796.7293

    xyouts, 0.15, pos2[3]-0.04, textoidl('\sigma_{*}=') $
      + string(vdans.vdisp,format='(i0)')+textoidl('\pm') $
      + string(vdans.vdisp_err,format='(i0)') + textoidl(' km s^{-1}'),/norm,charsize=charsize,charthick=thick,$
      color=cgcolor('red')
    xyouts, 0.65, pos2[3]-0.04, textoidl('redchi2=')+string(vdans.vdispchi2/vdans.vdispdof, $
      format='(f0.1)'),/norm, charsize=charsize,charthick=thick, color=cgcolor('red')
    xyouts, 0.15, pos2[3]-0.07, textoidl('\sigma_{*}=') $
      + string(vdans2.vdisp,format='(i0)')+textoidl('\pm') $
      + string(vdans2.vdisp_err,format='(i0)') + textoidl(' km s^{-1}'),/norm, $
      charsize=charsize,charthick=thick,color=cgcolor('blue')
    xyouts, 0.65, pos2[3]-0.07, textoidl('redchi2=')+string(vdans2.vdispchi2/vdans2.vdispdof, $
      format='(f0.1)'),/norm, charsize=charsize,charthick=thick, color=cgcolor('blue')

    ;-------------------------
    ; fit the AGN spectrum and output fits/err fits 
    ; I am fitting everything, regardless of the f_H value
    ; if f_H < 0.1 then all raw flux is AGN flux
    ; note that the spectrum has already been dereddened

    if result.f_H gt 0.1 then qsoflux_fit=result.flux-result.gal_recon else $
      qsoflux_fit=result.flux 

    qsoivar_fit=result.ivar
    qsoivar_use=qsoivar_fit
    ; mask the regions in the galaxy NLs
    
    emparfile='/home/yshen/products/Linux/idlrm/etc/qsoline_qsovar.par'
    fitsfile=fitsdir+objtag+'.fits'
    if file_test(fitsfile) eq 0 or keyword_set(refit) then begin ; fit for this object if not already fit
      rm_qsofit,wave[ind_fit], qsoflux_fit[ind_fit], dummy, z[i], ivar0=qsoivar_use[ind_fit], $
       ra=ra[i],dec=dec[i], emparfile=emparfile, deredden=0, /append,/psplot, $
       objtag=objtag,sdss_name=sdss_name, para=para
      ; output the fitsfile
      mwrfits,para,fitsfile, /create
    endif else para=mrdfits(fitsfile,1,/silent)
    output[i].qso_c_redchi2=para.conti_redchi2 
    ind_hb=where(para.linename eq 'Hbeta_br')
    output[i].qso_l_redchi2=(para.line_redchi2)[ind_hb[0]]
    get_qso_prop, para, linelist=['Hbeta_br', 'OIII5007', 'OIII5007c', 'OIII5007w'], $
        conti_prop=conti_prop, line_prop=line_prop,z=z[i]
    output[i].logL5100_QSO=conti_prop.logL5100
    output[i].FWHM_Hb=(line_prop.Hbeta_br)[1]
    output[i].FWHM_OIII5007=(line_prop.OIII5007)[1]
    output[i].FWHM_OIII5007c=(line_prop.OIII5007c)[1]
    output[i].FWHM_OIII5007w=(line_prop.OIII5007w)[1]

    ; estimate errors with MC trials
    errfile=errdir+objtag+'.fits'
    if file_test(errfile) eq 0 or keyword_set(refit) then begin
      if keyword_set(pop_qsofit_err) then begin
        ntrial=50L ; number of MC trials
        if n_elements(err_array) ne 0 then tmp=temporary(err_array) ; empty err_array
        if n_elements(logL5100_arr) ne 0 then tmp1=temporary(logL5100_arr)
        if n_elements(FWHM_Hb_arr) ne 0 then tmp2=temporary(FWHM_Hb_arr) 
        for jj=0L, ntrial-1 do begin
          rm_qsofit,wave[ind_fit], qsoflux_fit[ind_fit], dummy, z[i], ivar0=qsoivar_use[ind_fit], $
          /noplot,/silent,/diet,/add_noise,emparfile=emparfile,deredden=0,para=para
          tag_rej=['CONTI_FIT_ERR','LINE_FIT_ERR']
          remove_tags, para, tag_rej, para_keep
          splog, '  Error trial: ',jj+1
          ; Only add this trial if the continuum fit is OK
          ; The large redchi2 is to keep fits on very high SN spectrum
          if para_keep.conti_redchi2 lt 100. and para_keep.conti_redchi2 gt 0 then begin
             if n_elements(err_array) eq 0 then err_array=para_keep $
             else err_array=[err_array,para_keep]
             get_qso_prop, para_keep, linelist=['Hbeta_br'],conti_prop=conti_prop,line_prop=line_prop,z=z[i]
             if n_elements(logL5100_arr) eq 0 then logL5100_arr=conti_prop.logL5100 $
             else logL5100_arr=[logL5100_arr, conti_prop.logL5100]
             if n_elements(FWHM_Hb_arr) eq 0 then FWHM_hb_arr=(line_prop.Hbeta_br)[1] $
             else FWHM_Hb_arr=[FWHM_Hb_arr, (line_prop.Hbeta_br)[1]]
          endif
        endfor
        output[i].logL5100_QSO_err=0.5*( quantile_1d(0.84,logL5100_arr) - quantile_1d(0.16,logL5100_arr) )
        output[i].FWHM_Hb_err=0.5*( quantile_1d(0.84,FWHM_Hb_arr) - quantile_1d(0.16,FWHM_Hb_arr) )
        splog, 'Finished: '+objtag,' ', i+1, '/', nobj
        print, output[i].logL5100_QSO, output[i].logL5100_QSO_err
        print, output[i].FWHM_Hb, output[i].FWHM_Hb_err
        ; message, 'stop and diagnose'
      
        ; output the error file
        if n_elements(err_array) ne 0 then mwrfits, err_array, errfile, /create
      endif
    endif else begin 
        err_array=mrdfits(errfile, 1, /silent)
        ntrial_good=n_elements(err_array)
        if n_elements(logL5100_arr) ne 0 then tmp1=temporary(logL5100_arr)
        if n_elements(FWHM_Hb_arr) ne 0 then tmp2=temporary(FWHM_Hb_arr)
        if n_elements(FWHM_OIII5007_arr) ne 0 then tmp3=temporary(FWHM_OIII5007_arr)
        if n_elements(FWHM_OIII5007c_arr) eq 0 then tmp4=temporary(FWHM_OIII5007c_arr)
        if n_elements(FWHM_OIII5007w_arr) eq 0 then tmp5=temporary(FWHM_OIII5007w_arr)
        for itrial=0L, ntrial_good - 1 do begin
           get_qso_prop, err_array[itrial], linelist=['Hbeta_br', 'OIII5007', 'OIII5007c', 'OIII5007w'], $
             conti_prop=conti_prop,line_prop=line_prop,z=z[i]
           if n_elements(logL5100_arr) eq 0 then logL5100_arr=conti_prop.logL5100 $
            else logL5100_arr=[logL5100_arr, conti_prop.logL5100]
           if n_elements(FWHM_Hb_arr) eq 0 then FWHM_hb_arr=(line_prop.Hbeta_br)[1] $
            else FWHM_Hb_arr=[FWHM_Hb_arr, (line_prop.Hbeta_br)[1]]
           if n_elements(FWHM_OIII5007_arr) eq 0 then FWHM_OIII5007_arr=(line_prop.OIII5007)[1] $
            else FWHM_OIII5007_arr=[FWHM_OIII5007_arr, (line_prop.OIII5007)[1]]
           if n_elements(FWHM_OIII5007c_arr) eq 0 then FWHM_OIII5007c_arr=(line_prop.OIII5007c)[1] $
            else FWHM_OIII5007c_arr=[FWHM_OIII5007c_arr, (line_prop.OIII5007c)[1]]
           if n_elements(FWHM_OIII5007w_arr) eq 0 then FWHM_OIII5007w_arr=(line_prop.OIII5007w)[1] $
            else FWHM_OIII5007w_arr=[FWHM_OIII5007w_arr, (line_prop.OIII5007w)[1]]
        endfor
        ; now get errors
        output[i].logL5100_QSO_err=0.5*( quantile_1d(0.84,logL5100_arr) - quantile_1d(0.16,logL5100_arr) )
        output[i].FWHM_Hb_err=0.5*( quantile_1d(0.84,FWHM_Hb_arr) - quantile_1d(0.16,FWHM_Hb_arr) )
        output[i].FWHM_OIII5007_err=0.5*(quantile_1d(0.84,FWHM_OIII5007_arr) - quantile_1d(0.16,FWHM_OIII5007_arr))
        output[i].FWHM_OIII5007c_err=0.5*(quantile_1d(0.84,FWHM_OIII5007c_arr) - quantile_1d(0.16,FWHM_OIII5007c_arr))
        output[i].FWHM_OIII5007w_err=0.5*(quantile_1d(0.84,FWHM_OIII5007w_arr) - quantile_1d(0.16,FWHM_OIII5007w_arr))
        splog, 'Finished: '+objtag,' ', i+1, '/', nobj

    endelse

  endfor

  endplot

  ; now populate Jenny Greene's sigma measurements
  file = '/data3/quasar/yshen/work/agn_host/greene/final/valdes_4000-5350.dat'
  readcol,file,format='a,d,d', sdss_name_g, sigma_greene, sigma_greene_err
  nobj_g=n_elements(sdss_name_g)
  for i_g=0,nobj_g-1 do begin

    ind_tmp=where(output.sdss_name eq sdss_name_g[i_g])
    output[ind_tmp].sigma_greene=sigma_greene[i_g]
    output[ind_tmp].sigma_greene_err=sigma_greene_err[i_g]

  endfor


  if ~keyword_set(nofits) then $
    mwrfits, output, outfile, /create

  ; only update the sigma measurements
  if file_test(outfile) eq 1 and keyword_set(update_sigma) then $
      mwrfits, output0, outfile, /create


end
