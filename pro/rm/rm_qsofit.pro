;+
; NAME:
;   rm_qsofit
;
; PURPOSE:
;   Globally fit the SDSS QSO spectrum from 1350A to 7180A
;
; CALLING SEQUENCE:
;   rm_qsofit, obs_wave, flux, err, z, [ra=,dec=,/psplot,/fits]
;
;
; INPUTS:
;   lam0       -  Wavelength array of the spectrum [NPIX]
;   flux0      -  Flux array of the spectrum in units of 1d-17 erg/s/cm^2/A [NPIX]
;   err0       -  Error array of the spectrum; same unites as flux0
;   z          -  Redshift of the quasar
;
; OPTIONAL INPUTS:
;   ivar0      -  Inverse variance array of the spectrum [NPIX]; will override err0
;   ra         -  RA of the object; required for dereddening
;   dec        -  Dec of the object; required for dereddening
;   deredden   -  /deredden to deredden the spectrum using Galactic reddening
;   emparfile  -  File path of the linefit Yanny parameter file
;   wave_range -  Trim the spectrum within [minwave,maxwave]
;   input_fitmask -  User-defined additional pixel mask [e.g., absorption mask]
;   add_noise  -  Set /add_noise to perturb the spectrum once using err array
;   fit_flag   -  A 5-element array to set fitting for different continuum components
;   outdir     -  Path of the output fits and QA plot files
;   output_name - Prefix of the output files 
;   psplot     -  Set /psplot to output a ps QA file to outdir/QA/
;   fits       -  Set /fits to write the output structure to outdir/fits/  
;   fit_line   -  If not set, then only fit the continuum
;   nsmooth    -  Number of pixels to boxcar-smooth for plotting
;
; OUTPUTS:
;   para       -  Output structure of the fitting results
;   linelist   -  List of lines
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
; 
;
; PROCEDURES CALLED:
; 
;
; Bugs:
;   1. auto_abs_rej iteration and outlier_mask tend to reject good pixels in high SN
;      spectra (primarily for Hbeta and Halpha fits). Fixed by excluding Hbeta/Halpha
;      regions in using auto_abs_rej and outlier_mask.
;   2. during outlier_mask to reject 5sigma lower pixels, some pixels near sharp narrow lines
;      (such as [NeV]) get rejected in high SN spectrum. Fixed by only apply outlier_mask to
;      lambda_rest<3000A
;
; REVISION HISTORY:
;   20-Feb-2014  Written by Yue Shen, Carnegie Obs
;   08-May-2014  Changed outlier_mask to apply to lambda_rest<3000A
;


function Fe_flux_mgii, xval, pp

   common Fe_temp, wave_Fe_civ, flux_fe_civ, wave_Fe_mgii, flux_Fe_mgii $
         , wave_Fe_balmer, flux_Fe_balmer

   if n_elements(wave_Fe_mgii) eq 0 then begin
     ; read in the Salvaider/Vestergaard template for MgII
     file = getenv('IDLRM_DIR')+'/template/feconv_uv_shen'
     ; wave_Fe1 [2200, 3090]; wave_Fe2 [4526, 6357]
     ; readcol, file, format = 'x,d,d,d,d', wave_1, flux_1, wave_2, flux_2, /silent
     readcol, file, format='d,d', logwave_1, flux_1
     wave_Fe_mgii = 10.D^logwave_1
     flux_Fe_mgii = flux_1*1d15
   endif

   ;c = 2.9979246d5
   Fe_FWHM = pp[1] ; broadening FWHM

   nnn = n_elements(xval)
   yval = dblarr(nnn)
   xval_new = xval*(1.0 + pp[2])

   ind = where(xval_new ge 1200. and xval_new le 3500.)
   if ind[0] ne -1 then begin

      if Fe_FWHM le 900.0 then sig_conv = sqrt(910.0^2 - 900.0^2)/2./sqrt(2.*alog(2.)) $
       else sig_conv = sqrt(Fe_FWHM^2 - 900.0^2)/2./sqrt(2.*alog(2.))

      ; Get sigma in pixel space
      sig_pix = sig_conv/103.6  ; where 550 km/s is the median dispersion
            ; for Salviander et al. (Vestergarrd's) Fe template from [2200, 3090]
            ; but note that the Vestergaard's template has a different median
            ; dispersion 55 km/s in wavelength for [2200,3090]
            ; For the Tsuzuki06 template, dispersion is ~50km/s within [2200,3500]
            ; and varies between 65km/s and 40km/s
            ; Here I am using feconv_uv_shen, which has constant v_disp=103.6 km/s

      khalfsz = round (4*sig_pix+1)
      xx= findgen(khalfsz*2+1) - khalfsz
      kernel = exp(-xx^2/(2*sig_pix^2))
      kernel = kernel/total(kernel)

      flux_Fe_conv = convol(flux_Fe_mgii, kernel, /center, /edge_truncate)
      yval[ind] = pp[0]*spline(wave_Fe_mgii, flux_Fe_conv, xval_new[ind])
   endif

   return, yval

end

function Fe_flux_balmer, xval, pp

   common Fe_temp, wave_Fe_civ, flux_fe_civ, wave_Fe_mgii, $
            flux_Fe_mgii, wave_Fe_balmer, flux_Fe_balmer

   if n_elements(wave_Fe_balmer) eq 0 then begin

    ; read in the Boroson&Green optical template for Balmer
      file = getenv('IDLRM_DIR')+'/template/irontemplate.dat'
      readcol, file, format = 'd,d', logwave_1, flux_1, /silent
      wave_Fe_balmer = 10.0D^logwave_1
      flux_Fe_balmer = flux_1*1d15
      ind = where(wave_Fe_balmer ge 3686. and wave_Fe_balmer le 7484.)
      wave_Fe_balmer = wave_Fe_balmer[ind]
      flux_Fe_balmer = flux_Fe_balmer[ind]

   endif


   Fe_FWHM = pp[1] ; broadening FWHM

   nnn = n_elements(xval)
   yval = dblarr(nnn)
   xval_new = xval*(1.0 + pp[2])

   ind = where(xval_new ge 3686. and xval_new le 7484.)
   if ind[0] ne -1 then begin
      ; Convolve the original Fe template with parameter Fe_FWHM 
      ; (I Zw1 has intrinsic FWHM 900 km/s)
      ; initial FWHM in the template is 900 km/s, 
      ; so force Fe_FWHM=910.0 if input Fe_FWHM <= 900 km/s
      if Fe_FWHM le 900.0 then sig_conv = sqrt(910.0^2 - 900.0^2)/2./sqrt(2.*alog(2.)) $
       else sig_conv = sqrt(Fe_FWHM^2 - 900.0^2)/2./sqrt(2.*alog(2.))  ; in km/s

      ; Get sigma in pixel space
      sig_pix = sig_conv/106.3 ; 106.3 km/s is the dispersion for the BG92 FeII template

      khalfsz = round (4*sig_pix+1)
      xx= findgen(khalfsz*2+1) - khalfsz
      kernel = exp(-xx^2/(2*sig_pix^2))
      kernel = kernel/total(kernel)
  
      flux_Fe_conv = convol(flux_Fe_balmer, kernel, /center, /edge_truncate)
      yval[ind] = pp[0]*spline(wave_Fe_balmer, flux_Fe_conv, xval_new[ind])
   endif

   return, yval

end

function f_poly_conti, xval, pp

  xval2 = xval - 3000.
  yval = 0*xval2
  for i=0L, n_elements(pp) - 1L do $
    yval = yval + pp[i]*xval2^(i+1)

  return, yval
end

function f_conti_only, xval, pp

   f_pl = pp[0]*(xval/3000.0)^pp[1]  ; power-law continuum
   f_conti_BC = balmer_conti(xval, pp[2:4])  ; Balmer continuum
   f_poly = f_poly_conti(xval, pp[5:*])
   yval = f_pl + f_conti_BC + f_poly

   return, yval
end

function f_fe_only, xval, pp

   yval = Fe_flux_mgii(xval, pp[0:2]) + Fe_flux_balmer(xval, pp[3:5])

   return, yval
end

function f_conti_all, xval, pp

   ; pp[0]: norm_factor for the MgII Fe_template
   ; pp[1]: FWHM for the MgII Fe_template
   ; pp[2]: small shift of wavelength for the MgII Fe template
   ; pp[3:5]: same as pp[0:2] but for the Hbeta/Halpha Fe template
   ; pp[6]: norm_factor for continuum f_lambda = (lambda/3000.0)^{-alpha}
   ; pp[7]: slope for the power-law continuum
   ; pp[8:10]: norm, Te and Tau_e for the Balmer continuum at <3646 A
   ; pp[11:*]: polynomial for the continuum

   ;common objdata, ra_c, dec_c, z_c

   nnn = n_elements(xval)
   yval = dblarr(nnn)
   f_Fe_MgII = Fe_flux_mgii(xval, pp[0:2])
   f_Fe_Balmer = Fe_flux_balmer(xval, pp[3:5])

   f_pl = pp[6]*(xval/3000.0)^pp[7]  ; power-law continuum
   f_conti_BC = balmer_conti(xval, pp[8:10])  ; Balmer continuum
   f_poly = f_poly_conti(xval, pp[11:*])
   ;f_pl = pp[6]*(xval/3000.0)^pp[7]+pp[8]*(xval/3000.0)^pp[9] ;double PL continuum
   ;f_conti_BC = balmer_conti(xval, pp[10:12])  ; Balmer continuum
   ;f_poly = f_poly_conti(xval, pp[13:*])

   yval = f_pl + f_Fe_MgII + f_Fe_Balmer + f_conti_BC + f_poly

   return, yval

end

;-----------------------------------------------------------------------------------
;   Main Routine 
;-----------------------------------------------------------------------------------

pro rm_qsofit, lam0, flux0, err0, z, ivar0=ivar0,ra=ra,dec=dec, deredden=deredden $
       , emparfile=emparfile, conti_fit=conti_fit, fit_flag=fit_flag $
       , line_fit_flag=line_fit_flag $
       , f_conti_model=f_conti_model, xrange=xrange, outdir=outdir $
       , psplot=psplot,fits=fits,output_name=output_name,silent=silent $
       , fit_line=fit_line, para=para, SDSS_name=SDSS_name,objtag=objtag,append=append $
       , linelist=linelist, wave_range=wave_range $
       , input_fitmask=input_fitmask $ ; input fitmasks, e.g., custom absorption masks
       , poly_ord=poly_ord, rej_iter=rej_iter $
       , add_noise = add_noise, local_CIV_CIII_fit = local_CIV_CIII_fit $
       , auto_abs_rej = auto_abs_rej, conti_auto_abs_rej = conti_auto_abs_rej $
       , nsmooth=nsmooth, diet=diet, sim_conti_fe_fit=sim_conti_fe_fit $
       , noplot=noplot, get_line_prop=get_line_prop $
       , xtol=xtol,ftol=ftol, subdir=subdir, more_anno=more_anno

   ;on_error, 2

   ;common objdata, ra_c, dec_c, z_c
   ;ra_c = ra & dec_c = dec & z_c = z

   ; Define line fitting tolerance in mpfitfun, 
   ; which helps for some objects in line fits (i.e., Hbeta and Halpha)
   if ~keyword_set(xtol) then xtol=1d-20
   if ~keyword_set(ftol) then ftol=1d-15

   ; Define constants
   cs = 2.9979246d5  ; speed of light, km/s
   if n_elements(z) eq 0 then z = 0. ; default is in restframe
   ; if n_elements(conti_auto_abs_rej) eq 0 then conti_auto_abs_rej = 1L
   
   ; default is to auto reject absorption from the first linefit
   if n_elements(auto_abs_rej) eq 0 then auto_abs_rej = 1L 
   if n_elements(fit_line) eq 0 then fit_line = 1L ; default is to fit the lines

   lam = lam0 & flux = flux0
   if keyword_set(ivar0) then begin
      err = 0.*ivar0
      ind = where(ivar0 NE 0.)
      err[ind] = 1./sqrt(ivar0[ind])
   endif else err = err0
   
   npix = n_elements(flux)
   
   if keyword_set(add_noise) then begin ; add gaussian noise
      flux = flux + randomn(seed, npix)*err
   endif

   ; fit_flag=[fit_balmer_conti,not_fit_UV_Fe,fix_mgii_fe_FWHM,
   ;           fix_balmer_fe_FWHM,fit_poly_conti]
   ; Set the flag to fit different components
   ; Note in some cases the Balmer continuum is poorly constrained, 
   ; and hence will degrade the global power-law fit
   if n_elements(fit_flag) eq 0 then fit_flag = [0,0,0,0,1]
   
   ; If not fitting UV FeII template, then do not set fix_mgII_fe_FWHM=1
   if fit_flag[1] eq 1 then fit_flag[2] = 0L 
   
   ; Default is to add a 3rd order polynomial to the continuum
   if n_elements(poly_ord) eq 0 then poly_ord = 3

   ; Number of linefit absorption rejection iteration
   if n_elements(rej_iter) eq 0 then rej_iter = 2L

   ; Default is to fit continuum+FeII simultaneously
   if n_elements(sim_conti_fe_fit) eq 0 then sim_conti_fe_fit=1L

   ; common block for the FeII template
   common Fe_temp, wave_Fe_civ, flux_fe_civ, wave_Fe_mgii, flux_Fe_mgii, $
            wave_Fe_balmer, flux_Fe_balmer
   if n_elements(wave_Fe_mgii) eq 0 then begin

      ; read in the Salvaider/Vestergaard template for MgII
      file = getenv('IDLRM_DIR')+'/template/feconv_uv_shen'
      ; wave_Fe1 [2200, 3090]; wave_Fe2 [4526, 6357]
      ; readcol, file, format = 'x,d,d,d,d', wave_1, flux_1, wave_2, flux_2, /silent
      readcol, file, format='d,d', logwave_1, flux_1
      wave_Fe_mgii = 10.D^logwave_1
      flux_Fe_mgii = flux_1*1d15

      ; read in the Boroson&Green optical template for Balmer
      file = getenv('IDLRM_DIR')+'/template/irontemplate.dat'
      readcol, file, format = 'd,d', logwave_1, flux_1, /silent
      wave_Fe_balmer = 10.0D^logwave_1
      flux_Fe_balmer = flux_1*1d15
      ind = where(wave_Fe_balmer ge 3686. and wave_Fe_balmer le 7484.)
      wave_Fe_balmer = wave_Fe_balmer[ind]
      flux_Fe_balmer = flux_Fe_balmer[ind]

   endif

   ; deredden Galactic extinction if asked
   if keyword_set(deredden) then begin
      if not keyword_Set(silent) then splog, $
         'Dereddening spectrum using the SFD map and CCM extinction curve'
      dereddening_spec, lam, flux, err=err, ra=ra, dec=dec $
         , dered_flux = dered_flux, dered_err = dered_err
      flux = dered_flux & err = dered_err
   endif

   ivar = dblarr(n_elements(err))
   indd = where(err gt 1d-6)
   ivar[indd] = 1./(err[indd])^2

   ; mask out 3-sigma lower outlier of the 20-pix smoothed spectrum, 
   ; to reduce the effects of absorption
   fitmask = ivar NE 0.0

   ; Mask outliers and absorption in the emission lines by finding 5-sigma
   ; lower outliers from 20 pix smoothed spectrum
   if n_elements(outlier_mask) eq 0 then outlier_mask = 1L
   if keyword_set(outlier_mask) then begin
      disp_vec = (alog10(lam)-shift(alog10(lam), 1))[1:*]
      disp = djs_median(disp_vec)
      NRES = 20L   ; 20L
      bkspace = NRES*disp
      spec_set = bspline_iterfit(alog10(lam), flux $
                             , invvar = ivar*fitmask $
                             , bkspace = bkspace $
                             , yfit = flux_spline, maxiter = 10 $
                             , upper = 0, lower = 5, nord = 3 $
                             , maxrej = 50, outmask = outmask, /silent)
      ; exclude redward of restframe 3000A, since little absorption expected there.
      ind_exclude = where(lam/(1.+z) ge 3000.)
      if ind_exclude[0] ne -1 then outmask[ind_exclude]=1
      fitmask = fitmask*outmask
   endif
   ; Let's mask the last 40 pixels in the red as bad
   ; fitmask[npix-40:*]=0
   
   ; Add additional mask bits from input
   if keyword_set(input_fitmask) then begin
      if n_elements(input_fitmask) eq n_elements(lam) then $
         fitmask = fitmask*input_fitmask else $
         splog, 'input_fitmask dimension does not match lam grid.'
   endif

   ; Define the continuum+Fe+Balmer continuum fitting windows
   window0 = [1350., 1360.]
   window1 = [1445., 1465.] & window2 = [1700., 1705.]
   window3_1 = [2155., 2400.] & window3_2 = [2480., 2675.]
   window4 = [2925., 3400.] & window5 = [4200., 4230.]
   window6 = [4435., 4700.] & window7 = [5100., 5535.]
   window8 = [6005., 6035.] & window9 = [6035., 6250.]
   window10 = [6800., 7000.] & window11 = [7160., 7180.]

   ; This is the continuum windows that are used in conti_fit.
   ; You may add new windows here
   window_all = [  [1150., 1170.], [1275., 1290.], [1350., 1360.], [1445., 1465.], $
                   [1700., 1705.], [1770., 1800.], [2155., 2400.], [2480., 2675.], $
                   [2925., 3400.], [3775., 3832.], [4000., 4050.], [4200., 4230.], $  ; [3775., 3832.]
                   [4435., 4700.], [5100., 5535.], [6005., 6035.], [6110., 6250.], $
                   [6800., 7000.], [7160., 7180.], [7500., 7800.], [8050., 8150.] ]

   ; Shift to restframe
   wave = lam/(1.0+z)

   ; Trim the spectrum in the windows specified by keywords wave_range
   if keyword_Set(wave_range) then begin
      ind = where(wave ge wave_range[0] and wave le wave_range[1])
      wave = wave[ind] & flux = flux[ind] & err = err[ind]
      fitmask = fitmask[ind]
   endif

   ; This is the pixels for simultaneous continuum+FeII fit
   if n_elements(ind) ne 0 then tmp_ind = temporary(ind)
   ;print, n_elements(ind)
   for jj=0L, (size(window_all))[2] - 1 do begin
      tmp = where( wave ge window_all[0,jj] and wave le window_all[1,jj] $
                  and fitmask ne 0., ntmp)
      if ntmp gt 0 then begin
         if n_elements(ind) eq 0 then ind = tmp else ind = [ind, tmp]
      endif
   endfor

   if n_elements(ind) lt 10 then begin
      splog, 'Fitting pixel < 10. return. '
      return
   endif

   ; ## now fit the psedu-continuum
   ; set limiting conditions on paramters
   parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0],tied:''}, $
              11+poly_ord)
   parinfo.value = [1.0D, 3000.0D, 0.0, 1.0, 3000.0D, 0.0, $
                    1.0D, -2.0D, $
                    1.0, 15000.D, 0.5, $
                    replicate(0.,poly_ord) ]

   ;-- limit the Fe_MgII component
   parinfo[0].limited = [1, 0] & parinfo[0].limits = [0., 1d10]
   parinfo[1].LIMITED = [1, 1] & parinfo[1].LIMITS = [1200, 1d4]
   parinfo[2].LIMITED = [1, 1] & parinfo[2].LIMITS = [-0.01, 0.01]
   if fit_flag[1] eq 1 then begin ; do not fit the UV iron template
      parinfo[0].limited = [0,0] & parinfo[0].value = 0. 
      parinfo[0:2].fixed = 1L
   endif
   if fit_flag[2] eq 1 then begin
      parinfo[1].fixed = 1L ; fix the FWHM of the MgII FeII
   endif
   ;-- limit the Fe_balmer component
   parinfo[3].limited = [1, 0] & parinfo[3].limits = [0., 1d10]
   parinfo[4].LIMITED = [1, 1] & parinfo[4].LIMITS = [1200, 1d4]
   parinfo[5].LIMITED = [1, 1] & parinfo[5].LIMITS = [-0.01, 0.01]
   if fit_flag[3] eq 1 then begin
      parinfo[4].fixed = 1L ; fix the FWHM of the Balmer FeII
   endif
   ;-- limit the PL component
   parinfo[6].limited = [1,0] & parinfo[6].limits = [0., 1d10]
   parinfo[7].LIMITED = [1, 1] & parinfo[7].LIMITS = [-5.0,3.0]
   ;parinfo[8].limited = [1,0] & parinfo[8].limits = [0., 1d10]
   ;parinfo[9].LIMITED = [1, 1] & parinfo[9].LIMITS = [-5.0,3.0]
   ;-- limit the Balmer continuum
   if fit_flag[0] eq 0 then begin ; do not fit for Balmer continuum
       parinfo[8:10].fixed = 1L & parinfo[8].value = 0.
      ;parinfo[10:12].fixed = 1L & parinfo[10].value = 0.
   endif else begin
      parinfo[8].limited = [1, 0] & parinfo[8].limits = [0., 1d10]
      parinfo[9].LIMITED = [1, 1] & parinfo[9].LIMITS = [10000., 50000.]
      parinfo[10].LIMITED = [1, 1] & parinfo[10].LIMITS = [0.1, 2.]
      ;parinfo[10].limited = [1, 0] & parinfo[10].limits = [0., 1d10]
      ;parinfo[11].LIMITED = [1, 1] & parinfo[11].LIMITS = [10000., 50000.]
      ;parinfo[12].LIMITED = [1, 1] & parinfo[12].LIMITS = [0.1, 2.]
   endelse
   ;-- limit the polynomial continuum component
   if fit_flag[4] eq 0 then begin
      parinfo[11:*].value = 0. & parinfo[11:*].fixed = 1L
      ;parinfo[13:*].value = 0. & parinfo[13:*].fixed = 1L
   endif
   ;-- enforce the UV FeII to be zero is no enough pixels near MgII
   ind_cover = where( ( (wave ge window3_2[0] and wave le window3_2[1]) or $
                        (wave ge window4[0] and wave le window4[1]) ) $
                     and fitmask ne 0, n_cover)
   if n_cover lt 100 then begin
      parinfo[0].value = 0. & parinfo[0:2].fixed = 1L
   endif
   ;-- enforce the optical FeII to be zero if not covered
   ind_cover = where( ( (wave ge window6[0] and wave le window6[1]) or $
                        (wave ge window7[0] and wave le window7[1]) ) $
                     and fitmask ne 0, n_cover)
   if n_cover lt 100 then begin
      parinfo[3].value = 0. & parinfo[3:5].fixed = 1L
   endif

   ; First try to fit the continuum-alone
   if n_elements(ind_conti_only) ge 50 and not keyword_set(sim_conti_fe_fit) then begin
      
      parinfo1 = parinfo[6:*]
      conti_only_fit = mpfitfun('f_conti_only', wave[ind_conti_only] $
                     , flux[ind_conti_only], err[ind_conti_only] $
                     , parinfo = parinfo1, perror = perror, yfit = yfit1, /quiet $
                     , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
      conti_only_err = perror
      ttt = where(parinfo1.fixed eq 1 or strlen(parinfo1.tied) gt 0, nfix)
      dof = n_elements(ind_conti_only) - (n_elements(parinfo1) - nfix)
      redchi2_conti_only = chi2/dof

      ; If this is successful, then proceed to fit the FeII
      if redchi2_conti_only gt 0 and redchi2_conti_only lt 100 then begin

         parinfo1 = parinfo[0:5]
         flux1 = flux[ind_Fe_only] - f_conti_only(wave[ind_Fe_only], conti_only_fit)
         fe_only_fit = mpfitfun('f_fe_only', wave[ind_Fe_only], flux1 $
                     , err[ind_Fe_only] $
                     , parinfo = parinfo1, perror = perror, yfit = yfit1, /quiet $
                     , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
         Fe_only_err = perror
         ttt = where(parinfo1.fixed eq 1 or strlen(parinfo1.tied) gt 0, nfix)
         dof = n_elements(ind_Fe_only) - (n_elements(parinfo1) - nfix)
         redchi2_Fe_only = chi2/dof

         if redchi2_Fe_only gt 0 and redchi2_Fe_only le 100 then begin
             conti_fit = [conti_only_fit, Fe_only_fit]
             conti_err = [conti_only_err, Fe_only_err]
         endif

      endif
   endif

   if n_elements(ind_conti_only) le 50 or keyword_set(sim_conti_fe_fit) then begin

      conti_fit = mpfitfun('f_conti_all', wave[ind], flux[ind], err[ind] $
                        , parinfo = parinfo, perror = perror, yfit = yfit1, /quiet $
                        , nfev = nfev, niter = niter, status = status, bestnorm = chi2 $
                        , MAXITER = 100L)
      conti_err = perror
      ttt = where(parinfo.fixed eq 1 or strlen(parinfo.tied) gt 0, nfix)
      dof = n_elements(ind) - (n_elements(parinfo) - nfix)
      redchi2 = chi2/dof
      if not keyword_set(silent) then splog, $
         'CONTI-FITTING: MPFIT nfev=', nfev, ' niter=', niter, $
          ' status=', status, ' redchi2=', redchi2

      ; Re-fit the continuum without the Balmer continuum,
      ; if the fit is better, replace the original fit
      ; only do this if fit_Balmer_continuum is set to 1
      ;if fit_flag[0] eq 1 then begin
      ;   parinfo[8:10].fixed = 1L & parinfo[8].value = 0. 
      ;
      ;   conti_fit2 = mpfitfun('f_conti_all', wave[ind], flux[ind], err[ind] $
      ;            , parinfo = parinfo, perror = perror2, yfit = yfit2, /quiet $
      ;            , nfev = nfev, niter = niter, status = status2, bestnorm = chi2)
      ;   dof = n_elements(ind) - (n_elements(parinfo) - nfix)
      ;   redchi2_2 = chi2/dof
      ;   if not keyword_set(silent) then splog, $
      ;     'CONTI-FITTING: MPFIT nfev=', nfev, ' niter=', niter, $
      ;     ' status=', status2, ' redchi2=', redchi2_2
      ;   if status2 gt 0 and redchi2_2 gt 0 and redchi2_2 lt redchi2 then begin
      ;      conti_fit = conti_fit2 & conti_err = perror2 & status = status2 
      ;      redchi2 = redchi2_2 & yfit1 = yfit2
      ;      if not keyword_set(silent) then splog, $
      ;       'CONTI-FITTING: Replace fit with NO BALMER CONTINUUM'
      ;   endif
      ;endif

      ; Perform one iteration to remove 3sigma pixel below the first continuum fit
      ; this is to avoid situations where one of the continuum windows falls within
      ; a BAL trough
      if keyword_Set(conti_auto_abs_rej) and redchi2 lt 100. then begin
         fitmask2 =  replicate(1L, n_elements(wave))
         ind_abs = where(flux[ind] lt yfit1 - 3.*err[ind] and wave[ind] lt 3500.)
         if ind_abs[0] ne -1 then fitmask2[ind[ind_abs]] = 0L
         fitmask2 = fitmask*fitmask2
         if n_elements(ind2) ne 0 then ind_tmp = temporary(ind2)
         for jj=0L, (size(window_all))[2] - 1 do begin
            tmp = where( wave ge window_all[0,jj] and wave le window_all[1,jj] $
                     and fitmask2 ne 0., ntmp)
            if ntmp gt 0 then begin
               if n_elements(ind2) eq 0 then ind2 = tmp else ind2 = [ind2, tmp]
            endif
         endfor

         if (n_elements(ind2) gt 10L) and (n_elements(ind2) lt n_elements(ind)) then begin
            conti_fit2 = mpfitfun('f_conti_all', wave[ind2], flux[ind2], err[ind2] $
                     , parinfo = parinfo, perror = perror2, yfit = yfit2, /quiet $
                     , nfev = nfev, niter = niter, status = status2, bestnorm = chi2 $
                     , MAXITER = 100L)
            ttt = where(parinfo.fixed eq 1 or strlen(parinfo.tied) gt 0, nfix)
            dof = n_elements(ind2) - (n_elements(parinfo) - nfix)
            redchi2_2 = chi2/dof
            if status2 gt 0 and redchi2_2 gt 0 and redchi2_2 lt redchi2 then begin
               if not keyword_set(silent) then splog, $
                'CONTI-FITTING: MPFIT nfev=', nfev, ' niter=', niter, $
                ' status=', status2, ' redchi2=', redchi2_2
               conti_fit = conti_fit2 & conti_err = perror2 & status = status2 
               redchi2 = redchi2_2 & yfit1 = yfit2
               ind = ind2 & fitmask = fitmask2
               if not keyword_set(silent) then splog, $
                'CONTI-FITTING: Replace fit with absorption rejection'
            endif
         endif
      endif

   endif

   ; Output the fitting results
   output = create_struct('z',z,'wave', wave, 'flux', flux, 'err', err, 'fitmask', fitmask, $
         'ind_conti_fit', ind, 'conti_fit', conti_fit, 'conti_fit_err', conti_err, $
         'conti_redchi2', redchi2, 'conti_status', status, 'conti_npix_fit', n_elements(ind))
   ; output = struct_addtags(specdata, output)

   f_fe_mgii_model = fe_flux_mgii(wave, conti_fit[0:2])
   f_fe_balmer_model = fe_flux_balmer(wave, conti_fit[3:5])
   f_pl_model = conti_fit[6]*(wave/3000.0)^conti_fit[7]
   f_bc_model = balmer_conti(wave, conti_fit[8:10])
   f_poly_model = f_poly_conti(wave, conti_fit[11:*]) 
   ;f_pl_model=conti_fit[6]*(wave/3000.0)^conti_fit[7] $
   ;  + conti_fit[8]*(wave/3000.0)^conti_fit[9]
   ;f_bc_model = balmer_conti(wave, conti_fit[10:12])
   ;f_poly_model = f_poly_conti(wave, conti_fit[13:*])
   f_conti_model = f_pl_model + f_fe_mgii_model + f_fe_balmer_model $
                   + f_bc_model + f_poly_model

   if keyword_set(local_CIV_CIII_fit) then begin ; locally fit CIV+CIII] in [1000, 2165]

      for jj=0L, (size(window_all))[2] - 1 do begin
         tmp = where( wave ge window_all[0,jj] and wave le window_all[1,jj] $
                     and fitmask ne 0. and wave le 2165., ntmp)
         if ntmp gt 0 then begin
            if n_elements(ind3) eq 0 then ind3 = tmp else ind3 = [ind3, tmp]
         endif
      endfor
      ; first turn off Balmer continuum and UV iron template fit
      parinfo[8].fixed = 1L & parinfo[8].value = 0.
      parinfo[9].fixed = 1L & parinfo[10].fixed = 1L  
      parinfo[2].limited = [0,0] & parinfo[2].value = 0. & parinfo[2].fixed = 1L
      parinfo[3].fixed = 1L & parinfo[4].fixed = 1L
      conti_fit3 = mpfitfun('f_conti_all', wave[ind3], flux[ind3], err[ind3] $
                 , parinfo = parinfo, perror = perror3, yfit = yfit3, /quiet $
                 , nfev = nfev, niter = niter, status = status3, bestnorm = chi2)
      if not keyword_set(silent) then splog, 'CONTI-FITTING: MPFIT nfev=', nfev, $
        ' niter=', niter, ' status=', status3
      ; update the psedu-continuum model, connected at 2165A
      ind_local = where(wave le 2165.)
      f_pl_model[ind_local] = conti_fit3[6]*(wave[ind_local]/3000.0)^conti_fit3[7]
      f_fe_mgii_model[ind_local] = 0.
      f_fe_balmer_model[ind_local] = 0.
      f_bc_model[ind_local] = 0.
      f_conti_model = f_pl_model + f_fe_mgii_model + f_fe_balmer_model $
                      + f_bc_model + f_poly_model

      temp = create_struct('conti_fit3', conti_fit3, 'conti_fit3_err', perror3)
      output = struct_addtags(output, temp)
   endif

   ; ----------------
   ; Start to fit the emission lines: CIV, CIII], MgII, Hbeta, Halpha
   ; first subtract the psedu-continuum to get emission line fluxs
   line_flux = flux - f_conti_model
   if keyword_set(local_CIV_CIII_fit) then $
     line_flux[ind_local] = flux[ind_local] - f_conti_all(wave[ind_local], conti_fit3)
   output = struct_addtags(output, {line_flux:line_flux})

   ; Setup line fitting parameters
   if keyword_set(emparfile) then emline_file=emparfile else $
    emline_file = getenv('IDLRM_DIR')+'/etc/qsoline.par'
   yanny_read, emline_file, pdata
   linelist = *pdata[0]
   yanny_free, pdata
   ngauss_tot = total(linelist.ngauss)
   parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0],tied:''}, $
          ngauss_tot*3)
   npar = n_elements(parinfo)
   line_fit_all = dblarr(npar) & perror_all = dblarr(npar)
   nline = n_elements(linelist)
   ;ncomp = n_elements(uniq(linelist.compname,/sort))
   linecomp = linelist.compname & linename = linelist.linename
   parindx = lonarr(2,nline)
   parfitflag = lonarr(npar)
   parlinename = replicate('',npar)
   line_redchi2 = dblarr(npar) & line_status = lonarr(npar)
   line_npix_fit = lonarr(npar)
   inext = 0
   ; assign initial values and limitations
   for i=0L, nline - 1L do begin         
      ngauss = linelist[i].ngauss
      parindx[*,i] = [inext, inext+3*ngauss-1]
     
      for j=0L, ngauss - 1L do begin

         ; add a small, random offset of each broad gaussian component
         if ngauss gt 1 then $
          lnlamini = alog(linelist[i].lambda) + (2*randomu(seed,1)-1.)*1d-3 else $
          lnlamini = alog(linelist[i].lambda)
         parlinename[inext:inext+2] = linelist[i].linename
         parinfo[inext:inext+2].value = [linelist[i].fvalue, lnlamini $
                , linelist[i].inisig]
         parinfo[inext].limited = [1,0] & parinfo[inext].limits=[0,1d10]
         parinfo[inext+1].limited = [1,1]
         parinfo[inext+1].limits=[lnlamini-linelist[i].voff, lnlamini+linelist[i].voff]
         parinfo[inext+2].limited = [1,1]
         parinfo[inext+2].limits = [linelist[i].minsig, linelist[i].maxsig]
         
         inext = inext + 3
      endfor
   endfor

   ; Only fit the line if the keyword FIT_LINE is set
   minwave = min(wave,max=maxwave)
   if keyword_set(fit_line) then begin

      ; Find which lines are covered
      ind1 = where(linelist.lambda ge minwave and linelist.lambda le maxwave, $
            complement = indd1)
      if indd1[0] ne -1 then linelist[indd1].fitflag = 0
      ind_fit = where(linelist.fitflag eq 1)
      if ind_fit[0] eq -1 then $
        splog, 'No line to fit, return' $
      else begin
         
         ; loop over each line complex
         uniq_linecomp = uniq(linelist[ind_fit].compname)
         ncomp = n_elements(uniq_linecomp)
         compname = (linelist[ind_fit].compname)[uniq_linecomp]

         for ii=0L, ncomp - 1L do begin
             ind_line = where(linelist.compname eq compname[ii] and $
                linelist.fitflag eq 1, nline_fit)
             linelist_fit = linelist[ind_line]
             ngauss_fit = linelist_fit.ngauss
             if not keyword_set(silent) then splog, 'Fit line complex: ', compname[ii]
             iline_start = ind_line[0] & iline_end = ind_line[nline_fit-1]
             parindx_sub = parindx[0,iline_start]+indgen(3*linelist[iline_start].ngauss)
             for iitmp=1, nline_fit - 1L do begin
                 parindx_sub = [parindx_sub, $
                  parindx[0,ind_line[iitmp]]+indgen(3*linelist[ind_line[iitmp]].ngauss)]
             endfor

             parinfo1 = parinfo[ parindx_sub ]
             ; Now setup constraints of parameters
             ; Tie velocity
             arr = linelist_fit.vindex
             arr = arr[sort(arr)]
             arr = arr[uniq(arr)]
             for iuniq=0L, n_elements(arr) - 1L do begin
                ind_fix = where(linelist_fit.vindex eq arr[iuniq],ntmp)
                if ntmp gt 0 and arr[iuniq] ne 0 then begin
                   
                   ; For a broad line, enforce symmetric profile if required
                   if ntmp eq 1 and linelist_fit[ind_fix[0]].ngauss gt 1 then begin
                      ind_fix = ind_fix[0]
                      pstr0='P(' + string(1+3*( total(ngauss_fit[0:ind_fix]) $
                       - ngauss_fit[ind_fix]),format='(i0)')+')'
                      for igauss=1, ngauss_fit[ind_fix]-1L do begin
                        parinfo1[1+3*(total(ngauss_fit[0:ind_fix]) $
                         - ngauss_fit[ind_fix])+igauss*3].tied = pstr0
                      endfor
                   endif

                   pstr0='(P(' + string(1+3*( total(ngauss_fit[0:ind_fix[0]])-1) $
                             , format='(i0)')+')'
                   for iitmp=1,ntmp-1 do begin
                      parinfo1[1+3*(total(ngauss_fit[0:ind_fix[iitmp]])-1)].tied=pstr0 $
                       + ' - alog('+string(linelist_fit[ind_fix[0]].lambda,format='(f7.2)') $
                       + ')) + alog('+string(linelist_fit[ind_fix[iitmp]].lambda,format='(f7.2)') $
                       + ')'
                   endfor
                endif
             endfor
             ; Tie line width
             arr = linelist_fit.windex
             arr = arr[sort(arr)]
             arr = arr[uniq(arr)]
             for iuniq=0L, n_elements(arr) - 1L do begin
                ind_fix = where(linelist_fit.windex eq arr[iuniq],ntmp)
                if ntmp gt 1 and arr[iuniq] ne 0 then begin
                   pstr0='P('+string(2+3*(total(ngauss_fit[0:ind_fix[0]])-1) $
                    ,format='(i0)')+')'
                   for iitmp=1,ntmp-1 do begin
                      parinfo1[2+3*( total(ngauss_fit[0:ind_fix[iitmp]])-1)].tied = pstr0
                   endfor
                endif
             endfor
             ; Tie line flux
             arr = linelist_fit.findex
             arr = arr[sort(arr)]
             arr = arr[uniq(arr)]
             for iuniq=0L, n_elements(arr) - 1L do begin
                ind_fix = where(linelist_fit.findex eq arr[iuniq],ntmp)
                if ntmp gt 1 and arr[iuniq] ne 0 then begin

                   pstr0='P('+string(3*( total(ngauss_fit[0:ind_fix[0]])-1) $
                    ,format='(i0)')+')'
                   for iitmp=1,ntmp-1 do begin
                      parinfo1[3*(total(ngauss_fit[0:ind_fix[iitmp]])-1)].tied=pstr0 $
                       + '*' + string(linelist_fit[ind_fix[0]].lambda,format='(f7.2)') $
                       + '/' + string(linelist_fit[ind_fix[iitmp]].lambda,format='(f7.2)') $
                       + '*' + string(linelist_fit[ind_fix[iitmp]].fvalue/linelist_fit[ind_fix[0]].fvalue $
                       , format='(f3.1)')
                   endfor
                endif
             endfor

             ;if compname[ii] eq 'CIV' then message, 'stop!'

             ; Do the fit
             ind = where(wave ge linelist_fit[0].minwav $
                     and wave le linelist_fit[0].maxwav $
                     and fitmask ne 0., nnn)           
             if nnn gt 10 then begin
                line_fit = mpfitfun('manygauss',alog(wave[ind]),line_flux[ind],err[ind] $
                     , MAXITER = 500L, parinfo = parinfo1, xtol=xtol, ftol=ftol $
                     , perror = perror, yfit = yfit2, /quiet $
                     , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
                ttt = where(parinfo1.fixed eq 1 or strlen(parinfo1.tied) gt 0, nfix)
                dof = n_elements(ind) - (n_elements(parinfo1) - nfix)
                redchi2 = chi2/dof
                if not keyword_set(silent) then splog,'LINE-FITTING: MPFIT nfev=',nfev $
                    , ' niter=', niter, ' status=', status, ' redchi2=', redchi2

                ; Do a further step to remove 3sigma absorption below the first fit
                ; Do not do this for the Halpha or Hbeta line complex
                ; first setup fitmask2 and the new line fitting windows
                if keyword_set(auto_abs_rej) and redchi2 lt 100. $
                   and compname[ii] ne 'Halpha' and compname[ii] ne 'Hbeta' then begin
                
                   n_iter = 1L

          iter1:
                   fitmask2 =  replicate(1L, n_elements(wave))
                   ind_abs = where(line_flux[ind] lt yfit2 - 3.*err[ind])
                   if ind_abs[0] ne -1 then fitmask2[ind[ind_abs]] = 0L
                   fitmask2 = fitmask*fitmask2

                   ind2 = where(wave ge linelist_fit[0].minwav $
                            and wave le linelist_fit[0].maxwav $
                            and fitmask2 ne 0.)

                   if (n_elements(ind2) gt 10L) and (n_elements(ind2) lt n_elements(ind)) then begin

                      line_fit2 = mpfitfun('manygauss', alog(wave[ind2]) $
                           , line_flux[ind2], err[ind2] $
                           , MAXITER = 500L, parinfo = parinfo1, xtol=xtol,ftol=ftol $
                           , perror = perror2, yfit = yfit2_2, /quiet $
                           , nfev = nfev, niter = niter, status = status2, bestnorm = chi2_2)
                      ttt = where(parinfo1.fixed eq 1 or strlen(parinfo1.tied) gt 0, nfix)
                      dof2 = n_elements(ind2) - (n_elements(parinfo1) - nfix)
                      redchi2_2 = chi2_2/dof2
                      if not keyword_set(silent) then splog $
                        , 'LINE-FITTING: MPFIT nfev=', nfev $
                        , ' niter=', niter, ' status=', status2, ' redchi2=', redchi2_2 $
                        , ' iter=', n_iter
                      ; replace the original fit if justified
                      if status2 gt 0 and redchi2_2 gt 0 and redchi2_2 lt redchi2 then begin
                         if not keyword_set(silent) then $
                          splog, 'Improved fits by rejecting absorptions. '
                         ind = ind2 & line_fit = line_fit2 & perror = perror2
                         redchi2 = redchi2_2 & status = status2 & yfit2 = yfit2_2
                         fitmask = fitmask2
                         if n_iter lt rej_iter then begin
                            n_iter = n_iter + 1L
                            goto, iter1
                         endif
                      endif

                   endif
                endif

                if n_elements(ind_line_fit) eq 0 then ind_line_fit=ind else $
                   ind_line_fit = [ind_line_fit, ind]
                line_fit_all[ parindx_sub ] = line_fit
                perror_all[ parindx_sub ] = perror
                parfitflag[parindx_sub] = 1L
                line_redchi2[parindx_sub] = redchi2
                line_status[parindx_sub] = status
                line_npix_fit[parindx_sub] = n_elements(ind)

            endif
         endfor
         output.fitmask = fitmask
         temp = create_struct('linename',parlinename, 'fitflag', parfitflag $
           , 'ind_line_fit', n_elements(ind_line_fit) gt 0 ? ind_line_fit : [-1L] $
           , 'line_fit', line_fit_all $
           , 'line_fit_err', perror_all, 'line_redchi2', line_redchi2 $
           , 'line_status', line_status, 'line_npix_fit', line_npix_fit)
         output = struct_addtags(output, temp)
      endelse
 
   endif
   ; --------------------------- End of fitting emission lines ---------------

   if keyword_set(diet) then begin
      tag_rej=['WAVE','FLUX','ERR','FITMASK','IND_CONTI_FIT','LINE_FLUX','IND_LINE_FIT']
      remove_tags, output, tag_rej, output1
      output = output1
   endif
   para = output

   ; Write the fit results
   if n_elements(outdir) eq 0 then cd, current=outdir
   if n_elements(output_name) eq 0 then output_name = 'qsofit'
   if keyword_Set(fits) then begin ; output the fits result
      if ~keyword_set(subdir) then fitsfile = outdir + '/fits/' + output_name + '.fits' $
        else begin 
          fitsfile = outdir + '/fits/' + subdir + '/' + output_name + '.fits'
          if file_test(outdir + '/fits/' + subdir) eq 0 then spawn, 'mkdir ' + outdir + '/fits/' + subdir
        endelse
      mwrfits, output, fitsfile, /create
      spawn, 'gzip -f ' + fitsfile
   endif

   if keyword_set(noplot) then return

   ; ------------------------  Plotting block ----------------------------
   ; now make a plot
   if keyword_set(psplot) then begin
      if not keyword_set(append) then begin
         if ~keyword_set(subdir) then figfile = outdir + '/QA/' + output_name + '.ps' $
            else begin 
               figfile = outdir + '/QA/' + subdir + '/' + output_name + '.ps'
               if file_test(outdir + '/QA/' + subdir) eq 0 then spawn, 'mkdir ' + outdir + '/QA/' + subdir
            endelse
         ;begplot, name=figfile, xsize = 40, ysize = 25,/color
         begplot, name=figfile, /landscape, /color
      endif
      charsize = 1. & thick = 4. & xticks = 2L & xminor = 5L
      linethick = 0.1 & symsize = 3.
      ang = string(197B)
   endif else begin
      linethick = 1. & symsize = 3.
      xticks = 2L & xminor = 5L
      ang = textoidl('\AA')
   endelse

   ind_s = where(lam0 ge 3800. and lam0 le 9900.)
   if ind_s[0] ne -1 then yrange = [-1.0, max(median(flux[ind_s],20))*1.05]
   if not keyword_set(xrange) then xrange = [minwave, maxwave]
   len = 0.04
   pos = [0.08, 0.5, 0.97, 0.98]  ;else pos = [0.12, 0.14, 0.96, 0.98]
   ytitle = textoidl('Flux Density f_\lambda (10^{-17} erg s^{-1} cm^{-2} ') $
      + ang + textoidl('^{-1})')
   xtitle = textoidl('Rest Wavelength (')+ang+')'
   if not keyword_set(nsmooth) then nsmooth = 1L
   plot, wave, smooth(flux,nsmooth), xrange=xrange, yrange = yrange, /ystyle $
      , charsize=charsize, xthick=thick, ythick=thick, charthick=thick, _extra=extra $
      , /xsty,xticklen=len, yticklen = len/2, pos = pos, thick = linethick
   ind_bad = where(fitmask eq 0, nbad)
   if nbad gt 0 then oplot, [wave[ind_bad]], [flux[ind_bad]] $
               , psym=4, thick =thick, color=cgcolor('cyan'),symsize=0.5

   ;oplot, lam, flux, psym=2, color=cgcolor('black'), symsize=0.6
   oplot, wave, f_conti_model, color = cgcolor('red')
   oplot, wave, f_pl_model+f_poly_model, color=cgcolor('brown')
   oplot, wave, f_fe_mgii_model, color=cgcolor('blue')
   oplot, wave, f_fe_balmer_model, color=cgcolor('blue')
   oplot, wave, f_bc_model, color=cgcolor('dark gray'), line=1

   for jj=0L, (size(window_all))[2] - 1 do begin
     oplot,[window_all[0,jj],window_all[0,jj]],[0.9, 0.8]*yrange[1],color=cgcolor('dark gray'),thick=5
     oplot,[window_all[1,jj],window_all[1,jj]],[0.9, 0.8]*yrange[1],color=cgcolor('dark gray'),thick=5
     oplot,[window_all[0,jj],window_all[1,jj]],[0.9, 0.9]*yrange[1],color=cgcolor('dark gray'),thick=5
   endfor

   ; now plot the windows of each emission line fitting
   ;ind = where( ((wave ge 6400. and wave le 6800.)   or $  ; Halpha
   ;              (wave ge 4700. and wave le 5100.)   or $  ; Hbeta
   ;              (wave ge 2700. and wave le 2900.)   or $  ; MgII
   ;              (wave ge 1820. and wave le 1970.)   or $  ; CIII] complex
   ;              (wave ge 1500. and wave le 1600.))     $  ; CIV
   ;
   ; pos = [0.84, 0.1, 0.97, 0.44] & d_pos = [0.19, 0., 0.19, 0.]  ; this is for 5 panels
   ; pos = [0.8, 0.1, 0.97, 0.44] & d_pos = [0.24, 0., 0.24, 0.]  ; this is for 4 panels
   uniq_linecomp = uniq(linelist.compname)
   ncomp = n_elements(uniq_linecomp)
   compname = (linelist.compname)[uniq_linecomp]
   ; determine how many line complex are fitted and should be plotted
   plotkey=lonarr(ncomp)
   for ii=0L, ncomp - 1L do begin
      ind_line = where(linelist.compname eq compname[ii] and linelist.fitflag eq 1, $
                nline_fit)
      if nline_fit gt 0 then begin
         linelist_fit = linelist[ind_line]
         pxrange=[linelist_fit[0].minwav,linelist_fit[0].maxwav]
         ind = where(wave ge pxrange[0] and wave le pxrange[1], nnn)
         if nnn gt 10 then plotkey[ii] = 1L ; plot this line complex
      endif
   endfor
   temp=where(plotkey eq 1, nplot)
   if nplot gt 0 then begin
      ww=(0.89 - 0.06*(nplot - 1.))/nplot
      dp=ww + 0.06
      pos = [0.97 - ww, 0.1, 0.97, 0.44] & d_pos = [dp, 0., dp, 0.]
   endif 

   ; now plot each fitted line complex
   for ii=0L, ncomp - 1L do begin

       ind_line = where(linelist.compname eq compname[ii] and linelist.fitflag eq 1, $
                nline_fit)

       if nline_fit gt 0 then begin
          linelist_fit = linelist[ind_line]
          ; splog, 'Plotting line complex: ', compname[ii]
          iline_start = ind_line[0] & iline_end = ind_line[nline_fit-1]
          parindx_sub = parindx[0,iline_start] + indgen(3*linelist[iline_start].ngauss)
          for iitmp=1, nline_fit - 1L do begin
             parindx_sub = [parindx_sub, parindx[0,ind_line[iitmp]] $
              + indgen(3*linelist[ind_line[iitmp]].ngauss)]
          endfor
          pxrange=[linelist_fit[0].minwav,linelist_fit[0].maxwav]
          ind = where(wave ge pxrange[0] and wave le pxrange[1], nnn)
          ind_bad=where(wave ge pxrange[0] and wave le pxrange[1] and fitmask eq 0, nbad)
          ind_good=where(wave ge pxrange[0] and wave le pxrange[1] and fitmask ne 0, ngood)
          if nnn gt 10 then begin
             ;yrange=[-0.5,1.05]*max((smooth(line_flux[where(flux*sqrt(ivar) gt 5.)],20))[ind])
             if ngood gt 10 then yrange=[min(median(line_flux[ind_good],5)), max(median(line_flux[ind_good],5))*1.05] else $
                yrange=[min(median(line_flux[ind],5)), max(median(line_flux[ind],5))*1.05]
             plot, wave[ind], smooth(line_flux[ind],nsmooth), /noerase $
              , xrange=pxrange, yrange=yrange $
              , pos = pos,/xsty,thick=linethick, charsize=charsize,xticklen=len $
              , yticklen=len, xthick=thick, ythick=thick, charthick=thick $
              , xticks = xticks, xminor = xminor
             oplot, wave[ind],err[ind],color=cgcolor('gray'),thick=linethick
             if nbad gt 0 then oplot, [wave[ind_bad]], [line_flux[ind_bad]] $
               , psym=4, thick =thick, color=cgcolor('cyan'),symsize=0.5
             if keyword_set(fit_line) then begin
                pp = line_fit_all[ parindx_sub ]
                pp_broad = pp[0:3*linelist_fit[0].ngauss-1]
                if n_elements(pp) gt n_elements(pp_broad) then begin
                   pp_narrow = pp[3*linelist_fit[0].ngauss:*]
                   for jj=1L, n_elements(pp_narrow)/3L do begin
                      pp_narrow_1gauss = pp_narrow[(jj-1)*3L:jj*3L-1]
                      oplot, wave[ind], onegauss(alog(wave[ind]), pp_narrow_1gauss) $
                        , color=cgcolor('cyan')
                   endfor
                endif
                oplot,wave[ind], manygauss(alog(wave[ind]),pp_broad) $
                 , color=cgcolor('green')
                ; now plot total
                oplot,wave[ind],manygauss(alog(wave[ind]),pp) $
                  , color=cgcolor('red')

                xyouts, pos[2]-0.05,pos[3]-0.03, compname[ii], /norm,charsize=charsize
                ; plot floating linenames
                if keyword_set(more_anno) then begin
                  for jline=0,n_elements(linelist_fit)-1 do begin
                    xyouts, linelist_fit[jline].lambda, yrange[1]*(0.9-(jline mod 2)*0.1), $
                       linelist_fit[jline].linename,color=cgcolor('magenta'),charsize=0.7,align=0.5
                  endfor
                endif
             endif
             pos = pos - d_pos
          endif
      endif
      ; pos = pos - d_pos  ; this will not start the line panels from the right edge
   endfor

   xyouts,0.5,0.01,xtitle,/norm,charsize=charsize,charthick=thick,align=0.5
   xyouts,0.02,0.5,ytitle,/norm,charsize=charsize,charthick=thick,align=0.5,orien=90

   if keyword_set(objtag) then $
     xyouts, 0.15, 0.942, objtag, /norm, charsize = charsize, charthick=thick
   if keyword_set(SDSS_name) then xyouts, 0.5, 0.942, /norm, SDSS_name $
     , charsize = charsize, charthick = thick
   xyouts, 0.85, 0.942, 'z='+string(z,format='(f5.3)'), /norm $
        , charsize=charsize, charthick=thick, color=cgcolor('red')

   if keyword_set(psplot) and not keyword_set(append) then begin
      endplot
      cgfixps, figfile
      spawn, 'gzip -f '+ figfile
   endif

   ;--------------
   ; Compile line properties (peak wave, FWHM, line flux)
   if keyword_set(get_line_prop) then begin
      struct = replicate({linewave:0.D,linefwhm:0.D,linearea:0.D}, nline)
      linelist = struct_addtags(linelist, struct)
      for i=0L, nline - 1L do begin
         if linelist[i].fitflag eq 1 then begin
            ind = where(output.linename eq linelist[i].linename )
            pp = (output.line_fit)[ind]
            prop = get_multi_gaussian_prop(pp)
            linelist[i].linewave = exp(prop[0]) ; peak wavelength
            linelist[i].linefwhm = prop[1]*cs  ; FWHM
            linelist[i].linearea = prop[2]*(1. + z) ; line flux, scale back to true flux
         endif
      endfor
   endif  

end
