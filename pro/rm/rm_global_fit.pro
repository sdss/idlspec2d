;+
; NAME:
;   rm_global_fit
;
; PURPOSE:
;   Globally fit the QSO spectrum from Lya to Halpha
;
; CALLING SEQUENCE:
;   rm_global_fit, obs_wave, flux, err, z, [ra=,dec=,/psplot,/fits]
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; PROCEDURES CALLED
; 
;
; REVISION HISTORY:
;
;
;


function Fe_flux_mgii, xval, pp

   common Fe_temp, wave_Fe_civ, flux_fe_civ, wave_Fe_mgii, flux_Fe_mgii $
         , wave_Fe_balmer, flux_Fe_balmer

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
        else sig_conv = sqrt(Fe_FWHM^2 - 900.0^2)/2./sqrt(2.*alog(2.))  ; in units of km/s

      ; Get sigma in pixel space
      sig_pix = sig_conv/106.3  ; where 106.3 km/s is the dispersion for the Boroson Fe template

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
    yval = yval + pp[i]*xval2^i

  return, yval
end

function f_conti_only, xval, pp

   f_conti_BC = balmer_conti(xval, pp[2:4])  ; Balmer continuum
   f_pl = pp[0]*(xval/3000.0)^pp[1]  ; power-law continuum
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
   ; pp[11:15]: 4th order polynomial for the continuum

   common objdata, ra_c, dec_c, z_c

   nnn = n_elements(xval)
   yval = dblarr(nnn)
   f_Fe_MgII = Fe_flux_mgii(xval, pp[0:2])
   f_Fe_Balmer = Fe_flux_balmer(xval, pp[3:5])
   f_conti_BC = balmer_conti(xval, pp[8:10])  ; Balmer continuum

   f_pl = pp[6]*(xval/3000.0)^pp[7]  ; power-law continuum
   f_poly = f_poly_conti(xval, pp[11:*])
   yval = f_pl + f_Fe_MgII + f_Fe_Balmer + f_conti_BC + f_poly

   return, yval

end

;############################### Main Routine ########################################

pro rm_global_fit, lam0, flux0, err0, z, ra=ra,dec=dec, deredden=deredden $
       , conti_fit=conti_fit, fit_flag = fit_flag, line_fit_flag = line_fit_flag $
       , f_conti_model = f_conti_model, xrange = xrange, outdir = outdir $
       , psplot = psplot, fits= fits, output_name = output_name, silent = silent $
       , fit_line = fit_line, para=para, SDSS_name = SDSS_name $
       , wave_range = wave_range $ ; setup the restframe wavelength to fit
       , input_fitmask = input_fitmask $ ; input additional fitmasks, e.g., user-defined absorption masks
       , poly_ord=poly_ord, rej_iter=rej_iter $
       , tie_CIII = tie_CIII $ ; tie AlIII and SiIII redshift to that of CIII]
       , add_noise = add_noise, local_CIV_CIII_fit = local_CIV_CIII_fit $
       , sig_uv_nl = sig_uv_nl, voff_uv_nl = voff_uv_nl $ ; [sig_uv_nl,voff_uv_nl] are constraints for the narrow line component 
       , ngauss_broad_Ha = ngauss_broad_Ha, ngauss_broad_Hb = ngauss_broad_Hb $
       , ngauss_broad_MgII = ngauss_broad_MgII, ngauss_broad_CIII = ngauss_broad_CIII $
       , ngauss_broad_CIV = ngauss_broad_CIV, no_narrow_line = no_narrow_line $
       , auto_abs_rej = auto_abs_rej, conti_auto_abs_rej = conti_auto_abs_rej $
       , HeII_fit = HeII_fit, nsmooth=nsmooth, diet=diet $
       , sim_conti_fe_fit=sim_conti_fe_fit

   common objdata, ra_c, dec_c, z_c

   ; Define constants
   c = 2.9979246d5  ; speed of light, km/s
   if n_elements(z) eq 0 then z = 0. ; default is in restframe
   if n_elements(conti_auto_abs_rej) eq 0 then conti_auto_abs_rej = 1L
   if n_elements(auto_abs_rej) eq 0 then auto_abs_rej = 1L ; default is to auto reject absorption from the first linefit
   if n_elements(fit_line) eq 0 then fit_line = 1L ; default is to fit the lines

   lam = lam0 & flux = flux0 & err = err0 
   if keyword_set(add_noise) then begin ; add gaussian noise
      npix = n_elements(flux)
      flux = flux + randomn(seed, npix)*err
   endif
   ra_c = ra & dec_c = dec & z_c = z

   ; fit_flag = [fit_balmer_conti,not_fit_UV_Fe,fix_mgii_fe_FWHM,fix_balmer_fe_FWHM,fit_poly_conti]
   ; set the flag to fit the component
   ; In some cases the Balmer continuum is poorly constrained, 
   ; and hence will degrade the global power-law fit
   if n_elements(fit_flag) eq 0 then fit_flag = [0,0,0,0,1]
   if fit_flag[1] eq 1 then fit_flag[2] = 0L ; if not fitting UV FeII template, then do not set fix_mgII_fe_FWHM=1

   ; Default is *NOT* to fit a narrow component for MgII, CIII] and CIV
   if n_elements(no_narrow_line) eq 0 then no_narrow_line = 1L

   ; Default is to fit the HeII1640 complex
   if n_elements(HeII_fit) eq 0 then HeII_fit = 1L

   ; Default is to add a 4th order polynomial to the continuum
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

   ; mask out 3-sigma outlier of the 20-pix smoothed spectrum, to reduce the effects of
   ; absorption
   fitmask = ivar NE 0.0

   ; Mask outliers and absorption in the emission lines by finding 5 sigma
   ; outliers from smoothed spectrum
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
                             , upper = 20, lower = 5, nord = 3 $
                             , maxrej = 50, outmask = outmask, /silent)
      fitmask = fitmask*outmask
   endif
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

   window_all = [[window0], [window1], [window2], [window3_1], $
              [window3_2], [window4], [window5], [window6], $
              [window7], [window8], [window9], [window10], [window11] ]

   window_conti = [[window0], [window1],[window2],[window5],[window8],[window11]]
   window_FeII = [[window3_1],[window3_2],[window4],[window6],[window7],[window9],[window10]]

   ; Shift to restframe
   wave = lam/(1.0+z)

   ; Cut the wavelength in the windows specified by keywords wave_range
   if keyword_Set(wave_range) then begin
      ind = where(wave ge wave_range[0] and wave le wave_range[1])
      wave = wave[ind] & flux = flux[ind] & err = err[ind]
      fitmask = fitmask[ind]
   endif

   ; This is the pixels for simultaneous continuum+FeII fit
   ind = where( ( (wave ge window0[0] and wave le window0[1]) or $
                  (wave ge window1[0] and wave le window1[1]) or $
                  (wave ge window2[0] and wave le window2[1]) or $
                  (wave ge window3_1[0] and wave le window3_1[1]) or $
                  (wave ge window3_2[0] and wave le window3_2[1]) or $
                  (wave ge window4[0] and wave le window4[1]) or $
                  (wave ge window5[0] and wave le window5[1]) or $
                  (wave ge window6[0] and wave le window6[1]) or $
                  (wave ge window7[0] and wave le window7[1]) or $
                  (wave ge window8[0] and wave le window8[1]) or $
                  (wave ge window9[0] and wave le window9[1]) or $
                  (wave ge window10[0] and wave le window10[1]) or $
                  (wave ge window11[0] and wave le window11[1]) ) $
                 AND fitmask ne 0. )
   ind_conti_only = where( ( (wave ge window0[0] and wave le window0[1]) or $
                             (wave ge window1[0] and wave le window1[1]) or $
                             (wave ge window2[0] and wave le window2[1]) or $
                             (wave ge window5[0] and wave le window5[1]) or $
                             (wave ge window8[0] and wave le window8[1]) or $
                             (wave ge window11[0] and wave le window11[1])) $
                           AND fitmask ne 0.)
   ind_Fe_only = where( ( (wave ge window3_1[0] and wave le window3_1[1]) or $
                  (wave ge window3_2[0] and wave le window3_2[1]) or $
                  (wave ge window4[0] and wave le window4[1]) or $
                  (wave ge window6[0] and wave le window6[1]) or $
                  (wave ge window7[0] and wave le window7[1]) or $
                  (wave ge window9[0] and wave le window9[1]) or $
                  (wave ge window10[0] and wave le window10[1]) ) $
                 AND fitmask ne 0. )

   if n_elements(ind) lt 10 then begin
      splog, 'Fitting pixel < 10. return. '
      return
   endif

   ; ## now fit the psedu-continuum
   ; set limiting conditions on paramters
   parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0],tied:''}, 11+poly_ord+1)
   parinfo.value = [1.0D, 3000.0D, 0.0, 1.0, 3000.0D, 0.0, 1.0D, -2.0D, 1.0, 15000.D, 0.5 $
        , replicate(0.,poly_ord+1) ]
   ; limit the PL component
   parinfo[6].limited = [1,0] & parinfo[6].limits = [0., 1d10]
   parinfo[7].LIMITED = [1, 1] & parinfo[7].LIMITS = [-5.0,3.0]
   ; limit the Fe_MgII component
   parinfo[0].limited = [1, 0] & parinfo[0].limits = [0., 1d10]
   parinfo[1].LIMITED = [1, 1] & parinfo[1].LIMITS = [900, 1d4]
   parinfo[2].LIMITED = [1, 1] & parinfo[2].LIMITS = [-0.01, 0.01]
   if fit_flag[1] eq 1 then begin ; do not fit the UV iron template
      parinfo[0].limited = [0,0] & parinfo[0].value = 0. 
      parinfo[0:2].fixed = 1L
   endif
   if fit_flag[2] eq 1 then begin
      parinfo[1].fixed = 1L ; fix the FWHM of the MgII iron
   endif
   ; limit the Fe_balmer component
   parinfo[3].limited = [1, 0] & parinfo[3].limits = [0., 1d10]
   parinfo[4].LIMITED = [1, 1] & parinfo[4].LIMITS = [900, 1d4]
   if fit_flag[3] eq 1 then begin
      parinfo[4].fixed = 1L ; fix the FWHM of the Balmer iron
   endif
   ; limit the Balmer continuum
   if fit_flag[0] eq 0 then begin ; if fit_flag not set then do not fit for Balmer continuum
      parinfo[8:10].fixed = 1L & parinfo[8].value = 0.
   endif else begin
      parinfo[8].limited = [1, 0] & parinfo[8].limits = [0., 1d10]
      parinfo[9].LIMITED = [1, 1] & parinfo[9].LIMITS = [10000., 50000.]
      parinfo[10].LIMITED = [1, 1] & parinfo[10].LIMITS = [0.1, 2.]
   endelse
   ; limit the polynomial continuum component
   if fit_flag[4] eq 0 then begin
      parinfo[11:*].value = 0. & parinfo[11:*].fixed = 1L
   endif
   ; enforce the optical FeII to be zero if not covered
   ind_cover = where( ( (wave ge window6[0] and wave le window6[1]) or $
                        (wave ge window7[0] and wave le window7[1]) ) $
                     and fitmask ne 0)
   if n_elements(ind_cover) lt 100 then begin
      parinfo[3].value = 0. & parinfo[3].fixed = 1L
   endif


   ; First try to fit the continuum-alone
   if n_elements(ind_conti_only) ge 50 and not keyword_set(sim_conti_fe_fit) then begin
      
      parinfo1 = parinfo[6:*]
      conti_only_fit = mpfitfun('f_conti_only', wave[ind_conti_only], flux[ind_conti_only] $
                     , err[ind_conti_only] $
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
                        , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
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
      ;   parinfo[8].fixed = 1L & parinfo[8].value = 0. 
      ;   parinfo[9].value = 15000.D & parinfo[10].value=0.5
      ;   parinfo[9].fixed = 1L & parinfo[10].fixed = 1L
      ;   parinfo[8].limited = [0,0] & parinfo[8].limits = [0.D,0]
      ;   parinfo[9].limited = [0,0] & parinfo[9].limits = [0.d,0]
      ;   parinfo[10].limited = [0,0] & parinfo[10].limits = [0.D,0]
      ;
      ;   conti_fit2 = mpfitfun('f_conti_all', wave[ind], flux[ind], err[ind] $
      ;            , parinfo = parinfo, perror = perror2, yfit = yfit2, /quiet $
      ;            , nfev = nfev, niter = niter, status = status2, bestnorm = chi2)
      ;   dof = n_elements(ind) - (n_elements(parinfo) - 3L - 3L*(fit_flag[1]) - fit_flag[2] - fit_flag[3] - 1L -  (1L-fit_flag[4]) - (1-fit_flag[4])*5L )
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
      ; this is to avoid situations where one of the continuum windows falls within a BAL trough
      if keyword_Set(conti_auto_abs_rej) and redchi2 lt 100. then begin
         fitmask2 =  replicate(1L, n_elements(wave))
         ind_abs = where(flux[ind] lt yfit1 - 3.*err[ind])
         if ind_abs[0] ne -1 then fitmask2[ind[ind_abs]] = 0L
         fitmask2 = fitmask*fitmask2
         ind2 = where( ( (wave ge window0[0] and wave le window0[1]) or $
                     (wave ge window1[0] and wave le window1[1]) or $
                     (wave ge window2[0] and wave le window2[1]) or $
                     (wave ge window3_1[0] and wave le window3_1[1]) or $
                     (wave ge window3_2[0] and wave le window3_2[1]) or $
                     (wave ge window4[0] and wave le window4[1]) or $
                     (wave ge window5[0] and wave le window5[1]) or $
                     (wave ge window6[0] and wave le window6[1]) or $
                     (wave ge window7[0] and wave le window7[1]) or $
                     (wave ge window8[0] and wave le window8[1]) or $
                     (wave ge window9[0] and wave le window9[1]) or $
                     (wave ge window10[0] and wave le window10[1]) or $
                     (wave ge window11[0] and wave le window11[1]) ) $
                     AND fitmask2 ne 0. )
         if (n_elements(ind2) gt 10L) and (n_elements(ind2) lt n_elements(ind)) then begin
            conti_fit2 = mpfitfun('f_conti_all', wave[ind2], flux[ind2], err[ind2] $
                     , parinfo = parinfo, perror = perror2, yfit = yfit2, /quiet $
                     , nfev = nfev, niter = niter, status = status2, bestnorm = chi2)
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
   output = create_struct('wave', wave, 'flux', flux, 'err', err, 'fitmask', fitmask, $
         'ind_conti_fit', ind, 'conti_fit', conti_fit, 'conti_fit_err', conti_err, $
         'conti_redchi2', redchi2, 'conti_status', status)
   output = struct_addtags(specdata, output)

   f_pl_model = conti_fit[6]*(wave/3000.0)^conti_fit[7]
   f_fe_mgii_model = fe_flux_mgii(wave, conti_fit[0:2])
   f_fe_balmer_model = fe_flux_balmer(wave, conti_fit[3:5])
   f_bc_model = balmer_conti(wave, conti_fit[8:10])
   f_poly_model = f_poly_conti(wave, conti_fit[11:*]) 
   f_conti_model = f_pl_model + f_fe_mgii_model + f_fe_balmer_model + f_bc_model + f_poly_model

   if keyword_set(local_CIV_CIII_fit) then begin ; if a local fit for CIV+CIII] is required [1000, 2165]
      ind = where( ( (wave ge window0[0] and wave le window0[1]) or $
               (wave ge window1[0] and wave le window1[1]) or $
               (wave ge window2[0] and wave le window2[1]) or $
               (wave ge window3_1[0] and wave le window3_1[1]) or $
               (wave ge window3_2[0] and wave le window3_2[1]) or $
               (wave ge window4[0] and wave le window4[1]) or $
               (wave ge window5[0] and wave le window5[1]) or $
               (wave ge window6[0] and wave le window6[1]) or $
               (wave ge window7[0] and wave le window7[1]) or $
               (wave ge window8[0] and wave le window8[1]) or $
               (wave ge window9[0] and wave le window9[1])) $
              AND fitmask ne 0.  AND wave le 2165.)
      ; first turn off Balmer continuum and UV iron template fit
      parinfo[8].fixed = 1L & parinfo[8].value = 0.
      parinfo[9].fixed = 1L & parinfo[10].fixed = 1L  
      parinfo[2].limited = [0,0] & parinfo[2].value = 0. & parinfo[2].fixed = 1L
      parinfo[3].fixed = 1L & parinfo[4].fixed = 1L
      conti_fit3 = mpfitfun('f_conti_all', wave[ind], flux[ind], err[ind] $
                 , parinfo = parinfo, perror = perror3, yfit = yfit3, /quiet $
                 , nfev = nfev, niter = niter, status = status3, bestnorm = chi2)
      if not keyword_set(silent) then splog, 'CONTI-FITTING: MPFIT nfev=', nfev, $
        ' niter=', niter, ' status=', status3
      ; update the psedu-continuum model, connected at 2165A
      ind_local = where(wave le 2165.)
      f_pl_model[ind_local] = conti_fit3[0]*(wave[ind_local]/3000.0)^conti_fit3[1]
      f_fe_mgii_model[ind_local] = 0.
      f_fe_balmer_model[ind_local] = 0.
      f_bc_model[ind_local] = 0.
      f_conti_model = f_pl_model + f_fe_mgii_model + f_fe_balmer_model + f_bc_model + f_poly_model

      temp = create_struct('conti_fit3', conti_fit3, 'conti_fit3_err', perror3)
      output = struct_addtags(output, temp)
   endif

   ; Start to fit the emission lines: CIV, CIII], MgII, Hbeta, Halpha
   ; first subtract the psedu-continuum to get emission line fluxs

   line_flux = flux - f_conti_model
   if keyword_set(local_CIV_CIII_fit) then $
     line_flux[ind_local] = flux[ind_local] - f_conti_all(wave[ind_local], conti_fit3)
   output = struct_addtags(output, {line_flux:line_flux})


   ; Setup line fitting windows
   win_civ = [1500., 1600.]
   if keyword_set(HeII_fit) then win_civ = [1500., 1700.]
   ind = where( ((wave ge 6400. and wave le 6800.)   or $  ; Halpha
              (wave ge 4700. and wave le 5100.)   or $  ; Hbeta
              (wave ge 2700. and wave le 2900.)   or $  ; MgII
              (wave ge 1820. and wave le 1970.)   or $  ; CIII] complex
              (wave ge win_civ[0] and wave le win_civ[1]))     $  ; CIV
              AND fitmask ne 0.   ) 

   if n_elements(ind) lt 10L then begin
      splog, 'Line fitting pixels<10. Return'
      return
   endif

   ; setup fitting parameters for emission line fitting
   if not keyword_set(ngauss_broad_Ha) then ngauss_broad_Ha = 3L
   if not keyword_set(ngauss_broad_Hb) then ngauss_broad_Hb = 3L
   if not keyword_set(ngauss_broad_MgII) then ngauss_broad_MgII = 3L
   if not keyword_set(ngauss_broad_CIII) then ngauss_broad_CIII = 2L
   if not keyword_set(ngauss_broad_CIV) then ngauss_broad_CIV = 3L
   ngauss = ngauss_broad_Ha + 5L + ngauss_broad_Hb + 3L $
       + ngauss_broad_MgII + 1L + ngauss_broad_CIII + 3L $
       + ngauss_broad_CIV + 1L

   ; adding two more Gaussians for HeII 1640 and OIII] 1663
   if keyword_set(HeII_fit) then ngauss = ngauss + 2L 

   parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0],tied:''}, ngauss*3L)

   ; default maxmium dispersion for the broad Gaussians
   if not keyword_set(max_sigma) then max_sigma = 0.05D

   ; only fit the line if the keyword FIT_LINE is set
   if keyword_set(fit_line) then begin

      ; For Halpha; 2 gaussians for [SII], 2 gaussians for [NII],
      ; 1 gaussian for narrow Ha and ngauss_broad_Ha gaussians for broad Ha
      for i_gau=0L, ngauss_broad_Ha - 1L do begin
          parinfo[0+i_gau*3L:2+i_gau*3L].value = [1.0D/ngauss_broad_Ha, alog(6564.61D), 0.005D]
          parinfo[0+i_gau*3L].limited = [1L, 0L]
          parinfo[0+i_gau*3L].limits = [0.0, 0.0]
          parinfo[1+i_gau*3L].limited = [1L,1L]
          parinfo[1+i_gau*3L].limits = [alog(6564.61) - 0.015D, alog(6564.61) + 0.015D]
          parinfo[2+i_gau*3L].limited = [1L, 1L]
          ; restrict each component to have 1400<FWHM<2.3548*c*max_sigma km/s
          parinfo[2+i_gau*3L].limits = [0.002D, max_sigma]

          ; restrict the centroid of each broad component is the same
          ; hence the overall broad line is symmetric
          ; if i_gau gt 0L then parinfo[1+i_gau*3L].tied = 'P(1)'
      endfor
      inext = (ngauss_broad_Ha - 1L)*3
      ; the Hbeta narrow component (FWHM < 1200 km/s)
      parinfo[3+inext:5+inext].value = [0.01D, alog(6564.61D), 0.001D]
      parinfo[3+inext].limited = [1L, 0L]
      parinfo[3+inext].limits = [0.0, 0.0]
      parinfo[4+inext].limited = [1L,1L]
      parinfo[4+inext].limits = [alog(6564.61D) - 5d-3, alog(6564.61D) + 5d-3]
      parinfo[5+inext].limited = [1L,1L]
       ; restrict the narrow component to have FWHM<1200 km/s but > 160km/s
      parinfo[5+inext].limits = [2.3d-4, 0.0017D]
      ; NII6549
      parinfo[6+inext:8+inext].value = [0.01D, alog(6549.85D), 0.001D]
      parinfo[6+inext].limited = [1L, 0L]
      parinfo[6+inext].limits = [0., 1d10]
      pstr1 = '(P(' + string(4+inext, format='(i0)') + ')'
      parinfo[7+inext].tied = pstr1 +  ' - alog(6564.61D)) + alog(6549.85D)'
      pstr2 = 'P(' + string(5+inext, format='(i0)') + ')'
      parinfo[8+inext].tied = pstr2
      ; NII6585
      parinfo[9+inext:11+inext].value = [0.01D, alog(6585.28D), 0.001D]
      parinfo[9+inext].limited = [1L, 0L]
      parinfo[9+inext].limits = [0., 1d10]
      parinfo[10+inext].tied = pstr1 + ' - alog(6564.61D)) + alog(6585.28D)'
      parinfo[11+inext].tied = pstr2
      pstr3 = 'P(' + string(6+inext, format='(i0)') + ')' + '*6549.85/6585.28*3'
      parinfo[9+inext].tied = pstr3
      ; SII 6718
      parinfo[12+inext:14+inext].value = [0.001D, alog(6718.29D), 0.001D]
      parinfo[12+inext].limited = [1L, 0L]
      parinfo[12+inext].limits = [0., 1d10]
      parinfo[13+inext].tied = pstr1 + ' - alog(6564.61D)) + alog(6718.29D)'
      parinfo[14+inext].tied = pstr2
      ; SII 6732
      parinfo[15+inext:17+inext].value = [0.001D, alog(6732.67D), 0.001D]
      parinfo[15+inext].limited = [1L, 0L]
      parinfo[15+inext].limits = [0., 1d10]
      parinfo[16+inext].tied = pstr1 + ' - alog(6564.61D)) + alog(6732.67D)'
      parinfo[17+inext].tied = pstr2
      ; Check if Halpha is covered. If not, then fix all Halpha parameters
      ind_cover = where(wave[ind] ge 6400. and wave[ind] le 6800., n_cover)
      if n_cover lt 100 then begin
         istart = 0 & iend = 17+inext
         ind1 = istart + 3*indgen((iend-istart+1)/3L)
         parinfo[ind1].value = 0. & parinfo[istart:iend].fixed = 1
      endif

      ;-------------------------End of Halpha setup ------
      ; for Hbeta
      inext = 17+inext + 1L ; advanced to Hb
      istart = inext
      for i_gau=0L, ngauss_broad_Hb - 1L do begin
          parinfo[inext+i_gau*3L:inext+2+i_gau*3L].value= $
            ; [1.0D/ngauss_broad_Hb, alog(4862.68D), 0.005D] 
            [0.1D/ngauss_broad_Hb, alog(4862.68D)+ 6.67d-4*(i_gau + 1 - ngauss_broad_Hb),0.005D]
          parinfo[inext+i_gau*3L].limited = [1L, 0L]
          parinfo[inext+i_gau*3L].limits = [0.0, 0.0]
          parinfo[inext+1+i_gau*3L].limited = [1L,1L]
          parinfo[inext+1+i_gau*3L].limits = [alog(4862.68D) - 0.015D, alog(4862.68D) + 0.015D]
          parinfo[inext+2+i_gau*3L].limited = [1L, 1L]
          ; restrict each component to have 1400<FWHM<2.3548*c*max_sigma km/s
          parinfo[inext+2+i_gau*3L].limits = [0.002D, max_sigma]
          ; restrict the centroid of each broad component is the same, 
          ; hence the overall broad line is symmetric
          ; if i_gau gt 0L then parinfo[1+i_gau*3L].tied = 'P(1)'
      endfor
      inext = inext + (ngauss_broad_Hb-1)*3
      ; the Hbeta narrow component (FWHM < 1200 km/s)
      parinfo[3+inext:5+inext].value = [0.01D, alog(4862.68D), 0.001D]
      parinfo[3+inext].limited = [1L, 0L]
      parinfo[3+inext].limits = [0.0, 0.0]
      parinfo[4+inext].limited = [1L,1L]
      parinfo[4+inext].limits = [alog(4862.68) - 0.015D, alog(4862.68) + 0.015D]
      parinfo[4+inext].tied = pstr1 + ' - alog(6564.61D)) + alog(4862.68D)'
      parinfo[5+inext].limited = [1L,1L]
      ; restrict the narrow component to have FWHM<1200 km/s
      parinfo[5+inext].limits = [2.3d-4, 0.0017D]
      parinfo[5+inext].tied = pstr2
      ; OIII 4959
      parinfo[6+inext:8+inext].value = [0.01D, alog(4960.30D), 0.001D]
      parinfo[6+inext].limited = [1L, 0L]
      parinfo[6+inext].limits = [0., 1d10]
      parinfo[7+inext].tied = pstr1 +  ' - alog(6564.61D)) + alog(4960.30D)'
      parinfo[8+inext].tied = pstr2
      ; the OIII 5007
      parinfo[9+inext:11+inext].value = [0.02D, alog(5008.24D), 0.001D]
      parinfo[9+inext].limited = [1L, 0L]
      parinfo[9+inext].limits = [0., 1d10]
      pstr3_OIII = 'P(' + string(6+inext, format='(i0)') + ')' + '*4960.30/5008.24*3.'
      parinfo[9+inext].tied = pstr3_OIII
      parinfo[10+inext].tied = pstr1 + ' - alog(6564.61D)) + alog(5008.24D)'
      parinfo[11+inext].tied = pstr2
      ; Check if Hbeta is covered. If not, then fix all Hbeta parameters
      ind_cover = where(wave[ind] ge 4700. and wave[ind] le 5100., n_cover)
      if n_cover lt 100 then begin
         iend = 11+inext
         ind1 = istart + 3*indgen((iend-istart+1)/3L)
         parinfo[ind1].value = 0. & parinfo[istart:iend].fixed = 1
      endif
      ;----------------------------Done with hbeta---------------------------
      ; for MgII
      inext = inext + 12L
      istart = inext
      for i_gau=0L, ngauss_broad_MgII - 1L do begin
         parinfo[inext+i_gau*3L:inext+2+i_gau*3L].value = $
            ; [1.0D/ngauss_broad_MgII, alog(2798.75D), 0.005D]
           [0.1D/ngauss_broad_MgII, alog(2798.75D)+ 6.67d-4*(i_gau + 1 - ngauss_broad_MgII), 0.005D]
         parinfo[inext+i_gau*3L].limited = [1L, 0L]
         parinfo[inext+i_gau*3L].limits = [0.0, 0.0]
         parinfo[inext+1+i_gau*3L].limited = [1L,1L]
         parinfo[inext+1+i_gau*3L].limits = [alog(2798.75D) - 0.015D, alog(2798.75D) + 0.015D]
         parinfo[inext+2+i_gau*3L].limited = [1L, 1L]
         ; restrict each component to have 1400<FWHM<2.3548*c*max_sigma km/s
         parinfo[inext+2+i_gau*3L].limits = [0.002D, max_sigma]
      endfor
      ; for the narrow MgII
      inext = inext + ngauss_broad_MgII*3L
      parinfo[inext:inext+2].value = [0.02D, alog(2798.75D), 0.001D]
      parinfo[inext].limited = [1L, 0L]
      parinfo[inext].limits = [0.0, 0.0] ; emission only
      parinfo[inext+1].tied = pstr1 + ' - alog(6564.61D)) + alog(2798.75D)'
      parinfo[inext+2].tied = pstr2
      if keyword_set(no_narrow_line) then begin
         parinfo[inext:inext+2].fixed = 1L & parinfo[inext].value = 0.D
      endif
      ; Check if MgII is covered. If not, then fix all MgII parameters
      ind_cover = where(wave[ind] ge 2700. and wave[ind] le 2900., n_cover)
      if n_cover lt 100 then begin
         iend = inext+2
         ind1 = istart + 3*indgen((iend-istart+1)/3L)
         parinfo[ind1].value = 0. & parinfo[istart:iend].fixed = 1
      endif
 
      ; --------------------------------------------------------------------
      ; for CIII]; becuase this is a complex of AlIII1857, SiIII]1892 and CIII] 1908, 
      ; only fit for broad AlIII1857, SiIII]1892 with a single Gaussian each, broad
      ; CIII] with 2 Gaussian with fixed centroids
      ; and another Gaussian for the narrow CIII], again with velocity and width tied
      ; to other NEls.
      inext = inext+3L
      istart = inext
      pstr_ciii = 'P(' + string(1+inext, format='(i0)')+')'
      pstr_ciii2 = '(P(' + string(1+inext, format='(i0)')+')'
      for i_gau=0L, ngauss_broad_CIII - 1L do begin
          ; [1.0D/ngauss_broad_CIII, alog(1908.73D), 0.005D]
          parinfo[inext+i_gau*3L:inext+2+i_gau*3L].value = $
          [0.01D/ngauss_broad_CIII, alog(1908.73D), 0.005D]
          parinfo[inext+i_gau*3L].limited = [1L, 0L]
          parinfo[inext+i_gau*3L].limits = [0.0, 0.0]
          parinfo[inext+1+i_gau*3L].limited = [1L,1L]
          parinfo[inext+1+i_gau*3L].limits = [alog(1908.73D) - 0.015D, alog(1908.73D) + 0.015D]
          if i_gau gt 0 then parinfo[inext+1+i_gau*3L].tied = pstr_ciii
          parinfo[inext+2+i_gau*3L].limited = [1L, 1L]
          ; restrict each component to have 1400<FWHM<2.3548*c*max_sigma km/s
          parinfo[inext+2+i_gau*3L].limits = [0.002D, max_sigma]
      endfor
      ; narrow CIII
      inext = inext + ngauss_broad_CIII*3L
      ; [0.02D, alog(1908.73D), 0.001D]
      parinfo[inext:inext+2].value = [0.002D, alog(1908.73D), 0.001D]
      parinfo[inext].limited = [1L, 0L]
      parinfo[inext].limits = [0.0, 0.0] ; emission only
      parinfo[inext+1].tied = pstr1 + ' - alog(6564.61D)) + alog(1908.73D)'
      parinfo[inext+2].tied = pstr2
      if keyword_set(no_narrow_line) then begin
         parinfo[inext:inext+2].fixed = 1L & parinfo[inext].value = 0.D
      endif
      ; two other broad Gaussians for ALIII1857 and SiIII]1892; 
      ; tie tie their velocity together but not to CIII]
      inext = inext + 3L
      ; [0.05D, alog(1892.03D), 0.005D], for SiIII]1892
      parinfo[inext:inext+2].value = [0.005D, alog(1892.03D), 0.005D]
      pstr_SiIII = '(P(' + string(1+inext, format='(i0)') + ')'
      pstr_SiIII2 = 'P(' + string(2+inext, format='(i0)') + ')'
      parinfo[inext].limited = [1L, 0L]
      parinfo[inext].limits = [0.0, 0.0] ; emission only
      parinfo[inext+1].limited = [1L,1L]
      parinfo[inext+1].limits = [alog(1892.03D) - 3.33d-3, alog(1892.03D) + 3.33d-3]
      if keyword_set(tie_CIII) then parinfo[inext+1].tied = $
           pstr_ciii2 + ' - alog(1908.73D)) + alog(1892.03D)'
      parinfo[inext+2].limited = [1L, 1L]
      ; for SiIII and AlIII, using slightly smaller FWHM upper limit
      parinfo[inext+2].limits = [0.002D, 0.02D]
      inext = inext + 3L
      ; [0.2D, alog(1857.40D), 0.005D], for AlIII]1857
      parinfo[inext:inext+2].value = [0.005D, alog(1857.40D), 0.005D]
      parinfo[inext].limited = [1L, 0L]
      parinfo[inext].limits = [0.0, 0.0] ; emission only
      if not keyword_set(tie_CIII) then parinfo[inext+1L].tied = pstr_SiIII $
           + ' - alog(1892.03D)) + alog(1857.40D)' $
      else parinfo[inext+1L].tied = pstr_ciii2 + ' - alog(1908.73D)) + alog(1857.40D)'
      parinfo[inext+2L].tied = pstr_SiIII2
      ;parinfo[inext+2].limited = [1L, 1L]
      ;parinfo[inext+2].limits = [0.0017D, max_sigma]
      ; Check if CIII is covered. If not, then fix all CIII parameters
      ind_cover = where(wave[ind] ge 1820. and wave[ind] le 1970., n_cover)
      if n_cover lt 100 then begin
         iend = inext+2
         ind1 = istart + 3*indgen((iend-istart+1)/3L)
         parinfo[ind1].value = 0. & parinfo[istart:iend].fixed = 1
      endif

      ;------------------------------
      ; now for CIV
      inext = inext+3L
      istart = inext
      for i_gau=0L, ngauss_broad_CIV - 1L do begin
         ; [1.0D/ngauss_broad_CIV, alog(1549.06D), 0.005D]
         parinfo[inext+i_gau*3L:inext+2+i_gau*3L].value = $
          [0.05D/ngauss_broad_CIV, alog(1549.06D) + 6.67d-4*(i_gau + 1 - ngauss_broad_CIV), $
           0.005D - 0.002*(i_gau + 1 - ngauss_broad_CIV)]
         parinfo[inext+i_gau*3L].limited = [1L, 0L]
         parinfo[inext+i_gau*3L].limits = [0.0, 0.0]
         parinfo[inext+1+i_gau*3L].limited = [1L,1L]
         parinfo[inext+1+i_gau*3L].limits = [alog(1549.06D) - 0.015D, alog(1549.06D) + 0.015D]
         parinfo[inext+2+i_gau*3L].limited = [1L, 1L]
         ; restrict each component to have 1400<FWHM<2.3548*c*max_sigma km/s
         parinfo[inext+2+i_gau*3L].limits = [0.002D, max_sigma]
      endfor
      inext = inext + ngauss_broad_CIV*3L ; narrow CIV, if any
      ;[0.05D, alog(1549.06D), 0.001D]
      parinfo[inext:inext+2].value = [0.01D, alog(1549.06D), 0.001D]
      parinfo[inext].limited = [1L, 0L]
      parinfo[inext].limits = [0.0, 0.0] ; emission only
      parinfo[inext+1].tied = pstr1 + ' - alog(6564.61D)) + alog(1549.06D)'
      parinfo[inext+2].tied = pstr2
      if keyword_set(no_narrow_line) then begin
         parinfo[inext:inext+2].fixed = 1L & parinfo[inext].value = 0.D
      endif
      ;--------------------------------------------------------------
      ; now for HeII 1640 and OIII] 1663
      if keyword_set(HeII_fit) then begin
         ; for HeII 1640
         inext = inext + 3L
         parinfo[inext:inext+2].value = [0.005D, alog(1640.42D), 0.005D]
         parinfo[inext].limited = [1L, 0L] & parinfo[inext].limits = [0., 0.] ; emission only
         parinfo[inext+1].limited = [1L, 1L]
         parinfo[inext+1].limits = [alog(1640.42D) - 0.008D, alog(1640.42D) + 0.008D]
         parinfo[inext+2].limited = [1L, 1L] & parinfo[inext+2].limits = [0.001D, max_sigma]
  
         ; for OIII] 1663
         HeII_cen_str = '(P(' + string(inext+1, format='(i0)') + ')' 
         HeII_sig_str = 'P(' + string(inext+2, format='(i0)') + ')'
         inext = inext + 3L
         parinfo[inext:inext+2].value = [0.005D, alog(1663.48D), 0.005D]
         parinfo[inext].limited = [1L, 0L] & parinfo[inext].limits = [0., 0.] ; emission only
         parinfo[inext+1].tied = HeII_cen_str + ' - alog(1640.42D)) + alog(1663.48D)'
         parinfo[inext+2].tied = HeII_sig_str
      endif
      ; Check if CIV is covered. If not, then fix all CIV parameters
      ind_cover = where(wave[ind] ge win_civ[0] and wave[ind] le win_civ[1], n_cover)
      if n_cover lt 100 then begin
         iend = inext+2
         ind1 = istart + 3*indgen((iend-istart+1)/3L)
         parinfo[ind1].value = 0. & parinfo[istart:iend].fixed = 1
      endif

      ; select the parameters that are being fitted
      indd_fit = where(parinfo.fixed eq 0)
      parinfo1 = parinfo[indd_fit]
      line_fit_all = replicate(0., n_elements(parinfo))
      perror_all = replicate(0., n_elements(parinfo))

      ; now proceed to fit the lines
      line_fit = mpfitfun('manygauss', alog(wave[ind]), line_flux[ind], err[ind] $
                     , MAXITER = 500L, parinfo = parinfo1 $
                     , perror = perror, yfit = yfit2, /quiet $
                     , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
      ttt = where(parinfo1.fixed eq 1 or strlen(parinfo1.tied) gt 0, nfix) 
      dof = n_elements(ind) - (n_elements(parinfo1) - nfix)
      redchi2 = chi2/dof
      if not keyword_set(silent) then splog, 'LINE-FITTING: MPFIT nfev=', nfev $
         , ' niter=', niter, ' status=', status, ' redchi2=', redchi2

      ; do a further step to remove 3sigma below the first fit
      ; first setup fitmask2 and the new line fitting windows
      if keyword_set(auto_abs_rej) and redchi2 lt 100. then begin
         splog, 'Start iterative absorption rejection: '
         n_iter = 1L

   iter1:
         fitmask2 =  replicate(1L, n_elements(wave))
         ind_abs = where(line_flux[ind] lt yfit2 - 3.*err[ind])
         if ind_abs[0] ne -1 then fitmask2[ind[ind_abs]] = 0L
         fitmask2 = fitmask*fitmask2

         ind2 = where( ((wave ge 6400. and wave le 6800.)   or $  ; Halpha
                (wave ge 4700. and wave le 5100.)   or $  ; Hbeta
                (wave ge 2700. and wave le 2900.)   or $  ; MgII
                (wave ge 1820. and wave le 1970.)   or $  ; CIII] complex
                (wave ge win_civ[0] and wave le win_civ[1]))     $  ; CIV
                AND fitmask2 ne 0.   )

         if (n_elements(ind2) gt 10L) and (n_elements(ind2) lt n_elements(ind)) then begin

            line_fit2 = mpfitfun('manygauss', alog(wave[ind2]), line_flux[ind2], err[ind2] $
                 , MAXITER = 500L, parinfo = parinfo1 $
                 , perror = perror2, yfit = yfit2_2, /quiet $
                 , nfev = nfev, niter = niter, status = status2, bestnorm = chi2_2)
            ttt = where(parinfo1.fixed eq 1 or strlen(parinfo1.tied) gt 0, nfix)
            dof2 = n_elements(ind2) - (n_elements(parinfo1) - nfix)
            redchi2_2 = chi2_2/dof2
            if not keyword_set(silent) then splog, 'LINE-FITTING: MPFIT nfev=', nfev $
               , ' niter=', niter, ' status=', status2, ' redchi2=', redchi2_2, ' iter=', n_iter
            ; replace the original fit if justified
            if status2 gt 0 and redchi2_2 gt 0 and redchi2_2 lt redchi2 then begin
               splog, 'Improved fits by rejecting absorptions. '
               ind = ind2 & line_fit = line_fit2 & perror = perror2 
               redchi2 = redchi2_2 & status = status2 & yfit2 = yfit2_2
               fitmask = fitmask2
               output.fitmask = fitmask
               if n_iter lt rej_iter then begin
                  n_iter = n_iter + 1L
                  goto, iter1
               endif
            endif

         endif
      endif

      line_fit_all[indd_fit] = line_fit
      perror_all[indd_fit] = perror
      temp = create_struct('ind_line_fit', ind, 'line_fit', line_fit_all $
        , 'line_fit_err', perror_all, 'line_redchi2', redchi2, 'line_status', status)
      output = struct_addtags(output, temp)

      ;splog, 'Output structure first created'

   endif
   ; --------------------------- End of fitting emission lines ---------------

   ; Write the fit results
   if n_elements(output_name) eq 0 then output_name = 'output'
   if keyword_Set(fits) then begin ; output the fits result
      fitsfile = outdir + '/fits/' + output_name + '.fits'
      mwrfits, output, fitsfile, /create
   endif

   ; now make a plot
   if keyword_set(psplot) then begin
      figfile = outdir + '/QA/' + output_name + '.eps'
      begplot, name=figfile, xsize = 40, ysize = 25
      charsize = 2. & thick = 4. & xticks = 2L & xminor = 5L
      linethick = 4. & symsize = 4.
   endif

   yrange = [-1.0, max(smooth(flux,10))*1.05]
   if not keyword_set(xrange) then xrange = [1000, 7100.]
   len = 0.04
   pos = [0.08, 0.5, 0.97, 0.98]  ;else pos = [0.12, 0.14, 0.96, 0.98]
   ytitle = textoidl('Flux Density f_\lambda (10^{-17} erg s^{-1} cm^{-2} \AA^{-1})')
   xtitle = textoidl('Rest Wavelength (\AA)')
   if not keyword_set(nsmooth) then nsmooth = 1L
   plot, wave, smooth(flux,nsmooth), xrange = xrange, yrange = yrange, /ystyle $
      , charsize = charsize, xthick = thick, ythick = thick, charthick = thick, _extra = extra $
      , /xsty, xticklen=len, yticklen = len/2, pos = pos, thick = linethick

   ;oplot, lam, flux, psym=2, color=fsc_color('black'), symsize=0.6
   if keyword_set(SDSS_name) then xyouts, 0.18, 0.85, /norm, SDSS_name $
     , charsize = 2.5, charthick = thick
   oplot, wave, f_conti_model, color = fsc_color('red')
   oplot, wave, f_pl_model+f_poly_model, color=fsc_color('brown')
   oplot, wave, f_fe_mgii_model, color=fsc_color('blue')
   oplot, wave, f_fe_balmer_model, color=fsc_color('blue')
   oplot, wave, f_bc_model, color=fsc_color('cyan')

   for jj=0L, (size(window_all))[2] - 1 do begin
     oplot, [window_all[0,jj], window_all[0,jj]], [0.9, 0.8]*yrange[1], color=fsc_color('gray')
     oplot, [window_all[1,jj], window_all[1,jj]], [0.9, 0.8]*yrange[1], color=fsc_color('gray')
     oplot, [window_all[0,jj], window_all[1,jj]], [0.9, 0.9]*yrange[1], color=fsc_color('gray')
   endfor

   ; now plot the windows of each emission line fitting
   ;ind = where( ((wave ge 6400. and wave le 6800.)   or $  ; Halpha
   ;              (wave ge 4700. and wave le 5100.)   or $  ; Hbeta
   ;              (wave ge 2700. and wave le 2900.)   or $  ; MgII
   ;              (wave ge 1820. and wave le 1970.)   or $  ; CIII] complex
   ;              (wave ge 1500. and wave le 1600.))     $  ; CIV
   pos = [0.84, 0.1, 0.97, 0.44] & d_pos = [0.19, 0., 0.19, 0.]
   ind = where( (wave ge 6400. and wave le 6800.), nnn )
   if nnn gt 10 then begin
      plot, wave[ind], smooth(line_flux[ind],nsmooth), /noerase, xrange=[6400., 6800.] $
        , pos = pos,/xsty,thick=linethick, charsize = charsize $
        , xthick = thick, ythick = thick, charthick = thick, xticks = xticks, xminor = xminor
      if keyword_set(fit_line) then begin
         pp = line_fit_all[0:(ngauss_broad_ha+5L)*3L-1L]
         pp_broad = line_fit_all[0:(ngauss_broad_ha)*3L-1L]
         pp_narrow = line_fit_all[(ngauss_broad_ha)*3L:(ngauss_broad_ha+5L)*3L-1L]
         for jj=1L, n_elements(pp_narrow)/3L do begin
           pp_narrow_1gauss = pp_narrow[(jj-1)*3L:jj*3L-1]
           oplot, wave[ind], onegauss(alog(wave[ind]), pp_narrow_1gauss), color=fsc_color('cyan')
         endfor
         oplot,wave[ind], manygauss(alog(wave[ind]),pp_broad,nline=n_elements(pp_broad)/3L) $
          , color=fsc_color('green')
         ; now plot total
         oplot,wave[ind],manygauss(alog(wave[ind]),pp,nline=n_elements(pp)/3L) $
          , color=fsc_color('red')
      endif
   endif
   pos = pos - d_pos & inext = (ngauss_broad_ha+5L)*3L
   ind = where( (wave ge 4700. and wave le 5100.), nnn)
   indx0 = inext
   indx1 = indx0 + (ngauss_broad_hb)*3L-1L
   indx2 = indx0 + (ngauss_broad_hb+3L)*3L-1L
   if nnn gt 10 then begin
      yrange1 = [-2, max(smooth(line_flux[ind],10)) ]
      plot, wave[ind], smooth(line_flux[ind],nsmooth), /noerase, xrange=[4700., 5100.] $
        , pos = pos, /xsty,thick=linethick, yrange = yrange1 $
        , charsize = charsize, xthick = thick, ythick = thick, charthick = thick $
        , xticks = xticks, xminor = xminor
      if keyword_set(fit_line) then begin
         pp = line_fit_all[indx0:indx2] & pp_broad = line_fit_all[indx0:indx1] 
         pp_narrow = line_fit_all[indx1+1:indx2]
         for jj=1L, n_elements(pp_narrow)/3L do begin
            pp_narrow_1gauss = pp_narrow[(jj-1)*3L:jj*3L-1]
          ;oplot, wave[ind], manygauss(alog(wave[ind]), pp_narrow, nline = n_elements(pp_narrow)/3L), color=fsc_color('cyan')
            oplot, wave[ind], onegauss(alog(wave[ind]), pp_narrow_1gauss), color=fsc_color('cyan')
         endfor
         oplot, wave[ind], manygauss(alog(wave[ind]), pp_broad, nline = n_elements(pp_broad)/3L) $
           , color=fsc_color('green')
         ; now plot total
         oplot, wave[ind], manygauss(alog(wave[ind]), pp, nline = n_elements(pp)/3L) $
           , color=fsc_color('red')
      endif
   endif
   pos = pos - d_pos & inext = indx2 + 1L
   ind = where( (wave ge 2600. and wave le 3000.), nnn)
   ind_badpix = where( wave ge 2600. and wave le 3000. and fitmask eq 0, nnn_bad)
   indx0 = inext
   indx1 = indx0 + (ngauss_broad_MgII)*3L-1L
   indx2 = indx0 + (ngauss_broad_MgII+1L)*3L-1L
   if nnn gt 10 then begin
      plot, wave[ind], smooth(line_flux[ind],nsmooth), /noerase, xrange=[2600., 3000.], pos = pos $
        , /xsty,thick=linethick $
        , charsize = charsize, xthick = thick, ythick = thick, charthick = thick $
        , xticks = xticks, xminor = xminor
      if nnn_bad gt 0 then oplot, [wave[ind_badpix]], [line_flux[ind_badpix]] $
        , psym=4, thick =thick, color=fsc_color('cyan')
      if keyword_set(fit_line) then begin
         pp = line_fit_all[indx0:indx2] & pp_broad = line_fit_all[indx0:indx1] 
         pp_narrow = line_fit_all[indx1+1:indx2]
         oplot, wave[ind], manygauss(alog(wave[ind]), pp_narrow), color=fsc_color('cyan')
         oplot, wave[ind], manygauss(alog(wave[ind]), pp_broad), color=fsc_color('green')
         oplot, wave[ind], manygauss(alog(wave[ind]), pp), color=fsc_color('red')
      endif
   endif
   pos = pos - d_pos & inext = indx2 + 1L
   ind = where( (wave ge 1800. and wave le 2000.), nnn)
   ind_badpix = where( wave ge 1800. and wave le 2000. and fitmask eq 0, nnn_bad)
   indx0 = inext
   indx1 = indx0 + (ngauss_broad_CIII)*3L-1L
   indx2 = indx0 + (ngauss_broad_CIII+1L)*3L-1L
   indx3 = indx0 + (ngauss_broad_CIII+3L)*3L-1L
   if nnn gt 10 then begin
      plot, wave[ind], smooth(line_flux[ind],nsmooth), /noerase, xrange=[1800., 2000.] $
         , pos = pos, /xsty,thick=linethick $
         , charsize = charsize, xthick = thick, ythick = thick, charthick = thick $
         , xticks = xticks, xminor = xminor
      if nnn_bad gt 0 then oplot, [wave[ind_badpix]], [line_flux[ind_badpix]] $
         , psym=4, thick =thick, color=fsc_color('cyan')
      if keyword_set(fit_line) then begin
         pp = line_fit_all[indx0:indx3] & pp_broad = line_fit_all[indx0:indx1] 
         pp_narrow = line_fit_all[indx1+1:indx2]
         pp_Si_Al = line_fit_all[indx2+1:indx3]
         oplot, wave[ind], manygauss(alog(wave[ind]), pp_narrow), color=fsc_color('cyan')
         oplot, wave[ind], manygauss(alog(wave[ind]), pp_broad), color=fsc_color('green')
         oplot, wave[ind], manygauss(alog(wave[ind]), pp_Si_Al), color=fsc_color('magenta')
         oplot, wave[ind], manygauss(alog(wave[ind]), pp), color=fsc_color('red')
      endif
   endif
   pos = pos - d_pos & inext = indx3 + 1L
   ind = where( wave ge 1500. and wave le 1700., nnn )
   ind_badpix = where( wave ge 1500. and wave le 1700. and fitmask eq 0, nnn_bad)
   indx0 = inext
   indx1 = indx0 + (ngauss_broad_CIV)*3L-1L
   indx2 = indx0 + (ngauss_broad_CIV+1L)*3L-1L
   if nnn gt 10 then begin
      plot, wave[ind], smooth(line_flux[ind],nsmooth), /noerase, xrange=[1500., 1700.], pos = pos $
       , /xsty,thick=linethick, charsize = charsize, xthick = thick, ythick = thick $
       , charthick = thick, xticks = xticks, xminor = xminor
      if nnn_bad gt 0 then oplot, [wave[ind_badpix]], [line_flux[ind_badpix]] $
       , psym=4, thick =thick, color=fsc_color('cyan')
      if keyword_set(fit_line) then begin
         pp = line_fit_all[indx0:indx2] & pp_broad = line_fit_all[indx0:indx1] 
         pp_narrow = line_fit_all[indx1+1:indx2]
         oplot, wave[ind], manygauss(alog(wave[ind]), pp_narrow), color=fsc_color('cyan')
         oplot, wave[ind], manygauss(alog(wave[ind]), pp_broad), color=fsc_color('green')
         oplot, wave[ind], manygauss(alog(wave[ind]), pp), color=fsc_color('red')
         if keyword_set(HeII_fit) then begin
            pp = line_fit_all[indx1+1+3:indx2+3]
            oplot, wave[ind],onegauss(alog(wave[ind]),pp),color=fsc_color('magenta')
            pp = line_fit_all[indx1+1+6:indx2+6]
            oplot, wave[ind],onegauss(alog(wave[ind]),pp),color=fsc_color('magenta')
         endif
      endif
   endif

   xyouts, 0.5, 0.01, xtitle, /norm, charsize = charsize, charthick = thick, align = 0.5
   xyouts, 0.02, 0.5, ytitle, /norm, charsize = charsize, charthick = thick, align =0.5, orien = 90

   ;xyouts, 0.7, 0.84, objtag, /norm, charsize = charsize, charthick = thick
   xyouts, 0.7, 0.8, 'z='+string(z,format='(f5.3)'), /norm, charsize=charsize, charthick=thick

   para = output

   if keyword_set(psplot) then endplot

end
