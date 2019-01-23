;------------------------------------------------------------------------
function gauss_conv, x, y, sig_conv = sig_conv
; x should be ln(wave), y -- flux
; this works much slower than the one used in get_Fe_flux

  dx = djs_mean(([x, 0] - [0, x])[1:n_elements(x)-1])

  y_conv = dblarr(n_elements(y))

  for i=0, n_elements(y) - 1L do begin
    gauss = exp(-(x - x[i])^2/2./sig_conv^2)/sig_conv/sqrt(2*!PI)
    y_conv[i] = total(y*gauss*dx)/total(gauss*dx)
  end

return, y_conv
end
;------------------------------------------------------------------------
function get_Fe_flux_hbeta, xval, Fe_FWHM
common Fe_temp, wave_Fe, flux_Fe  ; common block for the original Fe template

c = 2.9979246d5

; convolve the original Fe template with parameter Fe_FWHM (I Zw1 has intrinsic FWHM 900 km/s);
; initial FWHM in the template is 900 km/s, so force Fe_FWHM=910.0 if input Fe_FWHM <= 900 km/s
if Fe_FWHM le 900.0 then sig_conv = sqrt(910.0^2 - 900.0^2)/2./sqrt(2.*alog(2.)) $
else sig_conv = sqrt(Fe_FWHM^2 - 900.0^2)/2./sqrt(2.*alog(2.))  ; in units of km/s

; get sigma in pixcel space
sig_pix = sig_conv/106.3  ; where 106.3 km/s is the dispersion for the Boroson Fe template

khalfsz = round (4*sig_pix+1)
xx= findgen(khalfsz*2+1) - khalfsz
kernel = exp(-xx^2/(2*sig_pix^2))
kernel = kernel/total(kernel)

flux_Fe_conv = convol(flux_Fe, kernel, /center, /edge_truncate)

yval = spline(wave_Fe, flux_Fe_conv, xval)

return, yval
end

;--------------------------------------------------------------------------
function conti_Fe_hbeta, xval, pp
; pp[0]: norm_factor for continuum f_lambda = (lambda/3000.0)^{-alpha}
; pp[1]: slope for the power-law continuum
; pp[2]: norm_factor for the Fe_template
; pp[3]: FWHM for the Fe_template
; pp[4]: small shift of wavelength for the Fe template

yval = pp[0]*(xval/3000.0)^pp[1] + pp[2]*get_Fe_flux_hbeta(xval*(1.0 + pp[4]), pp[3])

return, yval

end
;------------------------------------------------------------------

pro Hbetafit_with_Fe, plate, fiber, mjd, z, psplot = psplot, noplot=noplot, para = para $
                    , fits = fits, SDSS_name = SDSS_name $
                    , deredden = deredden, ra = ra, dec = dec $ ; for dereddening purpose 
                    , spec_dir = spec_dir, QA_dir = QA_dir, out_dir = out_dir $
                    , ngauss_broad_Hb = ngauss_broad_Hb $ ; number of gaussians for the broad Hb line
                    , guess = guess, _extra=extra $ ; initial guess for the Hbeta + OIII fits
                    , xrange=xrange, oiii_plot=oiii_plot, silent=silent, print_para = print_para $
                    , gh = gh $ ; using Gaussian-Hermite function for each line instead 
                    , iterate = iterate $ ; add upto 3 Gaussians for the broad Hbeta
                    , double_oiii = double_oiii $ ; set this keyword to use two-gaussian for each [OIII] to account for blue wing and double-peaked [OIII]
                    , figfile_name = figfile_name, outfile_name = outfile_name, flux=flux, err=err, lam=lam $
                    , force = force $ ; proceed to linefit regardless of redchi^2 for the continuum fit
                    , max_sigma = max_sigma, no_plate_subdir = no_plate_subdir, dup_spec=dup_spec, diet = diet, conti_slope_limit = conti_slope_limit, no_iron = no_iron $
                    , oiii_hb_ratio = oiii_hb_ratio, input_fe_fwhm = input_fe_fwhm, conti_norm_guess = conti_norm_guess $
                    , fit_HeII = fit_HeII ; set to fit HeII 4687.02

max_ngauss_broad = 3L  ; maximun gaussians used for the broad Hbeta

if n_elements(z) eq 0 then z=0.

if keyword_set(iterate) and keyword_set(gh) then begin
  splog, 'ERROR: keywords iterate only works without gh'
  return
endif

if not keyword_set(max_sigma) then max_sigma = 0.05D  ; default maxmium dispersion for the broad Gaussians

; para struct contains fitted paramters and errors
if not keyword_set(diet) then begin
para = create_struct('zsys', z, 'conti_norm', 0.0D, 'conti_norm_err', -1.0D $
                   , 'conti_slope', 0.0D, 'conti_slope_err', -1.0D $
                   , 'Fe_norm', 0.0D, 'Fe_norm_err', -1.0D $
                   , 'Fe_FWHM', 0.0D, 'Fe_FWHM_err', -1.0D $
                   , 'Hbeta1_wave', 0.0D, 'Hbeta2_wave', 0.0D $
                   , 'OIII4959_wave', 0.0D, 'OIII5007_wave', 0.0D $
                   , 'Hbeta1_linesigma', 0.0D, 'Hbeta1_linesigma_err', -1.0D $
                   , 'Hbeta1_FWHM', 0.0D, 'Hbeta1_FWHM_err', -1.0D $
                   , 'Hbeta1_linearea', 0.0D, 'Hbeta1_linearea_err', -1.0D $
                   , 'Hbeta1_lineew', 0.0D, 'Hbeta1_lineew_err', -1.0D $
                   , 'Hbeta2_linesigma', 0.0D, 'Hbeta2_linesigma_err', -1.0D $
                   , 'Hbeta2_FWHM', 0.0D, 'Hbeta2_FWHM_err', -1.0D $
                   , 'Hbeta2_linearea', 0.0D, 'Hbeta2_linearea_err', -1.0D $
                   , 'Hbeta2_lineew', 0.0D, 'Hbeta2_lineew_err', -1.0D $
                   , 'OIII4959_linesigma', 0.0D, 'OIII4959_linesigma_err', -1.0D $
                   , 'OIII4959_FWHM', 0.0D, 'OIII4959_FWHM_err', -1.0D $
                   , 'OIII4959_linearea', 0.0D, 'OIII4959_linearea_err', -1.0D $
                   , 'OIII4959_lineew', 0.0D, 'OIII4959_lineew_err', -1.0D $
                   , 'OIII5007_linesigma', 0.0D, 'OIII5007_linesigma_err', -1.0D $
                   , 'OIII5007_FWHM', 0.0D, 'OIII5007_FWHM_err', -1.0D $
                   , 'OIII5007_linearea', 0.0D, 'OIII5007_linearea_err', -1.0D $
                   , 'OIII5007_lineew', 0.0D, 'OIII5007_lineew_err', -1.0D $
                   , 'redchi2_1', 0.0D, 'redchi2_2', 0.0D, 'status_1', 0L, 'status_2', 0L, 'line_npix', 0L $
                   , 'line_med_SN', 0.D, 'ngauss_broad_Hb', 0L)
endif else para = {zsys:z, redchi2_1:0.D, redchi2_2:0.D, status_1:0L, status_2:0L, line_npix:0L, line_med_SN:0.D, ngauss_broad_Hb:0L}

c = 2.9979246d5

if not keyword_set(spec_dir) then prefix1 = '/data1/quasar/yshen/DR7_spectra/' $
else prefix1 = spec_dir ; 1d spectra for DR7 QSOs
if not keyword_set(QA_dir) then prefix2 = '/peyton/scr/chimera0/yshen/DR7_QSO_fits/Hbeta/QA/' $
else prefix2 = QA_dir ; QA plots
if not keyword_set(out_dir) then prefix3 = '/peyton/scr/chimera0/yshen/DR7_QSO_fits/Hbeta/fits/' $
else prefix3 = out_dir  ; fitting results

if not keyword_set(ngauss_broad_Hb) then ngauss_broad_Hb = 1L ; default using a single Gaussian for the broad Hbeta line

if n_elements(flux) eq 0 then begin  ; using plate-fiber-mjd to locate the spectrum, unless the flux is specified
platestr = string(plate, format='(i4.4)')
mjdstr = string(mjd, format='(i5.5)')
fibstr = string(fiber, format='(i3.3)')
if not keyword_set(dup_spec) then infile =  prefix1 + platestr + '/' + platestr + '-' + fibstr + '-' + mjdstr + '.dat' $
else infile =  prefix1 + platestr + '-' + fibstr + '-' + mjdstr + '.dat'
if not keyword_Set(no_plate_subdir) then begin
  figfile = prefix2 + platestr + '/' + platestr + '-' + fibstr + '-' + mjdstr + '.eps'
  fitsfile = prefix3 + platestr + '/' + platestr + '-' + fibstr + '-' + mjdstr + '.fits'
endif else begin
  figfile  = prefix2 + platestr + '-' + fibstr + '-' + mjdstr + '.eps'
  fitsfile = prefix3 + platestr + '-' + fibstr + '-' + mjdstr + '.fits'
endelse

file_found = findfile(infile, count = n_file)
if n_file eq 0L then begin
  readspec, plate, fiber, mjd = mjd, wave = lam, flux = flux, invvar = ivar, /silent
endif else begin
  readcol, infile, format = 'd,d,d', lam, flux, ivar, /silent
endelse
err = dblarr(n_elements(lam))
indd = where(ivar gt 1e-6)
err[indd] = 1.0/sqrt(ivar[indd])

; deredden Galactic extinction if asked
if keyword_set(deredden) then begin
  if not keyword_set(silent) then splog, 'dereddening spectrum using the SFD map and CCM extinction curve'
  dereddening_spec, lam, flux, err=err, ra=ra, dec=dec, dered_flux = dered_flux, dered_err = dered_err
  flux = dered_flux & err = dered_err
endif

endif else begin ; input directly the 1d spectrum
 ivar = dblarr(n_elements(err))
 indd = where(err gt 1d-6)
 ivar[indd] = 1./(err[indd])^2
 if keyword_set(psplot) then figfile = figfile_name
 if keyword_set(fits) then fitsfile = outfile_name
endelse

; now scale to rest wavelength
if not keyword_set(silent) then begin
  splog, 'shift to restframe:'
endif
lam = lam/(1.0 + z)

;plot, lam, flux

;+++++++++++++++++++++++++++ Continuum + Fe fit ++++++++++++++++++++++++++++++++++++
; wave range for fitting Fe + continuum
window1 = [4435.0, 4700.0]
window2 = [5100.0, 5535.0]

common Fe_temp, wave_Fe, flux_Fe   ; common block for the original Fe template

if n_elements(wave_Fe) le 10 then begin
;file = '/u/yshen/Research/IDL/lib/Projects/quasar/linefit/Fe_fit/irontemplate.dat'
file = '/home/yshen/Research/IDL/lib/Projects/quasar/linefit/Fe_fit/irontemplate.dat'
; wave_Fe1 [2200, 3090]; wave_Fe2 [4526, 6357]
readcol, file, format = 'd,d', logwave_1, flux_1
wave_Fe = 10.0D^logwave_1
flux_Fe = flux_1*1d15
ind = where(wave_Fe ge 4400 and wave_Fe le 5600)
wave_Fe = wave_Fe[ind]
flux_Fe = flux_Fe[ind]
endif

ind = where( ((lam ge window1[0] and lam le window1[1]) or (lam ge window2[0] and lam le window2[1])) $
             and err gt 1e-6)

if n_elements(ind) gt 5 then begin

; set limiting conditions on paramters
parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0]}, 5)
if not keyword_set(conti_norm_guess) then conti_norm_guess = 1.D
parinfo.value = [conti_norm_guess, -2.0D, 1.0D, 3000.0D, 0.0]
parinfo[0].limited = [1,0]
parinfo[0].limits = [0., 1d10]
parinfo[1].LIMITED = [1, 1]
if not keyword_set(conti_slope_limit) then parinfo[1].LIMITS = [-5.0,3.0] else parinfo[1].LIMITS = conti_slope_limit
parinfo[2].limited = [1, 0]
parinfo[2].limits = [0., 1d10]
parinfo[3].LIMITED = [1, 1]
parinfo[3].LIMITS = [1000.0, 1d4]  ; [500.0, 2d4]
parinfo[4].LIMITED = [1, 1]
parinfo[4].LIMITS = [-0.005, 0.005]

; fix the iron template FWHM to the input value
if keyword_set(input_fe_fwhm) then begin
  parinfo[3].value = input_fe_fwhm
  parinfo[3].fixed = 1L
endif

if keyword_set(no_iron) then begin ; do not fit iron template
  parinfo[2].value = 0.D & parinfo[2].fixed = 1L
endif

conti_fit = mpfitfun('conti_Fe_hbeta', lam[ind], flux[ind], err[ind] $
                     , parinfo = parinfo, perror = perror, yfit = yfit, /quiet $
                     , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
if not keyword_set(silent) then splog, 'FE+CONTI FITTING: MPFIT nfev=', nfev, ' niter=', niter, ' status=', status
;splog, conti_fit

temp = create_struct('conti_fit', conti_fit, 'conti_fit_err', perror)
para = struct_addtags(para, temp)

if not keyword_set(diet) then begin
  para.conti_norm = conti_fit[0]
  para.conti_norm_err = perror[0]
  para.conti_slope = conti_fit[1]
  para.conti_slope_err = perror[1]
  para.Fe_norm = conti_fit[2]
  para.Fe_norm_err = perror[2]
  para.Fe_FWHM = conti_fit[3]
  para.Fe_FWHM_err = perror[3]
endif
para.redchi2_1 = chi2/(n_elements(ind) - 5.0)
para.status_1 = status

; if the Fe + power-law continuum failed, try to use a power-law contiuum alone
if para.redchi2_1 gt 10.0 or para.redchi2_1 lt 1d-5 then begin
  parinfo[2:4].fixed = 1L
  parinfo[2:4].value = [0.0D, 3000.0D, 0.0] ; force Fe flux=0
  conti_fit_new = mpfitfun('conti_Fe_hbeta', lam[ind], flux[ind], err[ind] $
                     , parinfo = parinfo, perror = perror, yfit = yfit, /quiet $
                     , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
if not keyword_set(silent) then  splog, 'POWER-LAW CONTI ONLY: MPFIT nfev=', nfev, ' niter=', niter, ' status=', status
  if chi2/(n_elements(ind) - 2.0) lt para.redchi2_1 then begin ; replace the previous Fe+conti fit
     conti_fit = conti_fit_new
     if not keyword_set(diet) then begin
       para.conti_norm = conti_fit[0]
       para.conti_norm_err = perror[0]
       para.conti_slope = conti_fit[1]
       para.conti_slope_err = perror[1]
       para.Fe_norm = conti_fit[2]
       para.Fe_norm_err = -1.
       para.Fe_FWHM = conti_fit[3]
       para.Fe_FWHM_err = -1.
     endif
     para.redchi2_1 = chi2/(n_elements(ind) - 2.0)
     para.status_1 = status
     para.conti_fit = conti_fit & para.conti_fit_err = perror
  endif
endif

;+++++++++++++++++++++++++++ Hbeta + [OIII] line fit ++++++++++++++++++++++++++++++++++++
; if redchi2_1 < 10 then proceed to fit the lines, otherwise return
; now proceed to fit the gaussians for Hbeta (narrow + broad), [OIII]4959, [OIII]5007
; substract the continuum and Fe
if keyword_set(force) then redchi2_cut = 1d10 else $ ; no cut on redchi2
redchi2_cut = 10.0

if para.redchi2_1 le redchi2_cut and para.redchi2_1 ge 1d-5 then begin
window3 = [4700.0, 5100.0] ; fitting range for the Hbeta + [OIII] doublets
if keyword_set(fit_HeII) then window3 = [4600, 5100] ; enlarge the fitting range by including HeII
; recall the fitting windows for Fe + contiuum is 
; window1 = [4435.0, 4700.0]; window2 = [5100.0, 5535.0]
ind = where(lam ge 4435.0 and lam le 5535.0 and err gt 1e-6)

if ind[0] ne -1 then begin
wave_fit = lam[ind]
conti_flux = conti_fit[0]*(wave_fit/3000.0)^conti_fit[1]
Fe_flux = conti_fit[2]*get_Fe_flux_hbeta(wave_fit*(1.0 + conti_fit[4]), conti_fit[3])
flux_fit = flux[ind] - conti_flux - Fe_flux  ; continuum-subtracted flux
err_fit = err[ind]


temp = create_struct('wavefit', wave_fit, 'conti_flux', conti_flux $
                     , 'Fe_flux', Fe_flux)
if not keyword_set(diet) then para = struct_addtags(para, temp)
para.ngauss_broad_hb = ngauss_broad_Hb

if not keyword_set(double_oiii) then ngauss_narrow_line = 3L else ngauss_narrow_line = 5L

if not keyword_set(gh) then begin
; fit the gaussians for Hbeta + [OIII]
parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0],tied:''}, 3*(ngauss_narrow_line+ngauss_broad_Hb))
functargs = {nline: ngauss_narrow_line+ngauss_broad_Hb}
if keyword_set(fit_HeII) then begin ; add two gaussians (one broad and one narrow) for HeII4687.02
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0],tied:''}, 3*(ngauss_narrow_line+ngauss_broad_Hb+2L))
  functargs = {nline: ngauss_narrow_line+ngauss_broad_Hb+2L}
endif
if ngauss_broad_Hb eq 1L then begin  ; assuming 4 gaussians (i.e., only 1 gaussian for the broad Hbeta)
; the Hbeta broad component
parinfo[0:2].value = [1.0D, alog(4862.68D), 0.01D]
parinfo[0].limited = [1L, 0L]
parinfo[0].limits = [0.0, 0.0]
parinfo[1].limited = [1L,1L]
parinfo[1].limits = [alog(4862.68) - 0.015D, alog(4862.68) + 0.015D]
parinfo[2].limited = [1L, 1L]
parinfo[2].limits = [0.0017D, max_sigma] ; restrict the broad component to have 1400<FWHM<35322 km/s
; the Hbeta narrow component (FWHM < 1200 km/s)
parinfo[3:5].value = [0.25D, alog(4862.68D), 0.001D]  ; [0.5D, alog(4862.68D), 0.001D]
parinfo[3].limited = [1L, 0L]
parinfo[3].limits = [0.0, 0.0]
parinfo[4].limited = [1L,1L]
parinfo[4].limits = [alog(4862.68) - 0.015D, alog(4862.68) + 0.015D] ; narrow offset [-1000,1000]km/s
parinfo[5].limited = [1L,1L]
parinfo[5].limits = [1d-4, 0.0017D] ; restrict the narrow component to have FWHM<1200 km/s, 1d-4
; the OIII 4959 and OIII 5007
parinfo[6:8].value = [0.25D, alog(4960.30), 0.001D]
parinfo[6].limited = [1L, 0L]
parinfo[6].limits = [0., 1d10]  ; emission line only
parinfo[7].tied = '(P(4) - alog(4862.68D)) + alog(4960.30D)'
parinfo[8].tied = 'P(5)'
; the OIII 5007
parinfo[9:11].value = [0.5D, alog(5008.24D), 0.001D]
parinfo[9].limited = [1L, 0L]
parinfo[9].limits = [0., 1d10]  ; emission line only
parinfo[10].tied = '(P(4) - alog(4862.68D)) + alog(5008.24D)'
parinfo[11].tied = 'P(5)'
; HeII 4687.02
if not keyword_set(double_oiii) and keyword_set(fit_HeII) then begin
  parinfo[12:14].value = [0.1D, alog(4687.02D), 0.001D]
  parinfo[12].limited = [1L,0L] & parinfo[12].limits = [0., 1d10]
  parinfo[13].limited = [1L,1L] & parinfo[13].limits = [alog(4687.02D) - 2d-3, alog(4687.02D) + 2d-3] ; within +-600km/s
  parinfo[14].limited = [1L,1L] & parinfo[14].limits = [2d-4, 0.0017D]
  parinfo[15:17].value = [0.1D, alog(4687.02D), 0.01D]
  parinfo[15].limited = [1L,0L] & parinfo[15].limits = [0., 1d10]
  parinfo[16].limited = [1L,1L] & parinfo[16].limits = [alog(4687.02D) - 2d-3, alog(4687.02D) + 2d-3] ; within +-600km/s
  parinfo[17].limited = [1L,1L] & parinfo[17].limits = [0.0017D, 0.015D]
endif

if keyword_set(double_Oiii) then begin
  ; additional two gaussians for the blue wings of the OIII doublets
  parinfo[12:14].value = [0.1D, alog(4960.30) - 1d-3, 0.001D]
  parinfo[12].limited = [1L, 0L]
  parinfo[12].limits = [0., 1d10]  ; emission line only
  parinfo[13].limited = [1L,1L]
  parinfo[13].limits = [alog(4960.30) - 3.3d-3, alog(4960.30)]
  parinfo[14].limited = [1L,1L]
  parinfo[14].limits = [1d-4, 0.0017D]
  parinfo[15:17].value = [0.2D, alog(5008.24D), 0.001D]
  parinfo[15].limited = [1L, 0L]
  parinfo[15].limits = [0., 1d10]  ; emission line only
  parinfo[16].tied = '(P(13) - alog(4960.30D)) + alog(5008.24D)'
  parinfo[17].tied = 'P(14)'
  if keyword_set(fit_HeII) then begin
    ; HeII 4687.02
    parinfo[18:20].value = [0.1D, alog(4687.02D), 0.001D]
    parinfo[18].limited = [1L,0L] & parinfo[18].limits = [0., 1d10]
    parinfo[19].limited = [1L,1L] & parinfo[19].limits = [alog(4687.02D) - 2d-3, alog(4687.02D) + 2d-3] ; within +-600km/s
    parinfo[20].limited = [1L,1L] & parinfo[20].limits = [2d-4, 0.0017D]
    parinfo[21:23].value = [0.1D, alog(4687.02D), 0.01D]
    parinfo[21].limited = [1L,0L] & parinfo[21].limits = [0., 1d10]
    parinfo[22].limited = [1L,1L] & parinfo[22].limits = [alog(4687.02D) - 2d-3, alog(4687.02D) + 2d-3] ; within +-600km/s
    parinfo[23].limited = [1L,1L] & parinfo[23].limits = [0.0017D, 0.015D]
  endif
endif
endif else begin  ; if using more than one gaussians for the broad Hbeta
  for i_gau=0L, ngauss_broad_Hb - 1L do begin
    parinfo[0+i_gau*3L:2+i_gau*3L].value = [1.0D/ngauss_broad_Hb, alog(4862.68D), 0.005D]
    parinfo[0+i_gau*3L].limited = [1L, 0L]
    parinfo[0+i_gau*3L].limits = [0.0, 0.0]
    parinfo[1+i_gau*3L].limited = [1L,1L]
    parinfo[1+i_gau*3L].limits = [alog(4862.68) - 0.015D, alog(4862.68) + 0.015D]
    parinfo[2+i_gau*3L].limited = [1L, 1L]
    parinfo[2+i_gau*3L].limits = [0.0017D, max_sigma] ; restrict each component to have 1400<FWHM<35322 km/s
    ; restrict the centroid of each broad component is the same, hence the overall broad line is symmetric
    ; if i_gau gt 0L then parinfo[1+i_gau*3L].tied = 'P(1)'
  endfor
  inext = (ngauss_broad_Hb - 1L)*3
  ; the Hbeta narrow component (FWHM < 1200 km/s)
  parinfo[3+inext:5+inext].value = [0.25D, alog(4862.68D), 0.001D]
  parinfo[3+inext].limited = [1L, 0L]
  parinfo[3+inext].limits = [0.0, 0.0]
  parinfo[4+inext].limited = [1L,1L]
  parinfo[4+inext].limits = [alog(4862.68) - 0.015D, alog(4862.68) + 0.015D]
  parinfo[5+inext].limited = [1L,1L]
  parinfo[5+inext].limits = [1d-4, 0.0017D] ; restrict the narrow component to have FWHM<1200 km/s
  ; the OIII 4959
  parinfo[6+inext:8+inext].value = [0.25D, alog(4960.30D), 0.001D]
  parinfo[6+inext].limited = [1L, 0L]
  parinfo[6+inext].limits = [0., 1d10]
  pstr1 = '(P(' + string(4+inext, format='(i0)')+')'
  parinfo[7+inext].tied = pstr1 +  ' - alog(4862.68D)) + alog(4960.30D)'
  pstr2 = 'P(' + string(5+inext, format='(i0)') + ')'
  parinfo[8+inext].tied = pstr2
  ; the OIII 5007
  parinfo[9+inext:11+inext].value = [0.5D, alog(5008.24D), 0.001D]
  parinfo[9+inext].limited = [1L, 0L]
  parinfo[9+inext].limits = [0., 1d10]
  parinfo[10+inext].tied = pstr1 + ' - alog(4862.68D)) + alog(5008.24D)'
  parinfo[11+inext].tied = pstr2
  if keyword_set(oiii_hb_ratio) then begin ; fix the oiii5007/hb flux ratio
    ratio = oiii_hb_ratio*4862.68D/5008.24D
    pstr_oiii_hb = '(P(' + string(3+inext, format='(i0)')+')*' + string(ratio,format='(f0)')+')'
    parinfo[9+inext].tied = pstr_oiii_hb
  endif
  if not keyword_set(double_oiii) and keyword_set(fit_HeII) then begin
    ; HeII 4687.02
    parinfo[12+inext:14+inext].value = [0.1D, alog(4687.02D), 0.001D]
    parinfo[12+inext].limited = [1L,0L] & parinfo[12+inext].limits = [0., 1d10]
    parinfo[13+inext].limited = [1L,1L] & parinfo[13+inext].limits = [alog(4687.02D) - 2d-3, alog(4687.02D) + 2d-3] ; within +-600km/s
    parinfo[14+inext].limited = [1L,1L] & parinfo[14+inext].limits = [2d-4, 0.0017D]
    parinfo[15+inext:17+inext].value = [0.1D, alog(4687.02D), 0.01D]
    parinfo[15+inext].limited = [1L,0L] & parinfo[15+inext].limits = [0., 1d10]
    parinfo[16+inext].limited = [1L,1L] & parinfo[16+inext].limits = [alog(4687.02D) - 2d-3, alog(4687.02D) + 2d-3] ; within +-600km/s
    parinfo[17+inext].limited = [1L,1L] & parinfo[17+inext].limits = [0.0017D, 0.015D]
  endif

  if keyword_set(double_oiii) then begin
    ; additional two gaussians for the blue wings of the OIII doublets
    parinfo[12+inext:14+inext].value = [0.1D, alog(4960.30) - 1d-3, 0.001D]
    parinfo[12+inext].limited = [1L, 0L]
    parinfo[12+inext].limits = [0., 1d10]  ; emission line only
    parinfo[13+inext].limited = [1L,1L]   ; limted was incorrectly used -- corrected 08/02/2010; note this doesn't matter as I was using keywords iterate+ngauss_broad_hb=1L
    parinfo[13+inext].limits = [alog(4960.30) - 3.3d-3, alog(4960.30)]
    parinfo[14+inext].limited = [1L,1L]
    parinfo[14+inext].limits = [1d-4, 0.0017D]
    parinfo[15+inext:17+inext].value = [0.2D, alog(5008.24D), 0.001D]
    parinfo[15+inext].limited = [1L, 0L]
    parinfo[15+inext].limits = [0., 1d10]  ; emission line only
    pstr3 = '(P(' + string(13+inext, format='(i0)') + ')'  ; added the first '(' -- corrected 08/02/2010; note this doesn't matter as I was using keywords iterate+ngauss_broad_hb=1L
    parinfo[16+inext].tied = pstr3 + ' - alog(4960.30D)) + alog(5008.24D)'
    parinfo[17+inext].tied = 'P(' + string(14+inext, format='(i0)') + ')'

    if keyword_set(fit_HeII) then begin
      ; HeII 4687.02
      parinfo[18+inext:20+inext].value = [0.1D, alog(4687.02D), 0.001D]
      parinfo[18+inext].limited = [1L,0L] & parinfo[18+inext].limits = [0., 1d10]
      parinfo[19+inext].limited = [1L,1L] & parinfo[19+inext].limits = [alog(4687.02D) - 2d-3, alog(4687.02D) + 2d-3] ; within +-600km/s
      parinfo[20+inext].limited = [1L,1L] & parinfo[20+inext].limits = [2d-4, 0.0017D]
      parinfo[21+inext:23+inext].value = [0.1D, alog(4687.02D), 0.01D]
      parinfo[21+inext].limited = [1L,0L] & parinfo[21+inext].limits = [0., 1d10]
      parinfo[22+inext].limited = [1L,1L] & parinfo[22+inext].limits = [alog(4687.02D) - 2d-3, alog(4687.02D) + 2d-3] ; within +-600km/s
      parinfo[23+inext].limited = [1L,1L] & parinfo[23+inext].limits = [0.0017D, 0.015D]
    endif
  endif
endelse
endif else begin; end assigning initial parameters for the Gaussian lines

; assign intial guess for the Gaussian-Hermite function fit, four GH function for broad Hbeta, narrow Hbeta, OIII4959 and OIII5007
parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0],tied:''}, 5*4L)
functargs = {nline:4L}
; for broad Hbeta
parinfo[0:4].value = [1.0D, alog(4862.68D), 0.01D, 0.D, 0.D]
parinfo[0].limited = [1L, 0L]
parinfo[0].limits = [0.0, 0.0]
parinfo[1].limited = [1L,1L]
parinfo[1].limits = [alog(4862.68D) - 0.015D, alog(4862.68D) + 0.015D]
parinfo[2].limited = [1L, 1L]
parinfo[2].limits = [0.0017D, max_sigma] ; restrict the broad component to have 1400<FWHM<35300 km/s
;parinfo[3].limited = [1L, 1L] & parinfo[3].limits = [-1., 1.]
;parinfo[4].limited = [1L, 1L] & parinfo[4].limits = [-1., 1.]
; for narrow Hbeta
parinfo[5:9].value = [0.5D, alog(4862.68D), 0.001D, 0.D, 0.D]
parinfo[5].limited = [1L, 0L]
parinfo[5].limits = [0.0, 0.0]
parinfo[6].limited = [1L,1L]
parinfo[6].limits = [alog(4862.68D) - 0.015D, alog(4862.68D) + 0.015D]
parinfo[7].limited = [1L, 1L]
parinfo[7].limits = [1d-4, 0.0017D]  ; restrict the narrow component to have FWHM<1200 km/s
;parinfo[8].limited = [1L, 1L] & parinfo[8].limits = [-1., 1.]
;parinfo[9].limited = [1L, 1L] & parinfo[9].limits = [-1., 1.]
; for OIII4959
parinfo[10:14].value = [0.5D, alog(4960.30D), 0.001D, 0.D, 0.D]
parinfo[10].limited = [1L, 0L]
parinfo[10].limits = [0., 1d10]  ; emission line only
parinfo[11].tied = '(P(6) - alog(4862.68D)) + alog(4960.30D)'
parinfo[12].tied = 'P(7)'
parinfo[13].tied = 'P(8)'
parinfo[14].tied = 'P(9)'
; for OIII5007
parinfo[15:19].value = [1.0D, alog(5008.24D), 0.001D, 0.D, 0.D]
parinfo[15].limited = [1L, 0L]
parinfo[15].limits = [0., 1d10]  ; emission line only
parinfo[16].tied = '(P(6) - alog(4862.68D)) + alog(5008.24D)'
parinfo[17].tied = 'P(7)'
parinfo[18].tied = 'P(8)'
parinfo[19].tied = 'P(9)'
endelse

; if specified fitting conditions for the Hbeta + OIII fits, then override parinfo
if keyword_set(guess) then parinfo = guess

; restrict the fitting range for lines to be within window3 = [4700, 5100]
indd = where(wave_fit ge window3[0] and wave_fit le window3[1])

if n_elements(indd) gt n_elements(parinfo) then begin

if not keyword_set(gh) then begin
  line_fit = mpfitfun('manygauss', alog(wave_fit[indd]), flux_fit[indd], err_fit[indd] $
                     , functargs = functargs, parinfo = parinfo, perror = perror, yfit = yfit2, /quiet $
                     , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
endif else begin
  line_fit = mpfitfun('manygh', alog(wave_fit[indd]), flux_fit[indd], err_fit[indd] $
                     , functargs = functargs, parinfo = parinfo, perror = perror, yfit = yfit2, /quiet $
                     , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
endelse
; test if the blueshift OIII component is indeed blueshifted wrt to the core oiii component
if keyword_set(double_oiii) then begin
  if line_fit[13+inext] gt line_fit[7+inext] then begin
    parinfo[13+inext].limits = [line_fit[7+inext] - 5.d-3, line_fit[7+inext]]
    parinfo.value = line_fit 
    parinfo[13+inext].value = line_fit[7+inext] - 2.d-3
    line_fit = mpfitfun('manygauss', alog(wave_fit[indd]), flux_fit[indd], err_fit[indd] $
                     , functargs = functargs, parinfo = parinfo, perror = perror, yfit = yfit2, /quiet $
                     , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
  endif

endif


if not keyword_set(silent) then splog, 'LINE-FITTING: MPFIT nfev=', nfev, ' niter=', niter, ' status=', status, ' chi2=',chi2
;splog, 'Hbeta FWHM=', line_fit[2]*2.3548D*c, line_fit[5]*2.3548D*c

if keyword_set(iterate) then begin  ; iterate using additional gaussians for the broad Hbeta
  nline = 1L + n_elements(parinfo)/3L ; 2 gaussian for broad Hb, and 3 gaussian for narrow Hb and OIII4959, 5007
  functargs_new = {nline:nline}
  parinfo_new = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0],tied:''}, 3*nline)
  iend = n_elements(parinfo_new) - 1L
  parinfo_new[3:iend] = parinfo
  parinfo_new[3:iend].value = line_fit & parinfo_new[6:iend].fixed = 1L & parinfo_new[0:iend].tied = ''
  parinfo_new[0:2].value = parinfo[0:2].value & parinfo_new[3:5].value = parinfo[0:2].value
  parinfo_new[1].value = alog(4862.68D) - 1d-3 & parinfo_new[4].value = alog(4862.68D) + 1d-3
  parinfo_new[0].limited = [1,0] & parinfo_new[0].limits = [0,0] ; limited to emission for the additional gaussian
  parinfo_new[1].limited = [1,1] & parinfo_new[1].limits = [alog(4862.68D) - 0.015D, alog(4862.68D) + 0.015D]
  parinfo_new[2].limited = [1,1]
  parinfo_new[2].limits = [0.0017D, max_sigma] & parinfo_new[5].limits = [0.0017D, max_sigma]
  line_fit_new = mpfitfun('manygauss', alog(wave_fit[indd]), flux_fit[indd], err_fit[indd] $
                     , functargs = functargs_new, parinfo = parinfo_new, perror = perror_new, yfit = yfit2_new, /quiet $
                     , nfev = nfev, niter = niter, status = status_new, bestnorm = chi2_new)
  if not keyword_set(silent) then splog, 'LINE-FITTING: MPFIT nfev=', nfev, ' niter=', niter, ' status=', status_new, ' chi2=', chi2_new
  if chi2_new gt 1d-6 and chi2_new lt chi2 and status_new ne 0 then begin ; replace the original fits with new fits
    parinfo = parinfo_new & functargs = functargs_new
    line_fit = line_fit_new & perror = perror_new & status = status_new & chi2 = chi2_new
    ngauss_broad_Hb = 2L
    para.ngauss_broad_Hb = ngauss_broad_Hb
    ; try not fixing the parameters for the narrow lines
    parinfo_new.value = line_fit & parinfo_new.fixed = 0L
    parinfo_new[10].tied = 'P(7) - alog(4862.68D) + alog(4960.30D)' & parinfo_new[11].tied = 'P(8)'
    parinfo_new[13].tied = 'P(7) - alog(4862.68D) + alog(5008.24D)' & parinfo_new[14].tied = 'P(8)'
    if keyword_set(double_OIII) then begin
      parinfo_new[19].tied = 'P(16) - alog(4960.30D) + alog(5008.24D)' & parinfo_new[20].tied = 'P(17)'
    endif
    line_fit_new = mpfitfun('manygauss', alog(wave_fit[indd]), flux_fit[indd], err_fit[indd] $
                  , functargs = functargs_new, parinfo = parinfo_new, perror = perror_new, yfit = yfit2_new, /quiet $
                  , nfev = nfev, niter = niter, status = status_new, bestnorm = chi2_new)
    if not keyword_set(silent) then splog, 'LINE-FITTING: MPFIT nfev=', nfev, ' niter=', niter, ' status=', status_new, ' chi2=', chi2_new
    if chi2_new gt 1d-6 and chi2_new lt chi2 and status_new ne 0 then begin
       parinfo = parinfo_new & functargs = functargs_new
       line_fit = line_fit_new & perror = perror_new & status = status_new & chi2 = chi2_new
    endif

    nline = 1L + n_elements(parinfo)/3L  ; add another gaussian
    functargs_new = {nline:nline}
    parinfo_new = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0],tied:''}, 3*nline)
    iend = n_elements(parinfo_new) - 1L
    parinfo_new[3:iend] = parinfo
    parinfo_new[3:iend].value = line_fit & parinfo_new[9:iend].fixed = 1L & parinfo_new[0:iend].tied = ''
    parinfo_new[0:2].value = parinfo[0:2].value & parinfo_new[3:5].value = parinfo[0:2].value & parinfo_new[6:8].value = parinfo[0:2].value
    parinfo_new[1].value = alog(4862.68D) - 1d-3 & parinfo_new[4].value = alog(4862.68D) & parinfo_new[7].value = alog(4862.68D) + 1d-3
    parinfo_new[0].limited = [1,0] & parinfo_new[0].limits = [0,0] ; limited to emission for the additional gaussian    
    parinfo_new[1].limited = [1,1] & parinfo_new[1].limits = [alog(4862.68D) - 0.015D, alog(4862.68D) + 0.015D]
    parinfo_new[2].limited = [1,1]
    parinfo_new[2].limits = [0.0017D, max_sigma] & parinfo_new[5].limits = [0.0017D, max_sigma] & parinfo_new[8].limits = [0.0017D, max_sigma]
    line_fit_new = mpfitfun('manygauss', alog(wave_fit[indd]), flux_fit[indd], err_fit[indd] $
                     , functargs = functargs_new, parinfo = parinfo_new, perror = perror_new, yfit = yfit2_new, /quiet $
                     , nfev = nfev, niter = niter, status = status_new, bestnorm = chi2_new)
    if not keyword_set(silent) then splog, 'LINE-FITTING: MPFIT nfev=', nfev, ' niter=', niter, ' status=', status_new, ' chi2=', chi2_new
    if chi2_new gt 1d-6 and chi2_new lt chi2 and status_new ne 0 then begin
       parinfo = parinfo_new & functargs = functargs_new
       line_fit = line_fit_new & perror = perror_new & status = status_new & chi2 = chi2_new
       ngauss_broad_Hb = 3L
       para.ngauss_broad_Hb = ngauss_broad_Hb
       ; try not fixing the parameters for the narrow lines
       parinfo_new.value = line_fit & parinfo_new.fixed = 0L
       parinfo_new[13].tied = 'P(10) - alog(4862.68D) + alog(4960.30D)' & parinfo_new[14].tied = 'P(11)'
       parinfo_new[16].tied = 'P(10) - alog(4862.68D) + alog(5008.24D)' & parinfo_new[17].tied = 'P(11)'
       if keyword_set(double_OIII) then begin
         parinfo_new[22].tied = 'P(19) - alog(4960.30D) + alog(5008.24D)' & parinfo_new[23].tied = 'P(20)'
       endif
       line_fit_new = mpfitfun('manygauss', alog(wave_fit[indd]), flux_fit[indd], err_fit[indd] $
                     , functargs = functargs_new, parinfo = parinfo_new, perror = perror_new, yfit = yfit2_new, /quiet $
                     , nfev = nfev, niter = niter, status = status_new, bestnorm = chi2_new)
       if not keyword_set(silent) then splog, 'LINE-FITTING: MPFIT nfev=', nfev, ' niter=', niter, ' status=', status_new, ' chi2=', chi2_new
       if chi2_new gt 1d-6 and chi2_new lt chi2 and status_new ne 0 then begin
          parinfo = parinfo_new & functargs = functargs_new
          line_fit = line_fit_new & perror = perror_new & status = status_new & chi2 = chi2_new
       endif

    endif
  endif
endif

;---------------------------------------------------------------------
temp = create_struct('line_fit', line_fit, 'line_fit_err', perror)
para = struct_addtags(para, temp)
dof = n_elements(indd) - n_elements(parinfo)
para.redchi2_2 = chi2/dof
para.status_2 = status
; set the line_npix, if the [4750-4950] region has less than 1/3 of good pixels then it should be regarded insufficent to fit the line
temp_ind = where(lam ge 4750.0 and lam le 4950.0 and err gt 1d-6, temp_npix)
para.line_npix = temp_npix
if temp_npix gt 0 then para.line_med_SN = median(flux[temp_ind]/err[temp_ind])


; set values for the broad-component of Hbeta
if not keyword_set(diet) then begin

if not keyword_set(gh) then begin
Hbeta1_linearea = dblarr(ngauss_broad_Hb) & Hbeta1_linearea_err = dblarr(ngauss_broad_Hb)
for i_gau=0L, ngauss_broad_Hb - 1L do begin
Hbeta1_linearea[i_gau] = line_fit[0+i_gau*3]*exp(line_fit[1+i_gau*3])
Hbeta1_linearea_err[i_gau] = perror[0+i_gau*3]*exp(line_fit[1+i_gau*3])
endfor
para.Hbeta1_linearea = total(Hbeta1_linearea,/double)
para.Hbeta1_linearea_err = sqrt(total(Hbeta1_linearea_err^2))

if total(Hbeta1_linearea,/double) gt 1d-6 then begin ; only measure the FWHM when the broad Hbeta is not zero-flux
if ngauss_broad_Hb eq 1L then begin ; if only one gaussian is used, it is easy to compute FWHM
wave_peak = line_fit[1]
para.Hbeta1_wave = exp(wave_peak)
para.Hbeta1_linesigma = line_fit[2]*c
para.Hbeta1_linesigma_err = perror[2]*c
para.Hbeta1_FWHM = line_fit[2]*2.3548D*c
para.Hbeta1_FWHM_err = perror[2]*2.3548D*c
endif else begin  ; if multiple gaussian is used, then find the non parametric FWHM
xgrid = 8.45532D + indgen(400L)*2d-4
ygrid = manygauss(xgrid, line_fit[0:3*ngauss_broad_Hb - 1L], nline=ngauss_broad_Hb)
; find the peak pixel
model_peak = max(ygrid, ind_peak)
wave_peak = xgrid[ind_peak]
para.Hbeta1_wave = exp(wave_peak)
ind_left = where(xgrid le wave_peak) & ind_right = where(xgrid ge wave_peak)
xgrid_left = xgrid[ind_left] & ygrid_left = ygrid[ind_left]
xgrid_right = xgrid[ind_right] & ygrid_right = ygrid[ind_right]
ind_sort = sort(ygrid_left)
half_left = interpol(xgrid_left[ind_sort], ygrid_left[ind_sort], 0.5*model_peak)
ind_sort = sort(ygrid_right)
half_right = interpol(xgrid_right[ind_sort],ygrid_right[ind_sort], 0.5*model_peak)
;print, half_left, half_right
para.Hbeta1_FWHM = (half_right - half_left)*c ; FWHM
endelse
; estimate lineew using the continuum level at the model peak (this is already the rest EW)
para.Hbeta1_lineew = para.Hbeta1_linearea/(conti_fit[0]*(exp(wave_peak)/3000.0)^conti_fit[1])
; estimate lineew err assuming it comes from the error in linearea
para.Hbeta1_lineew_err = para.Hbeta1_linearea_err/(conti_fit[0]*(exp(wave_peak)/3000.0)^conti_fit[1])
endif

inext = (ngauss_broad_Hb - 1L)*3
; set values for the narrow-component of Hbeta
para.Hbeta2_wave = exp(line_fit[4+inext])
para.Hbeta2_linearea = line_fit[3+inext]*exp(line_fit[4+inext])
para.Hbeta2_linearea_err = perror[3+inext]*exp(line_fit[4+inext])
para.Hbeta2_linesigma = line_fit[5+inext]*c
para.Hbeta2_linesigma_err = perror[5+inext]*c
para.Hbeta2_FWHM = line_fit[5+inext]*2.3548D*c
para.Hbeta2_FWHM_err = perror[5+inext]*2.3548D*c
; estimate lineew using the continuum level at the linecenter (this is already the rest EW)
para.Hbeta2_lineew = para.Hbeta2_linearea/(conti_fit[0]*(exp(line_fit[4+inext])/3000.0)^conti_fit[1])
; estimate lineew err assuming it comes from the error in linearea
para.Hbeta2_lineew_err = para.Hbeta2_linearea_err/(conti_fit[0]*(exp(line_fit[4+inext])/3000.0)^conti_fit[1])

; set values for the [OIII]4959
para.OIII4959_wave = exp(line_fit[7+inext])
if not keyword_set(double_oiii) then begin
para.OIII4959_linearea = line_fit[6+inext]*exp(line_fit[7+inext])
para.OIII4959_linearea_err = perror[6+inext]*exp(line_fit[7+inext])
endif else begin
;  if line_fit[12+inext] gt 1d-6 and perror[12+inext] gt 1d-6 then begin
    para.OIII4959_linearea = line_fit[6+inext]*exp(line_fit[7+inext]) + line_fit[12+inext]*exp(line_fit[13+inext])
    para.OIII4959_linearea_err = sqrt( (perror[6+inext]*exp(line_fit[7+inext]))^2 + (perror[12+inext]*exp(line_fit[13+inext]))^2)
;  endif
endelse
para.OIII4959_linesigma = line_fit[8+inext]*c
para.OIII4959_linesigma_err = perror[8+inext]*c
para.OIII4959_FWHM = line_fit[8+inext]*2.3548D*c
para.OIII4959_FWHM_err = perror[8+inext]*2.3548D*c
; estimate lineew using the continuum level at the linecenter (this is already the rest EW)
para.OIII4959_lineew = para.OIII4959_linearea/(conti_fit[0]*(exp(line_fit[7+inext])/3000.0)^conti_fit[1])
; estimate lineew err assuming it comes from the error in linearea
para.OIII4959_lineew_err = para.OIII4959_linearea_err/(conti_fit[0]*(exp(line_fit[7+inext])/3000.0)^conti_fit[1])

; set values for the [OIII]5007
para.OIII5007_wave = exp(line_fit[10+inext])
if not keyword_set(double_OIII) then begin
para.OIII5007_linearea = line_fit[9+inext]*exp(line_fit[10+inext])
para.OIII5007_linearea_err = perror[9+inext]*exp(line_fit[10+inext])
endif else begin
;  if line_fit[15+inext] gt 1d-6 and perror[15+inext] gt 1d-6 then begin
    para.OIII5007_linearea = line_fit[9+inext]*exp(line_fit[10+inext]) + line_fit[15+inext]*exp(line_fit[16+inext])
    para.OIII5007_linearea_err = sqrt( (perror[9+inext]*exp(line_fit[10+inext]))^2 +(perror[15+inext]*exp(line_fit[16+inext]))^2 )
;  endif
endelse
para.OIII5007_linesigma = line_fit[11+inext]*c
para.OIII5007_linesigma_err = perror[11+inext]*c
para.OIII5007_FWHM = line_fit[11+inext]*2.3548D*c
para.OIII5007_FWHM_err = perror[11+inext]*2.3548D*c
; estimate lineew using the continuum level at the linecenter (this is already the rest EW)
para.OIII5007_lineew = para.OIII5007_linearea/(conti_fit[0]*(exp(line_fit[10+inext])/3000.0)^conti_fit[1])
; estimate lineew err assuming it comes from the error in linearea
para.OIII5007_lineew_err = para.OIII5007_linearea_err/(conti_fit[0]*(exp(line_fit[10+inext])/3000.0)^conti_fit[1])
endif

endif

; now estimate the fitted line_flux in the wave_fit range ([4435,5535])
if not keyword_set(gh) then begin
  line_flux = manygauss(alog(wave_fit), line_fit)
endif else line_flux = manygh(alog(wave_fit), line_fit, nline=4L)

; add fitted line flux to the para structure
if not keyword_Set(diet) then begin
  temp = create_struct('line_flux', line_flux)
  para = struct_addtags(para, temp)
endif

fit_flux = conti_flux + Fe_flux + line_flux


if not keyword_set(noplot) then begin
if keyword_set(psplot) then begin
  ;set_plot, 'PS'
  ;device, filename = figfile, /color, xsize = 30, ysize = 20, /cmyk, /ENCAPSULATED
  ys_psopen, filename =figfile, /color, xsize=30,ysize=20
  thick = 5.
endif

yrange = [-1.0, max((smooth(flux,3))[ind])*1.05]
if not keyword_set(xrange) then xrange = [4400, 5600]
plot, lam, smooth(flux,1), xrange = xrange, yrange = yrange, /ystyle $
         , xtitle = textoidl('Rest Wavelength (\AA)'), ytitle = textoidl('Flux Density (10^{-17} erg s^{-1} cm^{-2} \AA^{-1})') $
         , charsize = 2.0, xthick = 3.0, ythick = 3.0, charthick = 3.0, _extra = extra
;oplot, lam, flux, psym=2, color=fsc_color('black'), symsize=0.6
         
if keyword_set(SDSS_name) then xyouts, 0.18, 0.85, /norm, SDSS_name, charsize = 2.5, charthick = 3.0
if keyword_set(print_para) then begin
  V1 = (line_fit[1]-line_fit[4])*c & FWHM1 = line_fit[2]*c*2.3548D
  V2 = (line_fit[7]-line_fit[4])*c & FWHM2 = line_fit[8]*c*2.3548D
  xyouts, 0.18, 0.8, textoidl('V_{off}=')+string(V1, format='(i0)') + ','+string(V2, format='(i0)'),/norm,charsize = 2.5, charthick = 3.0
  xyouts, 0.18, 0.75, textoidl('FWHM=')+string(FWHM1, format='(i0)') + ','+string(FWHM2, format='(i0)'),/norm,charsize = 2.5, charthick = 3.0
endif

if keyword_set(iterate) then xyouts, 0.18, 0.8, textoidl('N_{gauss,br}=')+string(ngauss_broad_Hb, format='(i0)'),/norm,charsize = 2.5, charthick = 3.0

oplot, lam, err, color = fsc_color('gray'), thick = thick
oplot, lam[ind], smooth(flux_fit,1)
;oplot, lam[ind], flux_fit, psym=2, color=fsc_color('black'), symsize=0.6
oplot, lam[ind], fit_flux, color = fsc_color('red'), thick = thick
oplot, lam[ind], conti_flux, color = fsc_color('orange'), thick = thick
oplot, lam[ind], Fe_flux, color = fsc_color('blue'), thick = thick

if not keyword_set(gh) then begin ; plot the multiple Gaussian components
for i_gau=0L, ngauss_broad_Hb - 1L do begin
  oplot, lam[ind], onegauss(alog(lam[ind]), line_fit[0+i_gau*3:2+i_gau*3]), color = fsc_color('green'), thick = thick
endfor
; plot the narrow Hbeta
oplot, lam[ind], onegauss(alog(lam[ind]), line_fit[3+inext:5+inext]), color = fsc_color('cyan'), thick = thick
if keyword_set(oiii_plot) then begin
  oplot, lam[ind], onegauss(alog(lam[ind]), line_fit[6+inext:8+inext]), color = fsc_color('cyan'), thick = thick
  oplot, lam[ind], onegauss(alog(lam[ind]), line_fit[9+inext:11+inext]), color = fsc_color('cyan'), thick = thick
  if not keyword_set(double_OIII) and keyword_set(fit_HeII) then begin
     oplot, lam[ind], onegauss(alog(lam[ind]), line_fit[12+inext:14+inext]), color = fsc_color('cyan'), thick = thick
     oplot, lam[ind], onegauss(alog(lam[ind]), line_fit[15+inext:17+inext]), color = fsc_color('cyan'), thick = thick
  endif
  if keyword_set(double_OIII) then begin
     oplot, lam[ind], onegauss(alog(lam[ind]), line_fit[12+inext:14+inext]), color = fsc_color('cyan'), thick = thick
     oplot, lam[ind], onegauss(alog(lam[ind]), line_fit[15+inext:17+inext]), color = fsc_color('cyan'), thick = thick
     if keyword_set(fit_HeII) then begin
       oplot, lam[ind], onegauss(alog(lam[ind]), line_fit[18+inext:20+inext]), color = fsc_color('cyan'), thick = thick
       oplot, lam[ind], onegauss(alog(lam[ind]), line_fit[21+inext:23+inext]), color = fsc_color('cyan'), thick = thick
     endif
  endif
endif
endif else begin ; plot the multiple Gaussian-Hermite components
  oplot, lam[ind], onegh(alog(lam[ind]), line_fit[0:4]), color=fsc_color('green'), thick = thick  ; broad Hbeta
  oplot, lam[ind], onegh(alog(lam[ind]), line_fit[5:9]), color=fsc_color('cyan'), thick = thick  ;  narrow Hbeta
  oplot, lam[ind], onegh(alog(lam[ind]), line_fit[10:14]), color=fsc_color('cyan'), thick = thick ; OIII4959
  oplot, lam[ind], onegh(alog(lam[ind]), line_fit[15:19]), color=fsc_color('cyan'), thick = thick ; OIII5007
endelse
oplot, lam[ind], line_flux, color = fsc_color('Magenta'), thick = thick

if keyword_set(psplot) then begin
  ys_psclose
  ; device, /close_file
  ; set_plot, 'X'
endif
endif ; noplot

endif
endif

endif else begin 
  if not keyword_set(silent) then splog, 'redchi2 > 10 for the continuum + Fe fits. Stop there.'
endelse

if keyword_set(fits) then $
mwrfits, para, fitsfile, /create, /silent

endif

; restore the original wavelengths
lam = lam*(1.0 + z)

return

end
