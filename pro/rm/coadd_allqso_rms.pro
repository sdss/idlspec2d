; coadd the rms spectra of all RM targets

pro coadd_allqso_rms, topdir=topdir

    if ~keyword_set(topdir) then topdir = '/data3/yshen/ftp/sdssrm/collab/prepspec/2014a/'

    file = topdir + 'var_img.fits'
    var_img = mrdfits(file,1)
    bh = {plate:0, fiber:0, mjd:0, redshift:0.D, RA:0.D, dec:0.D}
    nobj = n_elements(var_img.z)
    bh = replicate(bh, nobj)
    wave = var_img.wave
    npix = n_elements(wave)
    lam_in = dblarr(npix, nobj)
    ivar_in = dblarr(npix, nobj)
    ivar_in[*,*] = 1.D ; set all pixel to be OK in the stacking (although we could set using ivar from the original 56837 coadd)
    for i=0,nobj - 1 do lam_in[*,i] = wave

    ; get zsys
    file = '/data3/yshen/work/sdssrm_sample_char/sample_char_refined.fits'
    tt = mrdfits(file, 1)
    rest_z = tt.zsys ; this should be in the same order (RMID) as bh


    ; use the quick_coadd instead of coadd_spec (used for Vandenberk2001 composite) because
    ; the following doesn't work
    ;fitsfile = '/data3/yshen/work/sdssrm_sample_char/RMSx_medcoadd_zsys.fits'
    ; coadd_spec, bh=bh, type=1,deredden=0, fitsfile=fitsfile, rest_z=rest_z, z_HW=0, $
    ;    lam_in = lam_in, flux_in = flux_in, ivar_in = ivar_in


    ; coadd RMSc, the model RMS (error accounted for, pt corrected)
    flux_in = var_img.IMG_rmsc_norm
    quick_coadd, lam_in, flux_in, rest_z, result=result_rmsc
    fitsfile = topdir + 'coadd/RMSc_medcoadd_zsys.fits'
    mwrfits, result_rmsc, fitsfile, /create

    ; coadd the RMSx, maximum-likelihood estimate (error corrected, but pt corrected?)
    flux_in = var_img.IMG_rms_norm
    quick_coadd, lam_in, flux_in, rest_z, result=result_rmsx
    fitsfile = topdir + 'coadd/RMSx_medcoadd_zsys.fits'
    mwrfits, result_rmsx, fitsfile, /create

    ; get the median coadd of the avg spectra
    file = '/data3/yshen/work/sdssrm_sample_char/med_coadd_zsys.fits'
    medcoadd=mrdfits(file,1)

    ; make a plot
    figfile = topdir + 'coadd/avg_rms_coadd_spec.ps'
    begplot, name=figfile, /color, /landscape

    plot, result_rmsc.wave, result_rmsc.coadd_flux + 1.5, xrange=[500, 9000],/xsty, pos=[0.1, 0.1, 0.98, 0.98], $
      xtitle='Rest Wavelength', ytitle = 'Flux Density', yrange=[1.5,10], thick=4,/ylog,/ysty, /nodata
    oplot, result_rmsx.wave, result_rmsx.coadd_flux + 1.5, color=cgcolor('dark grey'),thick=4
    oplot, medcoadd.wave, medcoadd.coadd_flux/median(medcoadd.coadd_flux) + 1.5, color=cgcolor('red'), thick=4
    oplot, result_rmsc.wave, result_rmsc.coadd_flux + 1.5, thick=4

    label_qsoline, wave=result_rmsc.wave,flux=result_rmsc.coadd_flux,charsize=0.6,minREW=0.5, ypos=0.15, floaty=0
    xyouts, 0.5, 0.9, 'Avg spectra', /norm, color=cgcolor('red')
    xyouts, 0.7, 0.9, 'RMSc spectra', /norm, color=cgcolor('black')
    xyouts, 0.5, 0.85, 'RMSx spectra', /norm, color=cgcolor('dark grey')

    endplot
    cgfixps, figfile

end
   
