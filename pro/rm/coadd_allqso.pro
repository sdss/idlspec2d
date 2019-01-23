; coadd the spectra of all 849 qsos and test the difference between 
; the pipeline redshift (zfinal in target_fibermap) and systemic redshift (zsys)

pro coadd_allqso

    file = '/data3/yshen/work/sdssrm_sample_char/sample_char.fits'
    target = mrdfits(file,1)
    nobj = n_elements(target)

    ; first construct a struct containing plate, fiber, mjd, redshift, ra, dec
    str_zpip = {plate:0L, fiber:0L, mjd:0L, redshift:0.D, ra:0.d, dec:0.d}
    str_zpip = replicate(str_zpip, nobj)

    str_zpip.plate = 0L & str_zpip.fiber = lindgen(849) + 1L & str_zpip.mjd = 56837L
    str_zpip.redshift = target.zfinal & str_zpip.ra = target.ra & str_zpip.dec = target.dec

    str_zsys = str_zpip
    str_zsys.redshift = target.zsys
  

    fitsfile = '/data3/yshen/work/sdssrm_sample_char/med_coadd_zpip.fits'
    coadd_spec, bh=str_zpip, type = 1, /sdssrm, /boss, fitsfile = fitsfile, z_HW=0

    fitsfile = '/data3/yshen/work/sdssrm_sample_char/med_coadd_zsys.fits'
    coadd_spec, bh=str_zsys, type = 1, /sdssrm, /boss, fitsfile = fitsfile, z_HW=0


end

