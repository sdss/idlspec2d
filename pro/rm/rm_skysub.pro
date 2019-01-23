; NAME:
;   rm_skysub
; PURPOSE:
;   Use the Wild&Hewett algorithm to subtract sky residual
;

pro rm_skysub,plate,fiber,mjd,z=z,qso=qso,star=star,_Extra=extra, $
     wave=wave, flux=flux, err=err, $
     newflux=newflux, newerr=newerr
    

  calibdir='recalib/'
  if ~keyword_set(flux) or ~keyword_set(err) then $
    rm_readspec,plate,fiber,mjd=mjd,calibdir=calibdir,loglam=loglam,wave=wave,flerr=err,flux=flux

  ; set up eigenspec common block
  common subtractOH1
  if n_elements(espec) eq 0 then begin
    ;file='/data3/quasar/yshen/ftp/bossredux/v5_7_1/wh_skysub/pca.fits'
    file=getenv('IDLRM_DIR')+'/template/sky_residual_pca.fits'

    pca=mrdfits(file,1)
    lambda=pca.loglam & espec=pca.eigenspec & pix_sky=pca.ind_sky & pix_nosky=pca.ind_nonsky
    ; cull out a red-wavelength region
    ind_tmp=where(10.^lambda gt 5000., ncomplement=ncomp)
    lambda=lambda[ind_tmp]
    pix_sky=pix_sky - ncomp & pix_nosky=pix_nosky - ncomp
    ; remove the nosky pixels below 5000 A
    ind_tmp=where(pix_nosky ge 0)
    pix_nosky=pix_nosky[ind_tmp]
  endif

  ; setup weights_in
  common plate_weights, sarr, weights_in
  if n_elements(sarr) eq 0 then begin
    gen_err_rescaling, plate, mjd, narr=narr, carr=carr, sarr=sarr
    weights_in = 1./(narr*sarr)
  endif

  newerr=err*sarr

  ; note that the sky residual eigenvectors were created on the set of
  ; sky spectra unnormalized by the scaled errors, which is likely different
  ; from Wild&Hewett in which the eigenvectors were created on the normalized
  ; sky spectra. However, the two approaches give negligible difference, so
  ; I am sticking to my approach

  newflux=subtractOH(flux,err,loglam,z,dummyplate, $
    rms=rms,nrecon=nrecon,weights=weights_in, qso=qso,star=star, _Extra=extra)

end
