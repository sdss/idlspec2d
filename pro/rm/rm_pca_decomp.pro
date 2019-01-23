;+
; NAME:
;   rm_pca_decomp
;
; PURPOSE:
;   decompose quasar spectra to host+AGN using PCA templates
;
; INPUTS:
;   wave_in  --  input observed wavelength array [NPIX]
;   flux_in  --  input total flux array [NPIX]
;   err_in   --  input error array [NPIX]
;
; OPTIONAL INPUTS:
;   z        --  redshift
;   ra,dec   --  J2000 coordinates; required if /deredden
;   qso_npca --  number of QSO eigenspectra used; default 10
;   gal_npca --  number of galaxy eigenspectra used; default 5
;
; OUTPUTS:  
;   result   --  output structure
;-------------------------

pro rm_pca_decomp,wave_in,flux_in,err_in,z=z,ra=ra,dec=dec,deredden=deredden, $
     pca_lib=pca_lib,qso_npca=qso_npca,gal_npca=gal_npca, $
     result=result


  ; common block of pca templates
  common pca_data,wave_qso,pca_qso,wave_gal,pca_gal
  if ~keyword_set(qso_npca) then qso_npca=10
  if ~keyword_set(gal_npca) then gal_npca=5

  ; get default pca templates
  if n_elements(wave_qso) eq 0 then begin
     pca_file='/data1/quasar/yshen/Project/pca_decomp/Yip_pca_templates/qso_eigenspec_Yip2004_CZBIN1.fits'
     ttt=mrdfits(pca_file,1)
     wave_qso=ttt.wave
     pca_qso=ttt.pca
  endif
  if n_elements(wave_gal) eq 0 then begin
     pca_file='/data1/quasar/yshen/Project/pca_decomp/Yip_pca_templates/gal_eigenspec_Yip2004.fits'
     ttt=mrdfits(pca_file,1)
     wave_gal=ttt.wave
     pca_gal=ttt.pca
  endif
  
  ; process the input spectrum;
  if n_elements(ra) eq 0 then ra=0. 
  if n_elements(dec) eq 0 then dec=0.
  if keyword_set(deredden) then dereddening_spec,wave_in,flux_in,err=err_in $
         , ra=ra, dec=dec, dered_flux=flux, dered_err=err else begin
     flux=flux_in & err=err_in
  endelse
  npix=n_elements(flux)
  ivar=err*0
  ind=where(err gt 1d-6)
  if ind[0] ne -1 then ivar[ind]=1./(err[ind])^2

  if n_elements(z) eq 0 then z=0.
  wave=wave_in/(1.D + z)

  ; setup output struct
  result={wave:wave,flux_in:flux_in,err_in:err_in,flux:flux,err:err,ivar:ivar, z:z, $
          ra:ra, dec:dec, npca_qso:qso_npca, npca_gal:gal_npca, $
          acoeff_qso:dblarr(qso_npca), acoeff_gal:dblarr(gal_npca), $
          flux_recon:dblarr(npix), qso_recon:dblarr(npix), gal_recon:dblarr(npix), chi2:0.D, dof:0L, $
          f_H:0., f_H_5100:0.}

  ;++++++++++++++++++++++
  ; remap qso and gal eigenspec to the available wavelength
  pca_use_qso = dblarr(npix, qso_npca)
  for i=0L, qso_npca - 1 do begin
    combine1fiber,alog10(wave_qso),pca_qso[*,i],newloglam=alog10(wave), $
      newflux=newflux
    pca_use_qso[*,i] = newflux
  endfor
  pca_use_gal = dblarr(npix, gal_npca)
  for i=0L, gal_npca - 1 do begin
    combine1fiber,alog10(wave_gal),pca_gal[*,i],newloglam=alog10(wave), $
      newflux=newflux
    pca_use_gal[*,i] = newflux
  endfor
  pca_use = [ [pca_use_qso], [pca_use_gal] ]
  ;-------------------------------

  ;message, 'stop and diag'

  ; fit the spectrum in a chi^2 sense
  ; make sure that the input wavelength does not exceed the range in the eigenspec
  ; QSO PCA template in CZBIN1 has maxwave=7500A
  maxwave=min([max(wave_qso),max(wave_gal)]) & minwave=max([min(wave_qso),min(wave_gal)])
  ind=where(wave le maxwave and wave ge minwave)
  result=struct_addtags(result, {fit_ind:ind})

  chi2 = computechi2( flux[ind], sqrt(ivar[ind]), pca_use[ind,*], $
    acoeff=acoeff, dof=dof, yfit=yfit)

  result.flux_recon[ind]=yfit
  result.acoeff_qso=acoeff[0:qso_npca-1]
  result.acoeff_gal=acoeff[qso_npca:*]
  result.chi2=chi2 & result.dof=dof
  result.qso_recon=result.acoeff_qso ## pca_use_qso
  result.gal_recon=result.acoeff_gal ## pca_use_gal

  ; compute host fraction
  ind=where(wave gt 4160. and wave lt 4210. and err gt 1d-6, nn)
  if nn gt 10 then begin
     f_gal = int_tabulated(wave[ind], result.gal_recon[ind], /double)
     f_agn = int_tabulated(wave[ind], result.qso_recon[ind], /double)
     if f_agn + f_gal gt 0 then result.f_H=f_gal/(f_agn + f_gal)
  endif
  ind=where(wave gt 5080. and wave lt 5130. and err gt 1d-6, nn)
  if nn gt 10 then begin
     f_gal = int_tabulated(wave[ind], result.gal_recon[ind], /double)
     f_agn = int_tabulated(wave[ind], result.qso_recon[ind], /double)
     if f_agn + f_gal gt 0 then result.f_H_5100=f_gal/(f_agn + f_gal)
  endif

end











