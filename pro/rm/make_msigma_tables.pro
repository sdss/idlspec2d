; make the fits data table for the m-sigma paper

pro make_msigma_tables

file='/data3/quasar/yshen/work/agn_host/decomp_final.fits'
result=mrdfits(file,1)


; make External Database 1 for the Science paper
;tags=['SDSS_name', 'RA', 'DEC', 'Z', 'WAVE', 'FLUX', 'ERR', 'Flux_gal', $
;    'MEDSN_sigma', 'MEDSN_all', 'flux_recon', 'F_H', 'F_H_5100', 'sigma_SHEN', $
;    'sigma_shen_err', 'sigma_err_' ]
;ind=where(result.sigma_OK_shen eq 1)
;result=result[ind]

; make table 1 for the ApJ version
tags=['SDSS_name', 'Fiber', 'RA', 'DEC', 'Z', 'WAVE', 'FLUX', 'ERR', 'Flux_gal', $
    'MEDSN_sigma', 'MEDSN_all', 'flux_recon', 'F_H', 'F_H_5100', 'sigma_SHEN', $
    'sigma_shen_err', 'sigma_err_', 'SIGMA_GREENE_3780_4100', 'SIGMA_GREENE_3780_4100_ERR']


output={RMID:0L,sdss_name:'', ra:0.D, dec:0.D, z:0., wave:dblarr(4649), flux:dblarr(4649), $
   err:dblarr(4649), flux_gal:dblarr(4649), flux_qso:dblarr(4649), medsn_tot:0., $
   medsn_gal:0., f_h:0., f_H_5100:0., sigma:0.D, sigma_err:-1.D, sigma_err_warning:0L, $
   sigma_HK:0.D, sigma_HK_err:-1.D, $
   logL5100_tot:0.D, logL5100_tot_err:-1.D, logL5100_qso:0.D, logL5100_qso_err:-1.D, $
   fwhm_HB:0.D, fwhm_hb_err:-1.D, logmbh_vp06:0.D, logmbh_vp06_err:-1.D }
output=replicate(output, n_elements(result))

output.rmid=result.fiber - 1
output.sdss_name=result.sdss_name & output.ra=result.ra & output.dec=result.dec
output.z=result.z & output.wave=result.wave & output.flux=result.flux
output.err=result.err 

; only populate the decomposed spectra and sigma_HK if f_H>0
ind=where(result.f_h gt 0,nnn)
output[ind].flux_gal=result[ind].flux_gal
output[ind].flux_qso = result[ind].flux - (result[ind].flux_recon - (result[ind].flux-result[ind].flux_gal))
output[ind].sigma_HK=result[ind].SIGMA_GREENE_3780_4100
output[ind].sigma_HK_err=result[ind].SIGMA_GREENE_3780_4100_err
for i=0L, nnn-1 do begin
  ; set the pixels at which flux_recon=0 to have zero flux_qso
  flux_qso=output[ind[i]].flux_qso
  flux_recon=result[ind[i]].flux_recon
  indd=where(flux_recon eq 0)
  if indd[0] ne -1 then flux_qso[indd]=0.
  output[ind[i]].flux_qso = flux_qso
endfor
output.medsn_tot=result.medsn_all & output.medsn_gal=result.medsn_sigma
output.f_h=result.f_h & output.f_H_5100=result.f_H_5100
; fix the 2 objects with bad sigma_HK measurements
ind=where(output.sigma_hk_err eq 0)
output[ind].sigma_hk_err = -1. & output[ind].sigma_hk = 0.

; only populate sigma with good measurements
ind=where(result.sigma_ok_shen eq 1)
output[ind].sigma=result[ind].sigma_shen & output[ind].sigma_err=result[ind].sigma_shen_err
ind=where(result.sigma_ok_final eq 0 and result.sigma_OK_shen eq 1)
output[ind].sigma_err_warning=1L
output.logL5100_tot=result.logL5100_tot & output.logL5100_tot_err=result.logL5100_tot_err
output.logL5100_qso=result.logL5100_qso & output.logL5100_qso_err=result.logL5100_qso_err
output.fwhm_hb=result.fwhm_hb & output.fwhm_hb_err=result.fwhm_hb_err
output.logmbh_vp06=result.logmbh_vp06 & output.logmbh_vp06_err=result.logmbh_vp06_err

; fix the one object with no 5100 coverage
ind=where(output.logL5100_tot eq 0)
output[ind].logL5100_tot_err = -1. & output[ind].logL5100_qso_err=-1.
ind=where(output.f_h gt 0 and output.f_h_5100 eq 0)
output[ind].logL5100_tot=0. & output[ind].logL5100_qso=0.
output[ind].logL5100_tot_err=-1. & output[ind].logL5100_qso_err=-1.
output[ind].logmbh_vp06=0. & output[ind].logmbh_vp06_err=-1.


;outfile='/data3/quasar/yshen/work/agn_host/tables/ED_S1.fits'
outfile='/data3/quasar/yshen/work/agn_host/tables/table1_apj.fits'
mwrfits, output, outfile, /create

end
