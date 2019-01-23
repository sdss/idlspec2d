; output a fits file of the host properties from
; Shen++2015 for Yoshiki Matsuoka

; adapt this to output additional columns

pro output_host_decomp


file='/data3/quasar/yshen/work/agn_host/decomp_final.fits'

column=['SDSS_NAME','RA','DEC','z','F_H_5100','fiber','sigma_shen','sigma_shen_err', $
  'sigma_ok_shen', 'SIGMA_GREENE_3780_4100', 'SIGMA_GREENE_3780_4100_err', $
   'logL5100_tot', 'logL5100_tot_err', 'FWHM_hb', 'FWHM_hb_err', $
   'LOGMBH_VP06', 'LOGMBH_VP06_ERR']
result=mrdfits(file,1,column=column)

column1=['SDSS_NAME','RA','DEC','z','F_H_5100', $
   'logL5100_tot', 'logL5100_tot_err', 'FWHM_hb', 'FWHM_hb_err', $
   'LOGMBH_VP06', 'LOGMBH_VP06_ERR']
output=mrdfits(file,1,column=column1)

nnn=n_elements(output)
struct=replicate({rmid:0L,sigma_4125_5350:0.D,sigma_4125_5350_err:-1.D, sigma_3780_4100:0.D, $
  sigma_3780_4100_err:-1.D},nnn)
output=struct_addtags(output, struct)
output.rmid=result.fiber - 1
ind=where(result.sigma_ok_shen eq 1)
output[ind].sigma_4125_5350=result[ind].sigma_shen
output[ind].sigma_4125_5350_err=result[ind].sigma_shen_err
output[ind].sigma_3780_4100=result[ind].SIGMA_GREENE_3780_4100
output[ind].sigma_3780_4100_err=result[ind].SIGMA_GREENE_3780_4100_err

outfile='/data3/quasar/yshen/work/agn_host/Yoshiki/pca_decomp.fits'
mwrfits, output, outfile, /create

end
