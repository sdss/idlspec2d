; add more columns and modify the decomp.fits of the AGN/host spectral decomposition

pro add_decomp_fits

file='/data3/quasar/yshen/work/agn_host/decomp.fits'
result=mrdfits(file,1)
nnn=n_elements(result)
struct={medsn_all:0.D, logMBH_vp06:0.D, logMBH_vp06_err:-1.D, $
    logMBH_F14:0.D, logMBH_F14_err:-1.D, sigma_OK_Shen:0L, sigma_OK_final:0L, $
    logL5100_tot:0.D, logL5100_tot_err:-1.D, $
    sigma_greene_3780_4000:0.D, sigma_greene_3780_4000_err:-1.D, $
    sigma_greene_3780_4100:0.D, sigma_greene_3780_4100_err:-1.D }
struct=replicate(struct,nnn)

result=struct_addtags(result, struct)

;; replace the Greene's sigma with the correct masks ??
;file='/data3/quasar/yshen/work/agn_host/greene/final/valdes_4000-5350_correct_mask.dat'
;readcol,file,format='a,d,d',sdss_name, vdisp, vdisp_err
;nn=n_elements(sdss_name)
;for i=0L,nn-1 do begin
;  ind=where(result.sdss_name eq sdss_name[i])
;  result[ind].sigma_greene = vdisp[i]
;  result[ind].sigma_greene_err = vdisp_err[i]
;endfor


; add H+K dispersions
file='/data3/quasar/yshen/work/agn_host/greene/final/valdes_3780-4000.dat'
readcol,file,format='a,d,d',sdss_name1,vdisp1,vdisp1_err
file='/data3/quasar/yshen/work/agn_host/greene/final/valdes_3780-4100.dat'
readcol,file,format='a,d,d',sdss_name2,vdisp2,vdisp2_err
nn2=n_elements(sdss_name1)
for i=0L,nn2-1 do begin
  ind=where(result.sdss_name eq sdss_name1[i])
  result[ind].sigma_greene_3780_4000 = vdisp1[i]
  result[ind].sigma_greene_3780_4000_err = vdisp1_err[i]
  ind=where(result.sdss_name eq sdss_name2[i])
  result[ind].sigma_greene_3780_4100 = vdisp2[i]
  result[ind].sigma_greene_3780_4100_err = vdisp2_err[i]
endfor

for i=0L, nnn-1 do begin
   sn=result[i].flux*sqrt(result[i].ivar)
   result[i].medsn_all = median(sn)
endfor

file = '/data3/quasar/yshen/work/agn_host/checksheet'
readcol,file,format='a,a',sdss_name, shen_OK, skip=1L
nobj=n_elements(shen_OK)
for i=0, nobj-1 do begin
  ind=where(result.sdss_name eq sdss_name[i])
  result[ind].sigma_ok_shen=1L
endfor

; fix the two objects with ID=121 and 153 for which 
; the Greene sigma values seem more correct upon visual inspection
result[121].sigma_shen = result[121].sigma_greene
result[121].sigma_shen_err = result[121].sigma_greene_err
result[153].sigma_shen = result[153].sigma_greene
result[153].sigma_shen_err = result[153].sigma_greene_err

; now assign the final sigma_OK_final flag
ind=where( result.sigma_OK_shen eq 1L   and $
           result.sigma_greene_err gt 0 and $  
           abs(result.sigma_shen - result.sigma_greene) lt $
            2.*sqrt(result.sigma_shen_err^2+result.sigma_greene_err^2) and $
           result.sigma_shen gt 3.*result.sigma_shen_err and $
           result.sigma_shen_err gt 0  )
result[ind].sigma_OK_final = 1L

; now compute BH mass
logL5100_q=result.logL5100_qso
logL5100=result.logL5100_QSO
f_H=result.f_H & f_H_5100=result.f_H_5100
; recover the total logL5100 from the qso fit, note that the qsofit only fit the
; gal-subtracted quasar spectrum if f_H>0.1
ind=where(f_H gt 0.1 and logL5100 gt 0)
logL5100[ind]= alog10(10.D^logL5100[ind]/(1. - f_H_5100[ind]))
; note that the logL5100_qso in result is already the total logL5100 if f_H>0.1
result.logL5100_tot=logL5100 & result.logL5100_tot_err=result.logL5100_qso_err
; The qsofit only compute pure qso logL5100 if f_H>0.1, but here we should 
; get the qso logL5100 as long as f_H>0.05, in other words, we will modified the
; logL5100_qso values
ind=where(f_H gt 0.05 and f_H le 0.1)
result[ind].logL5100_qso=alog10((1. - f_H_5100[ind])*10.D^result[ind].logL5100_tot)
; now clean all failed decomposition, i.e., f_H<0.05
ind=where(result.f_H le 0.05)
result[ind].f_h=0. & result[ind].f_h_5100=0. & result[ind].medsn_sigma=0.

; fix one quasar (RMID=764, fiber=765 in 0000-56837), indx=177 in decomp.fits
; I only have fitted coadded spectra 56783 (instead of 56837) but for logL5100 the difference
; is negligible
file1='/data3/quasar/yshen/spectro/bossredux/v5_7_1/0000/qsofit/qso_prop-0000-56783.fits'
qso=mrdfits(file1,1)
result[177].logL5100_tot=qso[764].logL5100
result[177].logL5100_tot_err=qso[764].logL5100_err
result[177].logL5100_qso_err=qso[764].logL5100_err
result[177].logL5100_qso=alog10(10.D^result[177].logL5100_tot*(1. - result[177].f_H_5100))
;------------------------

logL5100_q=result.logL5100_qso
logL5100=result.logL5100_tot
fwhm=result.fwhm_hb & fwhm_err=result.fwhm_hb_err
ind=where(logL5100_q gt 0)
result[ind].logmbh_vp06=0.91 + 0.5*alog10(10.D^logL5100[ind]/1d44) + 2.*alog10(fwhm[ind])
result[ind].logmbh_vp06_err=sqrt( (0.5*result[ind].logL5100_tot_err)^2 + (2.*fwhm_err[ind]/fwhm[ind]/alog(10.D) )^2  )
result[ind].logmbh_f14=3.602 + 0.504*alog10(10.D^result[ind].logL5100_qso/1d44) + 1.2*alog10(fwhm[ind])
result[ind].logmbh_f14_err=sqrt( (0.504*result[ind].logL5100_qso_err)^2 + (1.2*fwhm_err[ind]/fwhm[ind]/alog(10.D))^2 )

outfile='/data3/quasar/yshen/work/agn_host/decomp_final.fits'
mwrfits, result, outfile, /create

end
