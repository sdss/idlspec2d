; given RMID (scalar or array), print target information
; and additional info

pro print_target_info, rmid, line=line

red, omegalambda=0.7,omega0=0.3,h100=0.7

targetfile=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
target=mrdfits(targetfile,1)

; host decomposed fits
file='/data3/quasar/yshen/work/agn_host/decomp_final.fits'
decomp=mrdfits(file,1)

; qsofit
file='/data3/quasar/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/more_lines/qso_prop-0000-56837.fits'
qso=mrdfits(file,1)

nnn=n_elements(rmid)
if n_elements(line) eq 0 then line=replicate('hbeta',nnn)

print, 'RMID SDSSNAME z  morph  line f5100 logL5100_qso fhost FWHM logMSE_Hbeta'

fmt='(i3.3, A20, " ", f6.4, " ", a9, " ", a5, " ", f0.3, " ", f0.3, " ", 2(f6.3, " "), " ", f4.2, " ", i5, " ", i5, " ", f6.3, " ", f6.3)'
for i=0,nnn-1 do begin
  ind=where(decomp.fiber eq rmid[i]+1)
  if target[rmid[i]].objc_type eq 6 then mflag='point' else mflag='extended'

  if line[i] eq 'hbeta' then begin
     FWHM=decomp[ind].FWHM_HB & FWHM_err=decomp[ind].FWHM_HB_Err
  endif
  if line[i] eq 'mgii' then begin
     FWHM=(qso[rmid[i]].MgII)[1]
     FWHM_err=(qso[rmid[i]].MgII_err)[1]
  endif

  ; convert L5100 to flux
  const = 1./(4.*!PI)/(dluminosity(decomp[ind].z,/cm))^2/(1. + decomp[ind].z)*1d17/5100.D
  f5100 = 10.D^decomp[ind].LOGL5100_QSO*const
  f5100_err = f5100*alog(10.D)*decomp[ind].LOGL5100_QSO_err

  print, rmid[i], decomp[ind].sdss_name,decomp[ind].z,mflag,line[i], f5100, f5100_err, $
   decomp[ind].LOGL5100_QSO,decomp[ind].LOGL5100_QSO_err, decomp[ind].f_H_5100, $
   round(FWHM), round(FWHM_Err), decomp[ind].LOGMBH_VP06,decomp[ind].LOGMBH_VP06_err, format=fmt

endfor

end
