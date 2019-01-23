; Output an ASCII table of the qsofit results

pro rm_output_qsofit,linename=linename

file='/data3/quasar/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/more_lines/qso_prop-0000-56837.fits'

qso=mrdfits(file,1)
tags=strtrim(tag_names(qso))

if ~keyword_set(linename) then $
   linename=['Halpha_br', 'Hbeta_br', 'HeII4687','OIII5007', $
             'MgII', 'CIII', 'CIV','SiIV_OIV', 'Lya']


if ~keyword_set(outfile) then outfile= $
  '/data3/quasar/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/more_lines/ascii_qsofit.dat'
fmt='(i3.3, " ",f6.4, 9(" ", f9.4," ",f9.4, " ", i5, " ", i5)  )'
openw,lun,outfile,/get_lun
printf, lun, '# RMID z ', linename
printf, lun, '# line [wave, wave_err, FWHM, FWHM_err]'

nline=n_elements(linename)
for iobj=0L, n_elements(qso) - 1 do begin
  
  output=dblarr(4*nline)
  for i=0L, nline - 1 do begin
  
     ind_tag = where(tags eq strupcase(linename[i]) )
     ind_err = where(tags eq strupcase(linename[i]+'_ERR') )

     output[i*4:i*4+3] = [ (qso[iobj].(ind_tag))[0], (qso[iobj].(ind_err))[0], $
                           (qso[iobj].(ind_tag))[1], (qso[iobj].(ind_err))[1] ]
  
  endfor
  printf, lun, format=fmt, qso[iobj].rm_ID, qso[iobj].z, output

endfor


close,lun
free_lun, lun

end
