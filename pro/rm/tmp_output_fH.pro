; output f_H5100 for the 212 z<1.09 RM quasars

pro tmp_output_fH

file='/data3/quasar/yshen/work/agn_host/decomp_final.fits'
result=mrdfits(file,1)

fmt='(i3.3, " ", f0.4)'

outfile='/data3/quasar/yshen/work/agn_host/tables/f_H_5100'
openw, lun, outfile, /get_lun
printf, lun, '# RM_ID  f_H_5100=fgal/(fgal+fagn) at 5100A'

for i=0L, 211 do printf, lun, result[i].fiber - 1, result[i].f_H_5100,format=fmt

close, lun
free_lun, lun

end

; output the decomposed spectrum
pro tmp_output_decomp_spec

file='/data3/quasar/yshen/work/agn_host/decomp_final.fits'
result=mrdfits(file,1)
; keep only those with successful decomposition
ind=where(result.f_H gt 0, nnn)
result=result[ind]

outdir='/data3/quasar/yshen/work/agn_host/tables/host_spec/'
for i=0L, nnn-1 do begin
outfile=outdir+'RMID_'+string(result[i].fiber-1, format='(i3.3)') + '_decomp_spec'
openw, lun, outfile, /get_lun
printf, lun, '# RMID='+string(result[i].fiber-1, format='(i3.3)') + $
     '  z='+string(result[i].z,format='(f0.4)')
printf, lun, '# Dereddened spectra'
printf, lun, '# f_H='+string(result[i].f_H,format='(f0.3)') + $
   ' f_H_5100='+string(result[i].f_H_5100,format='(f0.3)')
printf, lun, '# Obs Wave   f_tot   f_host=f_tot-f_qso_recon  err  [1d-17 erg/s/cm2/A]'
wave=result[i].wave & flux=result[i].flux & flux_gal=result[i].flux_gal
err=result[i].err
fmt='(f10.3, e12.4, e12.4, e12.4)'
npix=n_elements(wave)
for jj=0L, npix-1 do printf, lun, wave[jj], flux[jj], flux_gal[jj], err[jj],format=fmt

close, lun
free_lun, lun

endfor


end
