; test the dependence of CIV FWHM on SNR, using 
; the coadded fits and the individual-epoch fits

pro test_fwhm_sn

file='/data3/yshen/work/lineshifts/lineshift.fits'
result0=mrdfits(file,1)
fwhm0=(result0.CIV_br)[1,*]

file='/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
target=mrdfits(file,1)
plate=(target.plate)[0:31,0]
fiberall=(target.fiberid)[0:31,*]
mjd=(target.mjd)[0:31,0]

for i=0, 31L do begin

  pstr=string(plate[i],format='(i4.4)')
  mstr=string(mjd[i],format='(i5.5)')
  fiber=reform(fiberall[i,*])
  file='/data3/yshen/ftp/sdssrm/collab/bossredux/v5_7_1/'+pstr+'/qsofit/qso_prop-'+pstr+'-'+mstr+'.fits'
  result1=mrdfits(file,1) 
  fwhm1=(result1.civ)[1,*]

  ; CIV measureable 
  ind=where(fwhm0 gt 0 and fwhm1 gt 0, nobj)
  ; get the corresponding median SNR
  rm_readspec, plate[i], fiber[ind], mjd=mjd[i], flux=flux,invvar=ivar
  snr_all = flux*sqrt(ivar)
  snr=median(snr_all, dim=1)

  if n_elements(ratio_arr) eq 0 then ratio_arr = (fwhm1[ind] - fwhm0[ind])/fwhm0[ind] else $
      ratio_arr=[ratio_arr, (fwhm1[ind] - fwhm0[ind])/fwhm0[ind]]
  if n_elements(snr_arr) eq 0 then snr_arr = snr else snr_arr=[snr_arr, snr]

  splog, 'Finished epoch:', i
endfor

output={dfwhm_frac:ratio_arr, med_SNR:snr_arr}

outfile='/data3/yshen/work/lineshifts/fwhm_coadd_se_test.fits'
mwrfits, output, outfile, /create

end

pro plot_test_fwhm_coadd_se

file='/data3/yshen/work/lineshifts/fwhm_coadd_se_test.fits'
result=mrdfits(file,1)
figfile='/data3/yshen/work/lineshifts/fwhm_sntest_coadd_and_SE.eps'
begplot, name=figfile, /landscape, /color
dfwhm_frac=result.dfwhm_frac & snr=result.med_snr

bin=0.05 & xrange=[-0.8, 0.4]
title='coadd versus single-epoch measurements'
plothist,dfwhm_frac, bin=bin, xtitle=textoidl('\DeltaFWHM/FWHM'),ytitle='N',xrange=xrange, /xsty, title=title
oplot, [median(dfwhm_frac), median(dfwhm_frac) ], [0, 3d3]
ind=where(snr gt 10)
plothist, dfwhm_frac[ind], bin=bin, /over, color=cgcolor('brown')
oplot, [median(dfwhm_frac[ind]), median(dfwhm_frac[ind]) ], [0, 3d3],color=cgcolor('brown')
ind=where(snr lt 10)
plothist, dfwhm_frac[ind], bin=bin, /over, color=cgcolor('cyan')
oplot, [median(dfwhm_frac[ind]), median(dfwhm_frac[ind]) ], [0, 3d3],color=cgcolor('cyan')
ind=where(snr lt 5)
plothist, dfwhm_frac[ind], bin=bin, /over, color=cgcolor('green')
oplot, [median(dfwhm_frac[ind]), median(dfwhm_frac[ind]) ], [0, 3d3],color=cgcolor('green')
colors=cgcolor(['black', 'brown', 'cyan', 'green'])
items=['All', 'S/N>10', 'S/N<10', 'S/N<5']
legend, pos=[0.2, 0.9], box=0, items, color=colors, textcolor=colors,/norm,line=[0,0,0,0]

endplot
cgfixps, figfile


end
