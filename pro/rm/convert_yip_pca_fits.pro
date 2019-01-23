; convert the PCA templates by Yip et al. (2004) to fits files

pro convert_yip_pca_fits

  outdir='/data1/quasar/yshen/Project/pca_decomp/Yip_pca_templates/'

  ; galaxy PCA
  pca_dir=outdir+'GALeSpecPUBLIC/'
  result={wave:dblarr(3839L), pca:dblarr(3839L, 10L)}
  pca=dblarr(3839L, 10L)
  for i=0L, 9 do begin
    file=pca_dir+'galaxyKL_eigSpec_' + string(i+1,format='(i0)') + '.dat'
    readcol,file,format='d,d',wave, flux,/silent
    pca[*,i] = flux
  endfor
  result.wave=wave
  result.pca=pca
  outfile = outdir + 'gal_eigenspec_Yip2004.fits'
  mwrfits, result, outfile, /create

  ; QSO PCA
  pca_dir=outdir+'dr1qso/global/'
  tag='qso_eigenspec_Yip2004_global'
  result={wave:dblarr(2582L), pca:dblarr(2582L, 50L)}
  pca=dblarr(2582L, 50L)
  for i=0,49 do begin
    file=pca_dir+'eigSpec_qso_'+string(i,format='(i0)')+'.dat'
    readcol,file,format='d,d',wave, flux,/silent
    pca[*,i] = flux
  endfor
  result.wave=wave
  result.pca=pca
  outfile=outdir+tag+'.fits'  
  mwrfits, result, outfile, /create

  tag=['AZBIN4','AZBIN5','BZBIN2','BZBIN3','BZBIN4','BZBIN5','CZBIN1','CZBIN2', $
       'CZBIN3','CZBIN4','DZBIN1','DZBIN2']
  for jj=0L, n_elements(tag) - 1 do begin
    pca_dir=outdir+'dr1qso/mizbinned/'+tag[jj] + '/'
    readcol, pca_dir+'eigSpec_qso_0.dat',format='d',wave,/silent
    npix=n_elements(wave)
    result={wave:dblarr(npix), pca:dblarr(npix, 50L)}
    pca=dblarr(npix, 50L)
    for i=0,49 do begin
      file=pca_dir+'eigSpec_qso_'+string(i,format='(i0)')+'.dat'
      readcol,file,format='d,d',wave, flux,/silent
      pca[*,i] = flux
    endfor
    result.wave=wave
    result.pca=pca
    outfile=outdir+'qso_eigenspec_Yip2004_'+tag[jj]+'.fits'
    mwrfits, result, outfile, /create
  endfor

end
