;+
; rm_prep_skysub
; Collection of routines to prepare for the improved sky subtraction
; using the Wild&Hewett scheme

;---- generate mask for sky pixels for the WH scheme
pro rm_gen_skymask, diag=diag

  ; set the default spPlate* dir
  calibdir='recalib/'

  ;file='/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  ind=where(target[0].plate gt 0, nep)
  indd=where(strtrim(target.sourcetype) eq 'SKY', nobj)
  sky=target[indd]
  plate=(sky.plate)[ind,*] & fiber=(sky.fiberid)[ind,*] & mjd=(sky.mjd)[ind,*]
  nsky=nobj*nep
  plate=reform(plate,nsky) & fiber=reform(fiber,nsky) & mjd=reform(mjd,nsky)

  ; read in all sky spectra
  rm_readspec,plate,fiber,mjd=mjd,calibdir=calibdir,loglam=loglam,flux=flux,flerr=flerr,invvar=ivar
  loglam=loglam[*,0]
  lam=10.D^loglam
  npix=n_elements(lam)

  ; get the 67% abs(flux) and med err; the former is the rms flux
  rms_flux=quantile_2d(0.67,abs(flux),dim=2)
  med_err=median(flerr,dim=2)

  ratio=(rms_flux/(med_err) ) > 0
  ; fix the bad pixels with err=0 at the red and blue edge
  ind=where(finite(ratio) ne 1)
  ratio[ind] = 1.

  ; find the sky pixels
  ; first fit a low-order spline function to it
  nord=2
  fullbkpt = bspline_bkpts(loglam, everyn=100, $
    bkpt=0, nord=nord)
  sset = bspline_iterfit(loglam, ratio, $
    invvar=replicate(100., npix), lower=0, upper=1, fullbkpt=fullbkpt, $
    maxrej=100, outmask=outmask1, nord=nord, yfit=yfit)
   
  ;ind_sky=where(lam gt 6300. and lam lt 10300. and ratio gt rmin, npix_sky)
  if ~keyword_set(wavemin) then wavemin=4000.
  if ~keyword_set(wavemax) then wavemax=10300.
  ind_sky=where(outmask1 eq 0 and lam ge wavemin and lam le wavemax, npix_sky)

  ; output the skymask [loglam at the sky pixels]
  outfile=getenv('IDLRM_DIR')+'/etc/wh-sky-mask.txt'
  fmt='(f6.4, " ", f9.3)'
  openw, lun, outfile, /get_lun
  printf, lun, '#loglam  lam'
  for i=0,npix_sky-1 do printf,lun,format=fmt,loglam[ind_sky[i]],lam[ind_sky[i]]
  close, lun
  free_lun, lun

  ; plot a diagnostic plot
  if keyword_set(diag) then begin
    figfile=getenv('IDLRM_DIR')+'/etc/wh-sky-mask.eps'
    begplot, name=figfile,/color,xsize=10,ysize=6,/encap 
    xrange = [3500,10600]
    plot, lam, rms_flux, xrange=xrange,/xsty,yrange=[0,4],xtitle='Obs Wavelength'
    oplot,lam, med_err, color=cgcolor('red')
    oplot, lam, ratio,psym=5
    oplot, lam, yfit, color=cgcolor('green'),thick=2
    oplot,lam[ind_sky],ratio[ind_sky],psym=5,color=cgcolor('cyan')
    xyouts, 4000, 3.5, 'rms sky flux'
    xyouts, 4000, 3.2, 'med sky err',color=cgcolor('red')
    xyouts, 6000, 3.5, 'spline fit of rms/med_err ratio', color=cgcolor('green')
    xyouts, 6000, 3.2, 'sky pixels', color=cgcolor('cyan')
    endplot
  endif

end

;------- rescaling of the sdss errors ---
; instead of using eqn. 1 of WH05, use the formula in getweights.pro in the package 
; distributed by V. Wild
pro gen_err_rescaling, plate, mjd, calibdir=calibdir, alpha=alpha, beta1=beta1, $
     narr=narr, carr=carr, sarr=sarr, max_n_c=max_n_c, diag=diag, _EXTRA=Extra

  if ~keyword_set(calibdir) then calibdir='recalib/'
  ; Two preset parameters used in the recaling
  if ~keyword_set(alpha) then alpha=1
  if ~keyword_set(beta1) then beta1=0.3
  
  ; get the median noise from sky spectra in this plate
  rm_readspec,plate,mjd=mjd,plugmap=plugmap,calibdir=calibdir,/silent, _EXTRA=Extra

  ; On 10/29/2015, changed sourcetype tag 'SKY' to 'SKY' or 'NA'
  ; as sky fibers in plugmap are typed 'NA'
  ; this bug caused all the 2014 plates to use the median error for the 1000 objects
  ; rather than for the 80 sky to estimate the error scaling, but it turns out that
  ; these two median errors are very similar. 
  fiber=where( (strtrim(plugmap.sourcetype) eq 'SKY') or (strtrim(plugmap.sourcetype) eq 'NA' ) ) + 1L
  rm_readspec,plate,fiber,mjd=mjd,wave=lam,loglam=loglam,flux=flux,flerr=flerr,calibdir=calibdir,/silent,_EXTRA=Extra
  ; all wavelength arrays are identical
  loglam=loglam[*,0] & lam=lam[*,0]
  npix=n_elements(lam)

  ; this is the median noise array
  narr=median(flerr,dim=2)
  carr=narr

  ; get the sky and non-sky pixels
  skyfile = getenv('IDLRM_DIR')+'/etc/wh-sky-mask.txt'
  mask=rm_skymask(lam, margin=1L, skyfile=skyfile,fmt='x,d')
  ind = where(mask eq 0, complement=indd)
  if ind[0] ne -1 then carr[ind] = interpol(carr[indd], lam[indd], lam[ind] )

  ind_good=where(narr gt 0 and lam lt 10300. and lam gt 4000.)
  max_n_c = max( (narr - carr)[ind_good] )

  ; this is to be multiplied by the pipeline reported noise
  sarr = 1. - ( (narr - carr) / max_n_c )^alpha*beta1  < 1;

  if keyword_Set(diag) then begin
    plot, narr, psym=5
    oplot, carr, thick=2,color=cgcolor('cyan')

  endif

end
; ---- make a fits file to store all the scaled error arrays for the RM plates
pro make_wh_rescale_err

  file='/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  plate=target[0].plate & mjd=target[0].mjd
  ind=where(plate gt 0, nep)
  plate=plate[ind] & mjd=mjd[ind]

  rm_readspec, plate[0],1,mjd=mjd[0],calibdir='recalib/',wave=wave
  npix=n_elements(wave)
  result={plate:0L, mjd:0L, alpha:0.D, beta1:0.d, narr:dblarr(npix), carr:dblarr(npix),sarr:dblarr(npix), $
        max_n_c:0.D }
  result=replicate(result, nep)
  result.plate=plate & result.mjd=mjd
  alpha=1 & beta1=0.3
  result.alpha=alpha & result.beta1=beta1
  for i=0L, nep - 1 do begin
    gen_err_rescaling, plate[i], mjd[i], alpha=alpha, beta1=beta1, narr=narr, carr=carr, sarr=sarr, $
        max_n_c=max_n_c
    result[i].narr=narr & result[i].carr=carr & result[i].sarr=sarr & result[i].max_n_c=max_n_c
  endfor

  outfile = '/data3/quasar/yshen/ftp/bossredux/v5_7_1/wh_skysub/rescal_err.fits'
  mwrfits, result, outfile, /create

end


; -------- PCA of the sky spectra
pro gen_sky_pca, result=result, loglam=loglam, eigenspec=eigenspec, $
  nkeep=nkeep,niter=niter

; set the default spPlate* dir
  calibdir='recalib/'

  file='/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  ; reject the two epochs with the lowest sn
  ind=where(target[0].plate gt 0 and target[0].mjd ne 56669 and target[0].mjd ne 56713, nep)
  platelist=(target[0].plate)[ind]
  mjdlist=(target[0].mjd)[ind]
  indd=where(strtrim(target.sourcetype) eq 'SKY', nobj)
  sky=target[indd]
  plate=(sky.plate)[ind,*] & fiber=(sky.fiberid)[ind,*] & mjd=(sky.mjd)[ind,*]
  nsky=nobj*nep
  plate=reform(plate,nsky) & fiber=reform(fiber,nsky) & mjd=reform(mjd,nsky)
  
  ; read in all sky spectra
  rm_readspec,plate,fiber,mjd=mjd,calibdir=calibdir,loglam=loglam,flux=flux,flerr=flerr,invvar=ivar
  ; store the original noise array
  ivar_ori = ivar
  loglam=loglam[*,0]
  lam=10.D^loglam
  npix=n_elements(lam)
  splog, 'Total # of sky spectra: ', nsky

  ; rescale the pipeline error
  sarr_all = dblarr(npix, nep)
  alpha=1 & beta1=0.3
  for i=0L, nep - 1L do begin
    gen_err_rescaling, platelist[i], mjdlist[i], alpha=alpha, beta1=beta1, sarr=sarr,max_n_c=max_n_c
    sarr_all[*,i] = sarr

    ;print, platelist[i], mjdlist[i], 'max_n_c=', max_n_c, max(sarr)
    ;plot, sarr
    ;pause

    ind = where(plate eq platelist[i] and mjd eq mjdlist[i], nnn)
    for j=0l, nnn-1 do begin
      flerr[*, ind[j]] = flerr[*, ind[j]]*sarr
      ivar[*,ind[j]] = ivar[*,ind[j]]/sarr^2
      ;print, ind[j], '  ivar [min,max]=', min(ivar[*,ind[j]]), max(ivar[*,ind[j]])
      ;pause
    endfor
  endfor
  ;message, 'stop'
  splog, 'Finished rescaling noise arrays'

  ; keep only the sky pixel
  skyfile = getenv('IDLRM_DIR')+'/etc/wh-sky-mask.txt'
  mask=rm_skymask(lam, margin=1L, skyfile=skyfile, fmt='x,d') 
  ind_sky = where(mask eq 0, complement=ind_nonsky, n_skypix)

  ; remove objects with ivar=0 at any of the sky pixels
  rejflag=lonarr(nsky)
  for i=0L, nsky - 1 do begin
     if total(ivar[ind_sky, i] gt 0) ne n_skypix then rejflag[i]=1
  endfor
  ind=where(rejflag eq 0, nsky_good)
  splog, 'rejected bad sky spectra: ', nsky - nsky_good
  plate=plate[ind] & fiber=fiber[ind] & mjd=mjd[ind]
  nsky=nsky_good
  flux=flux[*,ind] & ivar=ivar[*,ind] & ivar_ori=ivar_ori[*,ind]

  ;message, 'stop'

  ; do the PCA
  if ~keyword_set(nkeep) then nkeep=n_elements(ind_sky)
  if ~keyword_set(niter) then niter=1
  ; use pca_solve to get eigenspectra
  ;eigenspec = pca_solve(flux[ind_sky,*], ivar[ind_sky,*], $
  ;  nkeep=nkeep, eigenval=eigenval, niter=niter, acoeff=acoeff)

  ; alternative way to get eigenspectra
  ; should I replace flux[ind_sky,*] with flux[ind_sky,*]*sqrt(ivar[ind_sky,*])?
  ; my test shows that it doesn't really matter which eigenvector basis set is used. 
  ;pca, transpose(flux[ind_sky,*]), eigenval,eigenspec,percentages, proj_obj, /silent
  pca, transpose(flux[ind_sky,*]*sqrt(ivar[ind_sky,*])), eigenval,eigenspec,percentages, proj_obj, /silent
  acoeff=flux[ind_sky,*] ## eigenspec
  acoeff=transpose(acoeff)
  eigenspec=transpose(eigenspec)

  ; create an output structure
  result = {plate:plate,fiber:fiber,mjd:mjd,nsky:nsky, niter:niter, $
     eigenspec:eigenspec, eigenval:eigenval, acoeff:acoeff, $
     loglam:loglam, flux:flux, ivar_ori:ivar_ori, ivar:ivar, skymask:mask, $
     ind_sky:ind_sky, ind_nonsky:ind_nonsky}

  ; pca on flux
  ;outfile='/data3/quasar/yshen/ftp/bossredux/v5_7_1/wh_skysub/pca.fits'
  ; pca on flux/rescaled noise
  outfile='/data3/quasar/yshen/ftp/bossredux/v5_7_1/wh_skysub/pca_norm.fits'
  mwrfits, result, outfile, /create

end
;

