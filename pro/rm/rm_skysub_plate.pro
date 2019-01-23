; subtract sky residual using the Wild&Hewett05 method
; for a given plate
; this only works for the relibrated data with fixed wavelength coverage
; 
; .comp rm_prep_skysub
; to run for the 2015 eBOSS data with new reduction v5_10_10, do 
; rm_skysub_plate, topdir='/data3/yshen/spectro/bossredux/v5_10_10/', run2d='v5_10_10' 

pro rm_skysub_plate, plate, mjd, topdir=topdir, diag=diag, outdir_tag=outdir_tag, $
     pcafile=pcafile, _Extra=extra

  calibdir='recalib/'
  if ~keyword_set(topdir) then topdir=getenv('BOSS_SPECTRO_REDUX') + '/' + $
    getenv('RUN2D')+'/'

  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  platelist=target[0].plate
  mjdlist=target[0].mjd

  if not keyword_set(outdir_tag) then outdir_tag='wh_skysub/'

  ; set up eigenspec common block used in subtractoh.pro
  common subtractOH1, espec, lambda, pix_sky, pix_nosky
  ;file='/data3/quasar/yshen/ftp/bossredux/v5_7_1/wh_skysub/pca.fits'
  file=getenv('IDLRM_DIR')+'/template/sky_residual_pca.fits'
  if keyword_Set(pcafile) then file=pcafile
  pca=mrdfits(file,1)
  lambda=pca.loglam & espec=pca.eigenspec & pix_sky=pca.ind_sky & pix_nosky=pca.ind_nonsky
  ; cull out a red-wavelength region
  ind_tmp=where(10.^lambda gt 5000., ncomplement=ncomp)
  lambda=lambda[ind_tmp]
  pix_sky=pix_sky - ncomp & pix_nosky=pix_nosky - ncomp
  ; remove the nosky pixels below 5000 A
  ind_tmp=where(pix_nosky ge 0)
  pix_nosky=pix_nosky[ind_tmp]

  ; start the sky subtraction process
  common plate_weights,sarr, weights_in
  for i=0L, n_elements(plate) - 1 do begin

    ; generate the by-plate weights (to be used by subtractOH)
    gen_err_rescaling, plate[i], mjd[i], narr=narr, carr=carr, sarr=sarr,_Extra=extra
    weights_in = 1./(narr*sarr)
    ;weights_in = getweights(plate[i], './')

    ; Added on 10/29/2015 to fix the sky subtraction issue for 2015 SDSS plates
    ; namely, the median error on the plate at some sky pixels could be zero, which
    ; is different from the 2014 plates
    ; make sure there are no INF weights_in due to narr=0 (median error for the plate)
    ind_inf=where(finite(weights_in) eq 0,complement=ind_f)
    if ind_inf[0] ne -1 then weights_in[ind_inf]=median(weights_in[ind_f])

    pstr=string(plate[i],format='(i4.4)')
    mstr=string(mjd[i],format='(i5.5)')
    outdir=topdir+pstr+'/' + outdir_tag
    if file_test(outdir) eq 0 then spawn, 'mkdir ' + outdir
    
    ; get the fits headers of spPlate* 
    spfile=topdir+pstr+'/'+calibdir+'spPlate-'+pstr+'-'+mstr+'.fits'
    rm_readspec,plate[i],1,mjd=mjd[i],wave=wave,loglam=loglam,calibdir=calibdir,_Extra=extra
    HDU0=mrdfits(spfile,0,hd0) ; flux
    HDU1=mrdfits(spfile,1,hd1) ; invvar
    HDU2=mrdfits(spfile,2,hd2) ; AND-pixelmask
    HDU3=mrdfits(spfile,3,hd3) ; OR-puxelmask
    HDU4=mrdfits(spfile,4,hd4) ; dispersion map
    HDU5=mrdfits(spfile,5,hd5) ; plugmap
    HDU6=mrdfits(spfile,6,hd6) ; sky

    ; we will only replace the flux and invvar images
    npix=(size(HDU0,/dim))[0] & nobj=(size(HDU0,/dim))[1]
    flux_arr=fltarr(npix,nobj) & ivar_arr=fltarr(npix,nobj)
    err_arr_in=HDU1*0.
    ind_tmp=where(HDU1 gt 0)
    err_arr_in[ind_tmp]=1./sqrt(HDU1[ind_tmp])

    ; loop over all fibers
    iep=where(platelist eq plate[i] and mjdlist eq mjd[i])
    for jj=0L, nobj-1 do begin

       qso=0 & star=0 & sky=0
       type=strtrim(target[jj].sourcetype)
       if strmatch(type,'RM*') ne 0 then QSO=1 else begin
          if type eq 'STD' then STAR=1 else SKY=1
       endelse

       fiber=(target[jj].fiberid)[iep] & z=target[jj].zfinal
       flux=HDU0[*, fiber-1] & err=err_arr_in[*,fiber-1]

       ; call the modified subtractoh to get the cleaned spectrum
       ; note that the sky residual eigenvectors were created on the set of
       ; sky spectra unnormalized by the scaled errors, which is likely different
       ; from Wild&Hewett in which the eigenvectors were created on the normalized
       ; sky spectra. However, the two approaches give negligible difference, so
       ; I am sticking to my approach
       newflux=subtractOH(flux,err,loglam,z,dummyplate, $
        rms=rms,nrecon=nrecon,weights=weights_in, qso=qso,star=star, sky=sky,_Extra=extra)
       if keyword_set(diag) then pause else $
         splog, 'finished: Plate, mjd, iobj:', plate[i], mjd[i], jj

       newerr=err*sarr
       newivar=0.*newerr
       ind_tmp=where(newerr gt 0)
       if ind_tmp[0] ne -1 then newivar[ind_tmp] = 1./(newerr[ind_tmp])^2
       flux_arr[*,fiber-1]=newflux
       ivar_arr[*,fiber-1]=newivar
    endfor

    ; now output the new plate
    outfile=outdir+'spPlate-'+pstr+'-'+mstr+'.fits'
    splog, 'outfile:', outfile
    mwrfits, flux_arr, outfile, hd0, /create
    mwrfits, ivar_arr, outfile, hd1
    mwrfits, HDU2, outfile, hd2
    mwrfits, HDU3, outfile, hd3
    mwrfits, HDU4, outfile, hd4
    mwrfits, HDU5, outfile, hd5
    mwrfits, HDU6, outfile, hd6

  endfor


end
  
