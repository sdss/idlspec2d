; compile the gri synthetic magnitudes/flux for RM targets
; using pt corrections from PrepSpec
pro rm_get_spec_synmag_one, rmid, ep=ep, prepspec_dir=prepspec_dir, outdir=outdir
  
  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1, /silent)
  epoch_info=mrdfits(file,2)
  file='/data3/yshen/ftp/sdssrm/collab/science_data/sample_char/sample_char_final.fits'
  char=mrdfits(file,1,/silent)

  nnn=n_elements(rmid)
  calibdir = 'wh_skysub/'

  if ~keyword_set(outdir) then outdir = prepspec_dir + '/specmag/'
  if file_test(outdir) eq 0 then spawn, 'mkdir ' + outdir

  if ~keyword_set(ep) then ep=lindgen(78)+1
  nep = n_elements(ep)
  result={RMID:0L, zfinal:0.D, zsys:0.D,objc_type:0L, MJD:dblarr(nep), gri_mag:dblarr(3,nep), $
          gri_mag_err:dblarr(3,nep), $
          gri_fnu:dblarr(3,nep), gri_fnu_err:dblarr(3,nep), $
          pt:dblarr(nep), pt_err:dblarr(nep), $
          gri_mag_ori:dblarr(3,nep), gri_mag_ori_err:dblarr(3,nep), $
          gri_fnu_ori:dblarr(3,nep), gri_fnu_ori_err:dblarr(3,nep) }

  for i=0, nnn-1 do begin
     out = result
     out.zsys = char[rmid[i]].zsys
     out.rmid = rmid[i] & out.mjd = epoch_info[ep-1].mean_mjd
     out.zfinal=target[rmid[i]].zfinal & out.objc_type=target[rmid[i]].objc_type

     for j=0, nep-1 do begin
        plate=(target[rmid[i]].plate)[ep[j]-1]
        fiber=(target[rmid[i]].fiberid)[ep[j]-1]
        mjd=(target[rmid[i]].mjd)[ep[j]-1]
        if mjd gt 56837 then run2d='v5_10_10/' else run2d='v5_7_1/'
        rm_readspec,plate,fiber, mjd=mjd, $
          calibdir=calibdir, wave=wave,flux=flux,invvar=ivar, run2d=run2d
    
        ; get synthetic AB magnitude and flux fnu (erg/s/cm2/Hz)
        synmag=rm_spec2mag(wave,flux,ivar,synmag_err=synmag_err, $
          synflux=synflux, synfnu_err=synfnu_err)
        out.gri_mag(*,j)=synmag[1:3]
        out.gri_mag_err(*,j)=synmag_err[1:3]
        out.gri_fnu(*,j)=synflux[1:3]
        out.gri_fnu_err(*,j)=synfnu_err[1:3]
     endfor
     out.gri_mag_ori = out.gri_mag
     out.gri_mag_ori_err = out.gri_mag_err
     out.gri_fnu_ori = out.gri_fnu
     out.gri_fnu_ori_err = out.gri_fnu_err

     ; apply the pt corrections from prepspec
     pt=rm_get_pt(rmid[i],err=pt_err,mjd_all=epoch_info[ep-1].mean_mjd,topdir=prepspec_dir)
     out.pt = pt & out.pt_err=pt_err
     for j=0,nep-1 do begin
        out.gri_fnu(*,j)=(out.gri_fnu)[*,j]/pt[j]
        out.gri_fnu_err(*,j)=sqrt(((out.gri_fnu_err)[*,j]/pt[j])^2 + $
           ((out.gri_fnu)[*,j]*pt_err[j]/pt[j]^2)^2)
        out.gri_mag(*,j)=-2.5*alog10((out.gri_fnu)[*,j]) - 48.6
        out.gri_mag_err(*,j)=2.5*(out.gri_fnu_err)[*,j]/((out.gri_fnu)[*,j]*alog(10.D))
     endfor

    outfile = outdir + 'rm' + string(rmid[i],format='(i3.3)') + '.fits'
    mwrfits, out, outfile, /create
    print, 'finished rmid=',i
  endfor

end

pro rm_make_spec_synmag, nep=nep, outfile=outfile, prepspec_dir=prepspec_dir

  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  epoch_info=mrdfits(file,2)

  rmid=indgen(920)
  nnn=n_elements(rmid)
  calibdir='wh_skysub/'
  if ~keyword_set(outfile) then outfile='/data3/yshen/ftp/sdssrm/collab/photo_lc/spec_syn.fits'

  ; define the output spec-synthetic LC structure
  if ~keyword_set(nep) then nep=32L
  result={RMID:0L, zfinal:0.D, objc_type:0L, MJD:dblarr(nep), gri_mag:dblarr(3,nep), $
          gri_mag_err:dblarr(3,nep), $
          gri_fnu:dblarr(3,nep), gri_fnu_err:dblarr(3,nep), $
          pt:dblarr(nep), pt_err:dblarr(nep), $
          gri_mag_ori:dblarr(3,nep), gri_mag_ori_err:dblarr(3,nep), $
          gri_fnu_ori:dblarr(3,nep), gri_fnu_ori_err:dblarr(3,nep) }
  result=replicate(result,nnn)
  result.rmid=rmid
  result.zfinal=target[rmid].zfinal
  result.objc_type=target[rmid].objc_type
  for i=0L, nnn-1 do begin
    result[i].mjd=epoch_info[0:31].mean_mjd
    
    rm_readspec, (target[i].plate)[0:nep-1], (target[i].fiberid)[0:nep-1], mjd=(target[i].mjd)[0:nep-1], $
       calibdir=calibdir, wave=wave,flux=flux,invvar=ivar, run2d=run2d

    ; get synthetic AB magnitude
    for j=0L,nep-1 do begin
       synmag=rm_spec2mag(wave[*,j],flux[*,j],ivar[*,j],synmag_err=synmag_err, $
         synflux=synflux, synfnu_err=synfnu_err)
       result[i].gri_mag(*,j)=synmag[1:3]
       result[i].gri_mag_err(*,j)=synmag_err[1:3]
       result[i].gri_fnu(*,j)=synflux[1:3]
       result[i].gri_fnu_err(*,j)=synfnu_err[1:3]
    endfor
    result[i].gri_mag_ori = result[i].gri_mag
    result[i].gri_mag_ori_err = result[i].gri_mag_err
    result[i].gri_fnu_ori = result[i].gri_fnu
    result[i].gri_fnu_ori_err = result[i].gri_fnu_err

    ; apply the pt corrections from prepspec
    pt=rm_get_pt(i,err=pt_err,mjd_all=epoch_info[0:nep-1].mean_mjd,topdir=prepspec_dir)
    result[i].pt = pt & result[i].pt_err=pt_err
    for j=0,nep-1 do begin
       ;flux[*,j]=flux[*,j]/pt[j] & ivar[*,j]=ivar[*,j]*pt[j]^2
       
       result[i].gri_fnu(*,j)=(result[i].gri_fnu)[*,j]/pt[j]
       result[i].gri_fnu_err(*,j)=sqrt(((result[i].gri_fnu_err)[*,j]/pt[j])^2 + $
          ((result[i].gri_fnu)[*,j]*pt_err[j]/pt[j]^2)^2)
       result[i].gri_mag(*,j)=-2.5*alog10((result[i].gri_fnu)[*,j]) - 48.6
       result[i].gri_mag_err(*,j)=2.5*(result[i].gri_fnu_err)[*,j]/((result[i].gri_fnu)[*,j]*alog(10.D))

    endfor
    print, 'finished rmid=',i

    ;wavevec = wave
    ;flambda2fnu = wavevec^2 / 2.99792e18
    ;fthru = filter_thru(flux * flambda2fnu, waveimg=wavevec, /toair,mask=(ivar LE 0))
    ;thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)  ; ab magnitude
    ;result[i].gri_mag=thismag[*,1:3]
    ;result[i].gri_fnu=10.D^((result[i].gri_mag + 48.6)*(-0.4))
    
  endfor

  mwrfits, result, outfile, /create

end
