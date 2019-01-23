;PURPOSE:
;  plot comparison of spectra before and after the WH05 improved sky subtraction

pro plot_sky_impro,plate,mjd,rmid=rmid, psplot=psplot

  ; get the master file
  target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
  fibermap = mrdfits(target_file,1,/silent)

  if n_elements(rmid) eq 0 then rmid=0
  target=fibermap[rmid]  ; this contains target info

  if plate ne 0 then begin ; this is a single epoch
    plate_list=target.plate
    mjd_list=target.mjd
    fiber_list=target.fiberid
    ind=where(plate_list eq plate and mjd_list eq mjd)
    fiber=fiber_list[ind]
  endif else fiber=rmid+1 ; this is the coadded plate

  ;read in spectra
  if plate ne 0 then $
    rm_readspec, plate,fiber,mjd=mjd, calibdir='recalib/',wave=wave,flux=flux,flerr=err else $
    rm_readspec, plate,fiber,mjd=mjd, calibdir='',wave=wave,flux=flux,flerr=err
  rm_readspec, plate,fiber,mjd=mjd, calibdir='wh_skysub/',wave=wave1,flux=flux1,flerr=err1
  minflux=min(flux1[where(wave gt 3800 and wave le 10000)]) 
  maxflux=max(flux1[where(wave gt 3800 and wave le 10000)])
  med_sn1=median(flux1/err1)

  yrange=[-2, maxflux*1.2]

  ; get sky pixels
  skymask=rm_skymask(wave)

  ; now plot the figure
  ang=textoidl('\AA')
  if keyword_set(psplot) then begin
    figfile=getenv('IDLRM_DIR') + '/misc/wh_skysub_rmID_'+string(rmid,format='(i3.3)')+'.eps'
    begplot, name=figfile, /color,/encap,/cmyk,xsize=8,ysize=8
    ang=string(197B)
    thick=3.
  endif

  pos=[0.1,0.55,0.95,0.98]
  plot, wave,flux,xrange=[3600,10300],yrange=yrange,/xsty, /ysty, $
    pos=pos, xtickname=replicate(' ', 10L),thick=thick
  oplot, wave, err, color=cgcolor('red'),thick=thick
  xyouts, 0.15, 0.93, 'original spectrum',/norm
  ;xyouts, 0.15, 0.93, 'RMID='+string(rmid,format='(i0)'),/norm
  pstr=string(plate,format='(i4.4)')+'-'+string(fiber,format='(i4.4)')+'-'+string(mjd,format='(i5.5)')
  xyouts, 0.4, 0.93, pstr,/norm
  xyouts, 0.8, 0.93, 'i='+string(target.psfmag[3],format='(f6.3)'),/norm
  ind=where(skymask eq 0, nnn)
  oplot, wave[ind], replicate(0.5*yrange[1], nnn), psym=symcat(5),color=cgcolor('cyan'),symsize=0.1

  pos=[0.1,0.1,0.95,0.53]
  plot, wave1,flux1,xrange=[3600,10300],yrange=yrange,xtitle=textoidl('Wavelength [')+ang+']', $
    /noerase, pos=pos, /xsty, /ysty,thick=thick
  oplot, wave1, err1, color=cgcolor('red'),thick=thick
  xyouts, 0.03, 0.5, align=0.5, /norm, orient=90, textoidl('f_{\lambda} [10^{-17}erg s^{-1}cm^{-2}')+ang+textoidl('^{-1}')+']'
  xyouts, 0.15, 0.48, 'improved sky subtraction', /norm
  xyouts, 0.60, 0.48, 'median S/N/pixel='+string(med_SN1,format='(i0)'),/norm

  if keyword_set(psplot) then begin
    endplot
  endif

end
