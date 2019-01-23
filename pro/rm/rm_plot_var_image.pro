; plot a 2D image of the varbility of the 849 SDSS-RM sample

pro rm_plot_var_image, maxvalue=maxvalue, minvalue=minvalue, log=log, psplot=psplot, $
   tag=tag, type=type, figfile=figfile, infile=infile

  if n_elements(maxvalue) eq 0 then maxvalue=2.5 
  if n_elements(minvalue) eq 0 then minvalue=0.
  range=[minvalue, maxvalue]
  if ~keyword_set(tag) then tag='total RMSx spectrum'
  
  if keyword_set(infile) then file = infile else $
    file='/data3/yshen/ftp/sdssrm/collab/prepspec/working/nov24_2015/var_img.fits'
  result=mrdfits(file,1)
  zz=result.z & wave=result.wave
  if n_elements(type) eq 0 then type=0
  if type eq 0 then img=result.img_rms_norm
  if type eq 1 then img=result.img_avg_norm
  if type eq 2 then img=result.IMG_FRAC
  if type eq 3 then img=result.img_rmsc_norm

  trim_npix=100L 
  img=img[trim_npix:*, *] & wave=wave[trim_npix:*]
  ind=sort(zz)
  zz=zz[ind]
  img=img[*,ind]
  if keyword_set(log) then img=alog10(img)
  nx=n_elements(wave) & ny=n_elements(zz)

  angstr=textoidl('\AA')
  if keyword_set(psplot) then begin
    if keyword_set(figfile) then fig=figfile else $
    fig='/data3/yshen/ftp/sdssrm/collab/prepspec/working/nov24_2015/totrms_img.ps'
    begplot, name=fig,/color,/landscape
    angstr=string(197B)
  endif

  lines=textoidl(['Ly\alpha', 'SiIV', 'CIV', 'CIII]', 'MgII', 'HeII4687', 'H\beta', 'H\alpha']) ; , $
;  'H\gamma', 'H\delta','H8/HeI', 'FeII:','HeI'])
  linewave=[1215.67, (1396.76 + 1402.06)*0.5, 1549.06,1908.73,2798.75,4687.02, 4862.68,6564.61, $
  4341.68, 4102.89, 3890., 3800., 3588.30]
  linepix=[800,700,600,500, 400, 150, 80, 30, $
    200,250,200,300,300]
  zline=interpol(zz, indgen(ny), linepix)
  linewave_pix=interpol(indgen(nx), wave, linewave*(1+zline))
  nline=n_elements(lines)
  xrange=[0,nx-1] & yrange=[0,ny-1]
  pos=[0.08,0.1,0.92,0.9]
  title='RMS Variability'
  charsize=1.5
  plot,[0],[0], xrange=xrange, yrange=yrange, xsty=5,ysty=5, pos=pos, /norm, /nodata
  cgLoadCT, 63, reverse=reverse  ; color table: 55
  cgimage, img, maxvalue=maxvalue, minvalue=minvalue, pos=pos,/norm, /noerase
  cgColorbar, Position=[pos[0], pos[3], pos[2], pos[3]+0.05], range=range,/top
  cgLoadCT, 0

  for j=0, nline-1 do begin
    ;print, linewave_pix[j], linepix[j], lines[j]
    xyouts, linewave_pix[j], linepix[j], lines[j], color=cgcolor('magenta'),charsize=charsize
  endfor

  ; map the linear wavelength and z grid onto the pixelized image
  zgrid=[0, 1, 2, 3, 4]
  zgrid_pix=interpol(indgen(ny), zz, zgrid)
  wgrid=[4000, 5000, 6000, 7000, 8000, 9000, 10000.]
  wgrid_pix=interpol(indgen(nx), wave, wgrid)

  axis, xaxis=0, xrange=xrange, /xsty, xtitle=textoidl('Observed Wavelength (')+angstr+')', xtickv=wgrid_pix, $
     xtickname=string(wgrid, format='(i0)'), xticks=n_elements(wgrid) - 1
  axis, yaxis=0, yrange=yrange, /ysty, ytitle='Redshift', ytickv=zgrid_pix, $
     ytickname=string(zgrid, format='(i0)'), yticks=n_elements(zgrid) - 1
  axis, yaxis=1, yrange=yrange, /ysty, ytitle='# of Quasars'
  xyouts, pos[0], pos[1]-0.08, tag, /norm,color=cgcolor('red')

  if keyword_set(psplot) then begin
    endplot
    cgfixps, fig
  endif

end

pro rm_make_var_image, prepspec_dir=prepspec_dir, $
      result=result, check=check, outfile = outfile

  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1)

  target=target[0:848]
  zz=target.zfinal
  nobj=n_elements(target)
  
  if ~keyword_set(prepspec_dir) then $
     prepspec_dir='/data3/yshen/ftp/sdssrm/collab/prepspec/working/nov24_2015/'
  ; prefix of rms broad lines
  blrtag=['blr_w.dat', 'ha_w.dat', 'hb_w.dat', 'he2_4686_w.dat', 'mg2_w.dat', $
        'c3_w.dat', 'c4_w.dat', 'si4_w.dat', 'n5_w.dat', 'lya_w.dat']
  nblr=n_elements(blrtag)

  ; create a 2D array
  rm_readspec, 0, 1, mjd=56837,wave=wave
  npix=n_elements(wave)

  image_frac=dblarr(npix, nobj)
  image_avg=dblarr(npix, nobj)
  image_avg_norm=dblarr(npix, nobj)
  image_rms=dblarr(npix, nobj)
  image_rms_norm=dblarr(npix, nobj)
  image_rmsc=dblarr(npix, nobj)
  image_rmsc_norm=dblarr(npix, nobj)

  for i=0L, nobj-1 do begin

    rmtag='rm'+string(i,format='(i3.3)')
    rmdir=prepspec_dir + rmtag + '/'
    avgfile=rmdir + rmtag + '_avg_w.dat'
    rmsfile=rmdir + rmtag + '_rms_w.dat'
    splog, 'cycle through: ', rmtag

    if file_test(avgfile) eq 1 and file_test(rmsfile) eq 1 then begin
      readcol, avgfile, format='d,d', wave1, flux1, /silent

      if n_elements(wave1) gt 100. then begin
        wave1=wave1*(1. + zz[i])
        flux_avg=interpol(flux1, wave1, wave) > 0

        readcol, rmsfile, format='d,d,x,d', wave1, flux1, flux1_RMSc, /silent
        ; flux1_RMSc only include continuum model RMS, and add model BLR rms now
        for jj=0L, nblr-1 do begin
           extrafile=rmdir + rmtag + '_' + blrtag[jj]
           if file_test(extrafile) eq 1 then begin
              readcol, extrafile, format='x,d', flux_blr, /silent
              if n_elements(flux_blr) eq n_elements(flux1_RMSc) then $
                 flux1_RMSc=flux1_RMSc + flux_blr
              if keyword_set(check) then splog, 'adding BLR rms'
           endif
        endfor

        wave1=wave1*(1. + zz[i])
        flux_rms=interpol(flux1, wave1, wave) > 0  ; RMSx, maximum-likelihood estimate (error and pt corrected)
        flux_rmsc=interpol(flux1_RMSc, wave1, wave) > 0  ; RMSc, model, error-accounted, pt corrected

        frac=flux_rms / flux_avg > 0

        image_frac[*, i]=frac
        image_avg[*,i]=flux_avg
        image_avg_norm[*,i]=flux_avg/median(flux_avg)
        image_rms[*,i]=flux_rms
        image_rms_norm[*,i]=flux_rms/median(flux_rms)
        image_rmsc[*,i]=flux_rmsc
        image_rmsc_norm[*,i]=flux_rmsc/median(flux_rmsc)

        ;plot, wave1, flux1, xrange=[4000,8000]
        ;oplot, wave, flux_rms, color=cgcolor('red')
        ;message, 'diagonize'

        if keyword_set(check) then begin
          plot, image_rms_norm[*,i], yrange=[0, 3]
          oplot, image_rmsc_norm[*,i], color=cgcolor('red')
          pause
        endif

      endif  
    endif


  endfor
  ind=where(finite(image_frac) eq 0 or image_frac gt 100.)
  if ind[0] gt 0 then begin
     image_frac[ind] = 0.
     image_avg[ind] = 0.
     image_avg_norm[ind] = 0.
     image_rms[ind]=0.
     image_rms_norm[ind]=0.
     image_rmsc[ind]=0
     image_rmsc_norm[ind]=0
  endif

  result={z:zz, wave:wave, img_frac:image_frac, img_avg:image_avg, img_avg_norm:image_avg_norm,$
     img_rms:image_rms, img_rms_norm:image_rms_norm, img_rmsc:image_rmsc, img_rmsc_norm:image_rmsc_norm}
  if ~keyword_set(outfile) then $
      outfile='/data3/yshen/ftp/sdssrm/collab/prepspec/working/nov24_2015/var_img.fits'
  mwrfits, result, outfile,/create

end
