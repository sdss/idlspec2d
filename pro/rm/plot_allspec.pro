; compare the 2014 and 2015 SDSS-RM coadded spectra

pro plot_allspec, rmid=rmid, figfile=figfile, ep=ep, run2d=run2d, $
     coadd_mjd = coadd_mjd, years=years, normalize=normalize


  ; Read in the master file
  target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
  fibermap = mrdfits(target_file,1,/silent)

  if ~keyword_set(ep) then ep = indgen(32) + 1
  nep=n_elements(ep) ; # of epochs

  if n_elements(rmid) eq 0 then rmid=indgen(1000L)
  nobj=n_elements(rmid)

  if ~keyword_set(figfile) then $
    figfile='/data3/yshen/spectro/bossredux/v5_7_1/0000/wh_skysub/spec_2014-17.ps'
  begplot, name=figfile, /color
  ang = string(197B)
  pos1 = [0.1, 0.55, 0.97, 0.98]
  pos2 = [0.1, 0.08, 0.97, 0.51]
  xrange=[3650.,10200.]
  xtitle=textoidl('Wavelength [')+ang+']'
  ytitle = textoidl('Flux Density f_\lambda [10^{-17} erg s^{-1} cm^{-2} ') $
      + ang + textoidl('^{-1}]')
  charsize=0.8

  ; setup colors
  cgloadct, 13, /silent, ncolors=nep, bottom=1
  thick=2

  ; coadd spectrum colors
  colors = cgcolor(['black', 'red', 'green', 'cyan'])

  ;mjd1=56837L & mjd2=57196L & mjd3=57576L & mjd4=57934L
  if ~keyword_set(coadd_mjd) then begin
     coadd_mjd=[56837L,57196L,57576L,57934L]
     years = ['2014','2015','2016','2017','2018']
  endif
  for i=0L, nobj-1 do begin

    fiber=rmid[i]+1

    rm_readspec, 0, fiber, mjd=coadd_mjd[0], calibdir='wh_skysub/', wave=wave,flux=flux1,flerr=err1
    indgood=where(err1 gt 1d-5 and wave gt 3600. and wave lt 1d4)
    if indgood[0] ne -1 then $
      yrange=[min(median(flux1[indgood],55) ), max(median(flux1[indgood],55)) ]

    plot, [0],[0], /nodata,xrange=xrange,/xsty,yrange=yrange,pos=pos1

    for jj=0L, n_elements(coadd_mjd) - 1 do begin
      rm_readspec, 0, fiber, mjd=coadd_mjd[jj], calibdir='wh_skysub/', wave=wave,flux=flux1,flerr=err1
      oplot, wave, median(flux1,5),thick=thick, color=colors[jj]
      xyouts, 0.75, 0.93-jj*0.02, /norm, years[jj], color=colors[jj]
    endfor
    xyouts, 0.18, 0.93, /norm, 'RMID'+string(rmid[i],format='(i3.3)')
    xyouts, 0.6, 0.93, /norm, 'z='+string(fibermap[rmid[i]].zfinal,format='(f0.3)')
    xyouts, 0.39, 0.93, /norm, fibermap[rmid[i]].sourcetype

    if ~keyword_set(normalize) then yrange1=yrange else yrange1=[-2,10]
    plot, [0],[0], xrange=xrange,yrange=yrange1,/xsty,pos=pos2, /noerase, xtitle=xtitle,/nodata
    xpos=0.1 & ypos=pos2[3]-0.03
    for j=0,nep-1 do begin
       plate=(fibermap[rmid[i]].plate)[ep[j]-1] & fiber=(fibermap[rmid[i]].fiberid)[ep[j]-1] & mjd=(fibermap[rmid[i]].mjd)[ep[j]-1]
       str=string(plate,format='(i4.4)')+'-'+string(fiber,format='(i3.3)')+'-'+string(mjd,format='(i5.5)')
       rm_readspec, plate,fiber,mjd=mjd,calibdir='wh_skysub/',wave=wave,flux=flux,flerr=err, run2d=run2d
       flux_s = median(flux,15)
       if j ne 6 then begin 
         if ~keyword_set(normalize) then oplot, wave, flux_s, color=j +1L,thick=thick-1 $
          else oplot, wave, flux_s/median(flux_s), color=j +1L,thick=thick-1
       endif
       xpos_new = (j/12)*0.15 + xpos
       ypos_new = ypos - (j mod 12)*0.012
       xyouts, xpos_new, ypos_new, /norm, str, color=j+1L,charsize=charsize
       
    endfor
    xyouts, 0.05, 0.5, /norm, ytitle, align=0.5, orient=90

    splog, 'Finished plotting RMID=', rmid[i]
  endfor

  cgloadct, 0
  endplot
  ; spawn, 'gzip '  + figfile 

end
