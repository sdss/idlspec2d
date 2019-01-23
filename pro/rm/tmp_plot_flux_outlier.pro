; plot an object with catastropic flux calibration errors

pro tmp_plot_flux_outlier

target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
fibermap = mrdfits(target_file,1,/silent)

; this is the one with a catastropic flux calibration error
RM_ID=126
plate = fibermap[rm_id].plate
fiber = fibermap[rm_id].fiberid
mjd = fibermap[rm_id].mjd
ind = where(plate gt 0, nnn)
plate=plate[ind] & fiber=fiber[ind] & mjd=mjd[ind]
ind = where(mjd eq 56722)

rm_readspec, plate,fiber,mjd=mjd,wave=wave,flux=flux

figfile=getenv('IDLRM_DIR') + '/misc/flux_anomaly.ps'
begplot,name=figfile,/color
!p.multi=[0,1,2]

plot, [0], [0], xrange=[3600,1d4],/xsty, yrange=[0,30],xtitle='Observed Wavelength [A]', $
   ytitle=textoidl('Flux Density [erg s^{-1}cm^{-2}A^{-1}]'),/nodata
for i=0, nnn-1 do begin
  oplot, wave[*,i], flux[*,i], psym=3
endfor
oplot, wave[*,ind],flux[*,ind],color=cgcolor('red')
xyouts, 0.6,0.9,string(plate[ind],format='(i4.4)')+'-'+string(mjd[ind],format='(i5.5)')+ $
   '-' + string(fiber[ind],format='(i4.4)'),color=cgcolor('red'),/norm
xyouts, 0.2,0.9, 'RM_ID='+string(RM_ID,format='(i0)'),/norm

; find one subexp of the affect epoch (56722)
subexp='b1-00176329'
fiber0 = fiber[ind] ; affect object
; find another two object's fiberid that is closet
spherematch, fibermap[ind].ra,fibermap[ind].dec,fibermap.ra,fibermap.dec, 1., match0, match1, distance,maxmatch=0
indd = where(distance gt 0.)
match0=match0[indd] & match1=match1[indd] & distance=distance[indd]
indd=sort(distance)
match0=match0[indd] & match1=match1[indd] & distance=distance[indd]
ind1=match1[0] & ind2=match1[1]
fiber1 = (fibermap[ind1].fiberid)[ind]
fiber2 = (fibermap[ind2].fiberid)[ind]

; this is the flat field
flatfile=getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') $
       + '/7339/' +'spFlat-b1-00176328.fits.gz'
fiberflat=mrdfits(flatfile,0)
npix = (size(fiberflat))[1]

plot,[0],[0],xtitle='Pixel', ytitle = 'FiberFlat', xrange=[800, 3300.], yrange=[0.6,1.2],/xsty,/ysty,$
  title='spFlat-b1-00176328'
for i=0L,499 do oplot, findgen(npix), fiberflat[*,i],color=cgcolor('gray'),thick=0.2
oplot, findgen(npix), fiberflat[*,fiber0-1], color=cgcolor('red')
oplot, findgen(npix), fiberflat[*,fiber1-1], color=cgcolor('cyan')
oplot, findgen(npix), fiberflat[*,fiber2-1], color=cgcolor('blue')
xyouts, 0.4, 0.42, string(plate[ind],format='(i4.4)')+'-'+string(mjd[ind],format='(i5.5)')+ $
   '-' + string(fiber1,format='(i4.4)') + ', dist='+string(distance[0]*3600.,format='(f0.1)')+' arcsec',color=cgcolor('cyan'),/norm 
xyouts, 0.4, 0.39, string(plate[ind],format='(i4.4)')+'-'+string(mjd[ind],format='(i5.5)')+ $
   '-' + string(fiber2,format='(i4.4)') + ', dist='+string(distance[1]*3600.,format='(f0.1)')+' arcsec',color=cgcolor('blue'),/norm

!p.multi=0
endplot

end
