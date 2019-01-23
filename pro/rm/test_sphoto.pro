; Suite to test the spectrophotometry

;--------------
; Compare synmag to real-time photometry
pro comp_realtime_photo

   ; Find all the targets (quasar and stars)
   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)
   
   ; remove the one LRG and sky fibers
   ind = where(strmatch(fibermap.sourcetype,'LRG*') or strmatch(fibermap.sourcetype,'SKY*') $
             , complement = indd, ncomp=nobj)
   fibermap = fibermap[indd]

   ; compute the synmag for all targets
   ; We will only focus on epoch 7339-56686, for which we have overlap photometry from Bok
   specpath = getenv('BOSS_SPECTRO_REDUX') $
           + '/' + getenv('RUN2D') + '/7339/recalib/test20/'
   fiberid = fibermap.fiberid[4,*]
   synmag = fltarr(nobj, 5L)
   readspec, 7339, fiberid, mjd=56686, wave=wave, $
           flux=flux,invvar=invvar,andmask=mask,path=specpath,plugmap=plugmap,/silent
   flambda2fnu = wave^2 / 2.99792e18
   fthru = filter_thru(flux * flambda2fnu, waveimg=wave,/toair)
   synmag = -2.5 * alog10(fthru) - (48.6-2.5*17) >0 <99
   
   ; this is the default pipeline synmag
   synmag_pip = fltarr(nobj, 5L)
   readspec, 7339, fiberid, mjd=56686, wave=wave, $
           flux=flux,invvar=invvar,andmask=mask,plugmap=plugmap,/silent
   flambda2fnu = wave^2 / 2.99792e18
   fthru = filter_thru(flux * flambda2fnu, waveimg=wave,/toair)
   synmag_pip = -2.5 * alog10(fthru) - (48.6-2.5*17) >0 <99

   ; Readin the realtime photometry from Bok on the same night
   file = getenv('IDLRM_DIR') + '/etc/20140128.g.cat.fits'
   gmag = mrdfits(file,1)
   file = getenv('IDLRM_DIR') + '/etc/20140129.i.cat.fits'
   imag = mrdfits(file,1)

   ; now start the plot
   figfile = getenv('IDLRM_DIR')+'/misc/test_sphoto_realtime.ps'
   begplot, name=figfile,/color
   !p.multi = [0,2,2]

   minmag = 19. ; limit for bright objects

   ; there seems to be a constant offset between calibmag and the psfmag
   ; calibmag = 22.5 - alog10(calibflux)*2.5 = psfmag + offset
   ; offset = [-0.042, 0.036, 0.015, 0.013, -0.002]
   ; These offsets are 2002 version, see kcorrect/v4_1_4/pro/utils
   ; the calibflux in plugmap seems to be using the 2002 version of ab_offset

   ; First plot the gmag distribution
   spherematch, fibermap.ra,fibermap.dec,gmag.alpha_j2000,gmag.delta_j2000,2./3600.D, $
       match0,match1, distance
   magdiff = synmag[match0,1] - (gmag[match1].mag)  ; + 0.036)
   ind = where(gmag[match1].mag lt minmag and gmag[match1].mag gt 0, ntot )
   magdiff_bright = magdiff[ind]
   magdiff_pip = synmag_pip[match0,1] - (gmag[match1].mag) ; + 0.036)
   magdiff_bright_pip = magdiff_pip[ind]
   tmp = where(abs(magdiff_bright_pip) le 0.065, n_pip_good)
   tmp = where(abs(magdiff_bright) le 0.065, n_new_good)

   xrange = [-0.3,0.3] & bin = 0.02 & csize=1.0

   plothist, magdiff_pip, xrange = xrange, xtitle='SynMag - BokMag [g]' $
      , ytitle=textoidl('N_{obj}'), pos=[0.15,0.57,0.5,0.95], bin=bin
   plothist, magdiff, /over, color=fsc_color('red'),bin=bin
   xyouts, 0.35, 0.92, 'Default redux', /norm,charsize=csize
   xyouts, 0.35, 0.90, 'xyfit redux', /norm, color=fsc_color('red'),charsize=csize
   plothist, magdiff_bright_pip, xrange = xrange, xtitle='SynMag - BokMag [g<'+string(minmag,format='(i0)')+']' $
     , ytitle=textoidl('N_{obj}'), pos = [0.61,0.57, 0.96,0.95],bin=bin
   plothist, magdiff_bright, /over, color=fsc_color('red'),bin=bin
   xyouts, 0.65, 0.92, '<6%',/norm,charsize=csize
   xyouts, 0.65, 0.90, string(n_pip_good,format='(i0)')+'/'+string(ntot,format='(i0)'),charsize=csize,/norm
   xyouts, 0.65, 0.88, string(n_new_good,format='(i0)')+'/'+string(ntot,format='(i0)'),charsize=csize,color=fsc_color('red'),/norm

   ; now plot the imag distribution
   spherematch, fibermap.ra,fibermap.dec,imag.alpha_j2000,imag.delta_j2000,2./3600.D, $
       match0,match1, distance
   magdiff = synmag[match0,3] - (imag[match1].mag + 0.013)
   ind = where(imag[match1].mag lt minmag and imag[match1].mag gt 0, ntot)
   magdiff_bright = magdiff[ind]
   magdiff_pip = synmag_pip[match0,3] - (imag[match1].mag + 0.013)
   magdiff_bright_pip = magdiff_pip[ind]
   tmp = where(abs(magdiff_bright_pip) le 0.065, n_pip_good)
   tmp = where(abs(magdiff_bright) le 0.065, n_new_good)

   plothist, magdiff_pip, xrange = xrange,xtitle='SynMag - BokMag [i]' $
      , ytitle=textoidl('N_{obj}'),pos=[0.15, 0.1, 0.5, 0.48],bin=bin
   plothist, magdiff, /over, color=fsc_color('red'),bin=bin
   plothist, magdiff_bright_pip, xrange = xrange,xtitle='SynMag - BokMag [i<'+string(minmag,format='(i0)')+']', $
      ytitle=textoidl('N_{obj}'),pos=[0.61,0.1,0.96,0.48],bin=bin
   plothist, magdiff_bright, /over, color=fsc_color('red'),bin=bin
   xyouts, 0.65, 0.45, '<6%',/norm,charsize=csize
   xyouts, 0.65, 0.43, string(n_pip_good,format='(i0)')+'/'+string(ntot,format='(i0)'),charsize=csize,/norm
   xyouts, 0.65, 0.41, string(n_new_good,format='(i0)')+'/'+string(ntot,format='(i0)'),charsize=csize,color=fsc_color('red'),/norm

   !p.multi = 0
   endplot

   ;message, 'stop and diagnose'

end
