;+
; NAME:
;   rm_plotstdrms
;
; PURPOSE:
;   plot the rms spectra for the 70 standard stars and compare the pipeline 
;   reduction and the custom reduction;
;   Remove outliers with gri abs(log(synflux/calibflux))>0.1
;
;
;-------------------

pro rm_plotstdrms, interp_sky=interp_sky

   ; read in the master file
   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)

   ; setup the mjd list to search
   platelist = fibermap[0].plate
   mjdlist = fibermap[0].mjd
   ind=where(platelist gt 0)
   platelist=platelist[ind] & mjdlist=mjdlist[ind]
   ; reject 56669 and 56713
   ind = where(mjdlist ne 56669 and mjdlist ne 56713,nepoch)
   platelist=platelist[ind] & mjdlist=mjdlist[ind]
   nstd = 70L
   flag_pipe = lonarr(nepoch,nstd)
   flag_new = lonarr(nepoch,nstd)

   ; find the fiberid of the std and get synflux
   ; reject epoch with gri abs(log(synflux/calibflux))>0.1
   ; 1. for the pipeline reduction
   fiber = rm_findstd(platelist, mjdlist, calibdir='', synflux=synflux, $
           calibflux=calibflux)
   gri_ratio = (alog10(synflux/calibflux))[*,*,[1,2,3]]
   for i=0,nstd-1 do begin
      for j=0,nepoch-1 do begin
         arr=where(abs(gri_ratio[j,i,*]) le 0.1, nnn)
         if nnn eq 3 then flag_pipe[j,i]=1L
      endfor
   endfor
   ; 2. for the new reduction
   fiber = rm_findstd(platelist, mjdlist, calibdir='recalib/', synflux=synflux, $
           calibflux=calibflux)   
   gri_ratio = (alog10(synflux/calibflux))[*,*,[1,2,3]]
   for i=0,nstd-1 do begin
      for j=0,nepoch-1 do begin
         arr=where(abs(gri_ratio[j,i,*]) le 0.1, nnn)
         if nnn eq 3 then flag_new[j,i]=1L
      endfor
   endfor

   ; now compute the rms spectra for stars
   id_std = where(strtrim(fibermap.sourcetype) eq 'STD')
   rm_readspec,platelist[0],1,mjd=mjdlist[0],wave=wave,calibdir='wh_skysub/'  ; 'recalib/'
   npix=n_elements(wave)
   rmsspec_pipe=dblarr(npix,nstd) & meanspec_pipe=dblarr(npix,nstd)
   rmsspec_new=dblarr(npix,nstd) & meanspec_new=dblarr(npix,nstd)
   mederr_all=dblarr(npix,nstd)
   for i=0, nstd-1 do begin

      ; pipeline reduction
      ind = where(flag_pipe[*,i] eq 1)
      mjd_use = mjdlist[ind]
      rm_prepspec, mjdlist=mjd_use,calibdir='',id_coadd=id_std[i],$
       meanspec=meanspec,rmsspec=rmsspec
      meanspec_pipe[*,i] = meanspec
      rmsspec_pipe[*,i] = rmsspec

      ; new reduction
      ind = where(flag_new[*,i] eq 1)
      mjd_use = mjdlist[ind]
      rm_prepspec, mjdlist=mjd_use,calibdir='wh_skysub/',id_coadd=id_std[i],$
       meanspec=meanspec,rmsspec=rmsspec,mederr=mederr
      meanspec_new[*,i] = meanspec
      rmsspec_new[*,i] = rmsspec
      mederr_all[*,i] = mederr
   endfor

   result = {platelist:platelist, mjdlist:mjdlist, flag_pipe:flag_pipe, flag_new:flag_new, $
             meanspec_pipe:meanspec_pipe, rmsspec_pipe:rmsspec_pipe, meanspec_new:meanspec_new, $
             rmsspec_new:rmsspec_new,mederr:mederr_all}
   outfile = getenv('IDLRM_DIR') + '/misc/rmsspec_std.fits'
   mwrfits, result, outfile, /create

   ; make a plot
   rms_mean_ratio_pipe = rmsspec_pipe/meanspec_pipe
   rms_mean_ratio_new = rmsspec_new/meanspec_new
   skymask=rm_skymask(wave, margin=3) ;,skyfile = getenv('IDLRM_DIR')+'/etc/dr9-sky-mask.txt',fmt='d')
   ind=where(skymask gt 0, complement=indd)

   ang = string(197B)
   figfile = getenv('IDLRM_DIR') + '/misc/rmsspec_std.eps'
   begplot, name=figfile, /encap,/color,/cmyk, xsize=6,ysize=4
   thick=3
   xarr=wave & yarr=(median(rms_mean_ratio_pipe,dim=2))
   yarr[indd]=interpol(yarr[ind], xarr[ind], xarr[indd])
   plot, xarr[ind], yarr[ind], xrange=[3600, 1d4],/xsty,yrange=[0,0.2],/ysty, $
     xtitle = textoidl('Wavelength [')+ang+']', ytitle=textoidl('RMS Flux / Mean Flux'),thick=thick, pos=[0.15, 0.165, 0.94, 0.98]
   yarr=(quantile_2d(0.16, rms_mean_ratio_pipe,dim=2))
   yarr[indd]=interpol(yarr[ind], xarr[ind], xarr[indd])
   oplot, xarr[ind], yarr[ind],line=2,thick=thick
   yarr = (quantile_2d(0.84, rms_mean_ratio_pipe,dim=2))
   yarr[indd]=interpol(yarr[ind], xarr[ind], xarr[indd])
   oplot, xarr[ind], yarr[ind], line=2,thick=thick

   yarr=(median(rms_mean_ratio_new,dim=2))
   yarr[indd]=interpol(yarr[ind], xarr[ind], xarr[indd])
   oplot, xarr[ind], yarr[ind],color=cgcolor('cyan'),thick=thick
   yarr=(quantile_2d(0.16, rms_mean_ratio_new,dim=2))
   yarr[indd]=interpol(yarr[ind], xarr[ind], xarr[indd])
   oplot, xarr[ind], yarr[ind],line=2,color=cgcolor('cyan'),thick=thick
   yarr=(quantile_2d(0.84, rms_mean_ratio_new,dim=2))
   yarr[indd]=interpol(yarr[ind], xarr[ind], xarr[indd])
   oplot, xarr[ind], yarr[ind], line=2,color=cgcolor('cyan'),thick=thick

   yarr=(median(mederr_all/meanspec_new,dim=2))
   yarr[indd]=interpol(yarr[ind], xarr[ind], xarr[indd])
   oplot, xarr[ind], yarr[ind], color=cgcolor('red'),thick=thick
   yarr=(quantile_2d(0.16, mederr_all/meanspec_new,dim=2))
   yarr[indd]=interpol(yarr[ind], xarr[ind], xarr[indd])
   oplot, xarr[ind], yarr[ind],line=2,color=cgcolor('red'),thick=thick
   yarr=(quantile_2d(0.84, mederr_all/meanspec_new,dim=2))
   yarr[indd]=interpol(yarr[ind], xarr[ind], xarr[indd])
   oplot,xarr[ind],yarr[ind],line=2,color=cgcolor('red'),thick=thick

   xyouts, 0.18, 0.89, 'pipeline reduction',/norm
   xyouts, 0.18, 0.84, 'custom reduction', /norm, color=cgcolor('cyan')

   endplot

end
