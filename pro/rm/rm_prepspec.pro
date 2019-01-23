;+
; NAME:
;   rm_prepspec
; PURPOSE:
;   Generate RMS and mean spectra
;----------------------------------
pro rm_prepspec, mjdlist=mjdlist,calibdir=calibdir,id_coadd=id_coadd,$
    meanspec=meanspec,rmsspec=rmsspec,wave=wave,mederr=mederr,optimal=optimal,$
    maxerr=maxerr, outlier_rej=outlier_rej, percent=percent

   topdir = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') + '/'
   if n_elements(calibdir) eq 0 then calibdir='recalib/'

   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)
   ; If specified, only coadd these objects
   if n_elements(id_coadd) gt 0 then fibermap = fibermap[id_coadd]

   nobj=n_elements(fibermap)
   platearr = fibermap.plate
   fiberarr = fibermap.fiberid
   mjdarr = fibermap.mjd
   ind = where(platearr[*,0] gt 0, nnn)
   platearr = platearr[ind,0]
   fiberarr = fiberarr[ind,*]
   mjdarr = mjdarr[ind,0]
   ; If specified, only use these epochs
   if n_elements(mjdlist) gt 0 then begin
      flag = lonarr(nnn)
      for j=0L,n_elements(mjdlist) - 1L do begin
         indd = where(mjdarr eq mjdlist[j])
         if indd[0] ne -1 then flag[indd] = 1L
      endfor
      ind_mjd = where(flag eq 1)
      platearr = platearr[ind_mjd]
      fiberarr = fiberarr[ind_mjd,*]
      mjdarr = mjdarr[ind_mjd]
   endif
   nepoch = n_elements(platearr)

   ; Get the individual spectra
   ; Enforce the wave array is the same as in the new reduction
   rm_readspec,platearr[0],1,mjd=mjdarr[0],wave=wave,calibdir='recalib/'
   npix=n_elements(wave)

   wave_all = dblarr(npix,nobj,nepoch)
   for i=0,nepoch-1 do begin
      for j=0,nobj-1 do wave_all[*,j,i]=wave
   endfor
   spec_all = dblarr(npix,nobj,nepoch)
   ivar_all = dblarr(npix,nobj,nepoch)
   ferr_all = dblarr(npix,nobj,nepoch)
   meanspec = dblarr(npix,nobj)
   rmsspec = dblarr(npix,nobj)

   for i=0L, nepoch - 1L do begin

      rm_readspec,platearr[i], fiberarr[i,*],mjd=mjdarr[i],calibdir=calibdir, $
         flux=flux,invvar=ivar,flerr=flerr,wave=wave1

      if calibdir ne 'recalib/' and calibdir ne 'wh_skysub/' then begin
         for j=0L, nobj-1 do begin
            spec_all[*,j,i] = interpol(flux[*,j],wave1[*,0],wave)
            ivar_all[*,j,i] = interpol(ivar[*,j],wave1[*,0],wave)
            ferr_all[*,j,i] = interpol(flerr[*,j],wave1[*,0],wave)
         endfor
      endif else begin
         spec_all[*,*,i] = flux
         ivar_all[*,*,i] = ivar
         ferr_all[*,*,i] = flerr
      endelse
   endfor

   ; --------------------------------------
   ; Now compute the mean and rms spectra 
   ; This is the traditional way

   ; reject outliers if required using bspline_iterfit
   if keyword_set(outlier_rej) then begin
      lower=5. & upper=5.
      binsz = alog10(wave[1]) - alog10(wave[0])
      bkptbin = 1.2*binsz
      bkpt=0
      maxrej=ceil(nepoch*npix*nobj*0.1)
      sset = bspline_iterfit(alog10(wave_all),spec_all,invvar=ivar_all,nord=3, $
        bkspace=bkptbin,bkpt=bkpt,lower=lower,upper=upper,outmask=outmask, maxiter=10, $
        /silent,maxrej=maxrej)
      ivar_all = ivar_all*outmask
   endif

   sum1 = total(spec_all*ivar_all, 3, /double)
   sum2 = total(ivar_all, 3, /double)
   ; only deal with pixels where there is at least one good pixel
   ind = where(sum2 gt 0)
   meanspec[ind] = sum1[ind]/sum2[ind]

   devi_all = spec_all
   for i=0L, nepoch - 1 do devi_all[*,*,i]=spec_all[*,*,i] - meanspec
   sum1 = total(devi_all^2*ivar_all, 3, /double)
   rmsspec[ind] = sqrt(sum1[ind]/sum2[ind])

   if keyword_set(percent) then begin ; using 68% percentile to estimate rms
      for i=0L, nobj-1 do $
      rmsspec[*,i]=quantile_2d(0.68, abs(devi_all[*,i,*]), dim=2)
   endif

   mederr = median(ferr_all,dim=3)
   maxerr = max(ferr_all,dim=3)
   ;---------------------------------------

   ; --------------------------------------
   ; Optimal way of calculating mean and rms spectra and their errors
   ; Use Keith Horne's algorithm



   ; --------------------------------------

   ;message,'stop'
end
