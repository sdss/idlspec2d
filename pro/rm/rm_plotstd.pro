;+
; NAME:
;   rm_plotstd
;
; PURPOSE: 
;   Compare the final coadded std spectrum and the best-fit model during 
;   the spflux stage
;  
;   rm_plotstd,7338,0,56660 to plot all good std


pro rm_plotstd,plate,fiber,mjd,psplot=psplot,recalibdir=recalibdir $
      ,subexp=subexp,tag=tag,wave1=wave1,modelflux=modelflux,oldflux=oldflux,newflux=newflux

   ; default std to look at: 7338, 2, 56660

   ;if not keyword_set(recalibdir) then $
   ;  recalibdir='/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/recalib/test4/'

   ; First plot the default BOSS reduction results
   readspec,plate,mjd=mjd,wave=wave,flux=flux,plugmap=plugmap

   adderr = 0.03

   if n_elements(psplot) eq 0 then psplot=1L 
   if keyword_set(psplot) then begin
  
      plotfile = 'plotstdspec_' + string(plate,format='(i0)') + '-' + $
         string(mjd,format='(i5.5)') + '_' + tag + '.ps'
      set_plot, 'ps'
      dfpsplot, djs_filepath(plotfile, root_dir=getenv('IDLRM_DIR')+'/misc'),$
          /color, /landscape
      thick = 3
      
   endif

   ; find the best Kurucz model
   ; check spectrograph 1 and 2
   rm_expinfo, plate, mjd, expinfo=expinfo
   ind = where(expinfo.coadded[0,*] eq 1, nexp)
   exp_all = expinfo[ind].name
   calibfile_all = repstr(exp_all,'spFrame','spFluxcalib')
   exp12 = [(expinfo[ind[0]].name)[0], (expinfo[ind[0]].name)[1]]
   calibfile = repstr(exp12,'spFrame','spFluxcalib')

   topdir = getenv('BOSS_SPECTRO_REDUX')
   twoddir = getenv('RUN2D')
   platestr = string(plate,format='(i4.4)')   
   mjdstr = string(mjd,format='(i5.5)')

   for icam=0,1 do begin
      filename = lookforgzip(filepath(calibfile[icam], root_dir=topdir, $
         subdirectory=[twoddir,platestr]), count=ct)
      if ct eq 1 then filename=filename[0]
      splog, filename

      kindx=mrdfits(filename,2)
      if fiber gt 0 then ind_std = where(kindx.fiberid eq fiber) $
       else ind_std = indgen(n_elements(kindx))
   
      for i=0L, n_elements(ind_std) - 1L do begin

         fiberid = kindx[ind_std[i]].fiberid

         fiberstr = string(fiberid,format='(i4.4)')
         wave1=wave[*,fiberid-1] & flux1=flux[*,fiberid-1]
         oldflux = flux1
   
         ; estimate the default filterd flux (only use gri)
         wavevec = wave1
         flambda2fnu = wavevec^2 / 2.99792e18
         fthru = filter_thru(flux1 * flambda2fnu, waveimg=wavevec, /toair)
         thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)
         default_filter_flux = 10.^((22.5-thismag)/2.5)

         imodel = kindx[ind_std[i]].imodel
         npix = n_elements(wave1)
         tmploglam = alog10(wave1) & tmpdispimg = replicate(1d-4,npix)
         tmpflux = spflux_read_kurucz(tmploglam, tmpdispimg, $
             iselect=imodel)
         extcurve2 = ext_odonnell(10.^tmploglam, 3.1)
            tmpflux = tmpflux $
             * 10.^(-extcurve2 * 3.1 * plugmap[fiberid-1].ebv / 2.5)
         wavevec = 10.d0^tmploglam
         flambda2fnu = wavevec^2 / 2.99792e18
         fthru = filter_thru(tmpflux * flambda2fnu, waveimg=wavevec, /toair)
         thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)
         scalefac = plugmap[fiberid-1].calibflux[2] / 10.^((22.5-thismag[2])/2.5)  
         modelflux = tmpflux*scalefac
         
         ; this a somewhat larger wavelengh range
         tmploglam = 3.4780d0 + lindgen(5620) * 1.d-4
         tmpdispimg = 0 * tmploglam + 1.
         tmpflux = spflux_read_kurucz(tmploglam, tmpdispimg, $
             iselect=imodel)
         extcurve2 = ext_odonnell(10.^tmploglam, 3.1)
            tmpflux = tmpflux $
             * 10.^(-extcurve2 * 3.1 * plugmap[fiberid-1].ebv / 2.5)
         wavevec = 10.d0^tmploglam
         flambda2fnu = wavevec^2 / 2.99792e18
         fthru = filter_thru(tmpflux * flambda2fnu, waveimg=wavevec, /toair)
         thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)
         scalefac = plugmap[fiberid-1].calibflux[2] / 10.^((22.5-thismag[2])/2.5)
         model_filter_flux = 10.^((22.5-thismag)/2.5) * scalefac  

         stdflux = tmpflux * scalefac

         yrange = [0,max(stdflux)*1.2]
         plot, wave1, flux1, xrange=[3500,1d4],/xsty,xtitle='Wavelength', $
         ytitle=textoidl('Flux Density [10^{-17} ergs^{-1}cm^{-2}\AA^{-1}]'), $
         yrange=yrange

         oplot, 10.^tmploglam, stdflux, color=cgcolor('blue')

         ; now overplot the re-calibrated spectrum
         if keyword_set(recalibdir) then begin
            readspec,plate,fiberid,mjd=mjd,wave=wave2,flux=flux2,path=recalibdir
            oplot, wave2, flux2, color=cgcolor('cyan')
            xyouts,0.6,0.88, 'New redux',/norm, color=cgcolor('cyan')
            newflux = flux2

            ; calculate the recalibread flitered flux (only use gri)
            wavevec = wave2
            flambda2fnu = wavevec^2 / 2.99792e18
            fthru = filter_thru(flux2 * flambda2fnu, waveimg=wavevec, /toair)
            thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)
            recalib_filter_flux = 10.^((22.5-thismag)/2.5)

            ; plot the individual subexposures
            if keyword_set(subexp) then begin

               colors=replicate(cgcolor(['magenta']),nexp)

               blue_name = reform(calibfile_all[0, *])
               red_name = reform(calibfile_all[2,*])

               for iexp=0L, nexp - 1L do begin
               
                  ; do it for the blue spectrograph
                  bluefile = recalibdir+blue_name[iexp]+'.gz'
                  spframefile = lookforgzip(filepath(exp_all[0,iexp], root_dir=topdir, $
                     subdirectory=[twoddir,platestr]), count=ct)
                  spframe_read,spframefile,fiberid-1,loglam=loglam,objflux=objflux, $
                    wset=wset,objivar=objivar,dispimg=dispimg,adderr=adderr
                  calibfac=mrdfits(bluefile,0,/silent)
                  ; re-normalize the flux and dispersion
                  correct_dlam, objflux, objivar, wset, dlam=dloglam
                  correct_dlam, dispimg, 0, wset, dlam=dloglam, /inverse
                  ; flux calibration
                  minval = 0.05 * mean(calibfac)
                  divideflat,objflux,invvar=objivar,calibfac[*,fiberid-1],minval=minval
                  oplot,10.^loglam,objflux,color=colors[iexp mod n_elements(colors)],psym=3

                  ; do the same for the red spectrograph
                  redfile = recalibdir+red_name[iexp]+'.gz'
                  spframefile = lookforgzip(filepath(exp_all[2,iexp], root_dir=topdir, $
                     subdirectory=[twoddir,platestr]), count=ct)
                  spframe_read,spframefile,fiberid-1,loglam=loglam,objflux=objflux, $
                    wset=wset,objivar=objivar,dispimg=dispimg,adderr=adderr
                  calibfac=mrdfits(redfile,0,/silent)
                  ; re-normalize the flux and dispersion
                  correct_dlam, objflux, objivar, wset, dlam=dloglam
                  correct_dlam, dispimg, 0, wset, dlam=dloglam, /inverse
                  ; flux calibration
                  minval = 0.05 * mean(calibfac)
                  divideflat,objflux,invvar=objivar,calibfac[*,fiberid-1],minval=minval
                  oplot,10.^loglam,objflux,color=colors[iexp mod n_elements(colors)],psym=3

                endfor

            endif
  
         endif

         xyouts,0.6, 0.92, 'Default redux',/norm
         xyouts,0.6,0.9, 'Best-fit Kurucz model',/norm,color=cgcolor('blue')
         xyouts,0.2, 0.92, 'std: '+ platestr+'-'+fiberstr+'-'+mjdstr,/norm
         xyouts,0.2, 0.9, 'qgood='+string(kindx[ind_std[i]].qgood,format='(i0)'),/norm
         xyouts,0.18, 0.15, 'calibflux: ', /norm,color=cgcolor('gold'),charthick=thick
         xyouts,0.18,0.13,'model:',/norm,color=cgcolor('blue'),charthick=thick
         for ii=0L, 4L do begin
            xyouts,0.25+ii*0.06, 0.15, string(plugmap[fiberid-1].calibflux[ii],format='(f0.2)'),/norm $
              , color=cgcolor('gold'),charthick=thick
            xyouts,0.25+ii*0.06, 0.13, string(model_filter_flux[ii],format='(f0.2)'),/norm $
             ,color=cgcolor('blue'),charthick=thick
         endfor
         for ii=1,3 do begin
            xyouts,0.25+ii*0.06, 0.11, string(default_filter_flux[ii],format='(f0.2)'),/norm $
             , color=cgcolor('opposite'),charthick=thick
            xyouts,0.25+ii*0.06, 0.09, string(recalib_filter_flux[ii],format='(f0.2)'),/norm $
             , color=cgcolor('cyan'),charthick=thick
         endfor
         redchi2 = kindx[ind_std[i]].CHI2/kindx[ind_std[i]].dof
         xyouts,0.2,0.88,textoidl('\chi^2/dof=')+string(redchi2,format='(f0.1)'),/norm

         ;pause

      endfor
   endfor

   if keyword_set(psplot) then dfpsclose

end
