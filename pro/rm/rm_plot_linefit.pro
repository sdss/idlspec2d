;+
; NAME:
;   rm_plot_linefit
;
; PURPOSE:
;   Plot the spectal line windows along with the line fit
;
; CALLING SEQUENCE:
;   rm_plot_linefit, plate, fiber, mjd, [/psplot,figfile=figfile]
;   rm_plot_linefit, qsofit=qsofit

pro rm_plot_linefit, plate, fiber, mjd, zsys=zsys, qsofit=qsofit, $
      emparfile=emparfile, psplot=psplot, figfile=figfile, calibdir=calibdir

   cs=2.9979246d5

   if keyword_set(qsofit) then begin
      plate=qsofit.plate & fiber=qsofit.fiberid & mjd=qsofit.mjd
      zsys=qsofit.z
   endif
  
   nobj=n_elements(plate)
   ; default directory for the coadded RM spectra (plate=0)
   if ~keyword_Set(calibdir) then calibdir='wh_skysub/'

   ; read the line structure
   if keyword_set(emparfile) then emline_file=emparfile else $
    emline_file = getenv('IDLRM_DIR')+'/etc/qsoline_all.par'
   yanny_read, emline_file, pdata
   linelist_all = *pdata[0]
   yanny_free, pdata

   ; determine the plot xy positions for each line
   nline_all=n_elements(linelist_all)
   if keyword_set(psplot) then begin
     charsize=1. & thick = 3
   endif else begin
     charsize=1.5 & thick=1
   endelse

   if keyword_Set(psplot) then begin
      if ~keyword_set(figfile) then figfile1='lineplot.ps' else $
        figfile1=figfile
      begplot, name=figfile1,/color
   endif

   ; Start the for loop
   for i=0, nobj-1 do begin

      rm_readspec,plate[i], fiber[i],mjd=mjd[i],calibdir=calibdir, $
        wave=wave, flux=flux,flerr=flerr
      wave=wave/(1. + zsys[i])

      if keyword_set(qsofit) then begin
         ; find out which lines were fit
         alltag=tag_names(qsofit[i])
         fitflag=lonarr(nline_all)
         for jj=0, nline_all - 1 do begin
            ind_tmp=where(alltag eq strupcase(linelist_all[jj].linename))
            if ind_tmp ne -1 then begin
               if (qsofit[i].(ind_tmp))[0] gt 0 then fitflag[jj]=1
            endif
         endfor
         ind_fit=where(fitflag eq 1)
         if ind_fit[0] ne -1 then linelist=linelist_all[ind_fit]
         print, 'Line fit:', linelist.linename
      endif else begin
         linelist=linelist_all
      endelse
      nline=n_elements(linelist)

      plot_layout, nline, xypos=xypos, pmargin=[0.04,0.06], omargin=[0.08, 0.08, 0.98, 0.95]
      for j=0,nline - 1 do begin
         ind=where(wave ge linelist[j].minwav and wave le linelist[j].maxwav,nnn)
         if nnn gt 10 then begin
            xrange=linelist[j].voff*cs*[-1,1]/1000.
            flux_s=median(flux[ind],3)
            ; yrange=[min(flux_s), max(flux_s)]
            vwave=(wave[ind]-linelist[j].lambda)/linelist[j].lambda*cs/1000.
            ind_v = where(abs(vwave) le linelist[j].voff*cs/1000. )
            yrange=[min(flux_s[ind_v]), max(flux_s[ind_v])]

            if j eq 0 then noerase=0 else noerase=1
            plot, vwave,flux[ind], xrange=xrange,yrange=yrange,noerase=noerase,pos=xypos[*,j], $
              ytickname=replicate(' ', 10L),charsize=charsize, thick=thick, psym=10
            oplot, [0,0], !y.crange, line=1
            ind_tmp=where(alltag eq strupcase(linelist[j].linename))
            wave_fit=(qsofit[i].(ind_tmp))[0]
            vwave_fit=(wave_fit-linelist[j].lambda)/linelist[j].lambda*cs/1000.
            oplot, [vwave_fit,vwave_fit],!y.crange,line=2,color=cgcolor('green'),thick=thick
            xyouts, xypos[0,j]+0.02,xypos[3,j]-0.025,/norm,linelist[j].linename,charsize=charsize
          endif
      endfor
      xyouts, 0.5, 0.01, textoidl('Velocity [10^3 km s^{-1}]'), align=0.5, /norm, charsize=charsize
      xyouts, 0.5, 0.96, 'RMID-'+string(qsofit[i].rm_id,format='(i3.3)')+' z='+string(zsys[i],format='(f6.4)'), $
        /norm, align=0.5, charsize=charsize

      if ~keyword_set(psplot) then pause
   endfor

   if keyword_set(psplot) then endplot

end
