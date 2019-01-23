;+
; NAME:
;   rm_plot_LC_one
;
; PURPOSE:
;   Plot light curves (continuum+line) for one target
;   This is an old version to use my own fits, rather than PrepSpec lightcurves
;
;
; ---------------------------------------------------------
pro rm_plot_LC_one, rm_ID, topdir=topdir,LC_data=LC_data,tag=tag,mjdlist=mjdlist, $
     epoch_id=epoch_id, xrange=xrange, ccf_tau=ccf_tau, ccf_out=ccf_out,more_info=more_info, $
     vw_norm=vw_norm, sncut=sncut

   if n_elements(sncut) eq 0 then sncut=5. ; default SN cut to keep a LC point

   if not keyword_set(topdir) then $
    topdir='/data3/quasar/yshen/spectro/bossredux/v5_6_0/'

   ; readin the LC data
   if not keyword_set(LC_data) then begin
      rm_get_rms_prop,rm_ID,result=LC_data, /silent,epoch_id=epoch_id
   endif
   if not keyword_set(mjdlist) then begin
      mjdlist = LC_data[0].mjd
   endif
   tagname = tag_names(LC_data)
   if not keyword_set(xrange) then xrange=[56650,56850]
   if not keyword_set(yrange) then yrange=[0.2,1.8]

   ; Determine the Van Groningen & Wanders scaling if OIII is covered
   ind=where(strmatch(tagname, strupcase('OIII*') ))
   ngood_OIII=LC_data.(ind[2])
   LC_arr = LC_data.(ind[3])
   LC_err = LC_data.(ind[4])
   goodflag=LC_data.(ind[5])
   nepoch = n_elements(LC_arr)  
   vw_scaling = dblarr(nepoch) + 1.
   if ngood_OIII gt 0 and keyword_set(vw_norm) then begin
      ; find the best epoch as the reference epoch
      ind_good = where(goodflag eq 1)
      minsnr = min(LC_err[ind_good], ind_best)  ; error in dex
      if minsnr lt 0.2 then begin
         ref_ep = epoch_id[ind_good[ind_best]]
         vw_scaling=rm_vw_fluxing_obj(rm_ID,ref_ep=ref_ep,epoch_id=epoch_id)
      endif
   endif

   !P.multi = [0,1,3]
   charsize=2.5 & ticklen=0.05
   ; plot continuum and [OIII] light curves, if available
   pos = [0.11, 0.7, 0.95, 0.95]
   name=['OIII','L1350','L1700','L3000','L5100']
   meanmjd = []
   colors = cgcolor(['opposite','green','red','cyan','magenta'])
   psyms = [4,5,9,15,16]
   if not keyword_set(tag) then tag=''
   title='RM_ID='+string(rm_ID,format='(i0)') + ' ' + tag
   plot,[0],[0],xrange=xrange,yrange=yrange,xtitle='MJD',ytitle='Relative Flux', $
    title=title, pos=pos, charsize=charsize,xticklen=ticklen,xtickformat='(i0)'
   ypos=pos[3] - 0.02
   xpos1= pos[2]-0.8
   xyouts, xpos1,pos[3]-0.03,/norm,'RMS:',charsize=1
   for i=0L, n_elements(name) - 1 do begin
      ind = where(strmatch(tagname, strupcase(name[i])+'*') )
      ngood = LC_data.(ind[2])
      LC_arr = LC_data.(ind[3])
      LC_err = LC_data.(ind[4])
      goodflag = LC_data.(ind[5])
      if ngood gt 0 then begin
         ;LC_rms = LC_data.(ind[0])
         ;ind_good = where(LC_arr gt 0 and LC_err gt 0)
         ind_good = where(goodflag eq 1)
         LC_arr = 10.D^LC_arr[ind_good]*vw_scaling[ind_good] 
         LC_err = LC_err[ind_good] ; NB, err in dex
         mjd_plot = mjdlist[ind_good]
         mean_flux = mean(LC_arr)
         LC_rms = stddev(LC_arr)/mean_flux
         LC_err = LC_err*alog(10.D)*LC_arr/mean_flux ;
         LC_arr = LC_arr/mean_flux ; flux relative to the mean flux over the period
         oploterror,mjd_plot,LC_arr,LC_err,psym=symcat(psyms[i]),color=colors[i]
         legend,pos=[pos[2]-0.18,ypos],name[i],/norm,psym=symcat(psyms[i]),color=colors[i],$
          textcolor=colors[i],box=0,charsize=1.2
         xpos1=xpos1+0.07
         xyouts, xpos1, pos[3]-0.03,/norm,string(LC_RMS,format='(f0.2)'),charsize=1, $
           color=colors[i]
         ypos=ypos-0.02
      endif
   endfor

   ; plot the line LCs if available
   pos = [0.11,0.39, 0.95,0.64]
   name=['Halpha','Hbeta','MgII','CIII','CIV']
   plot,[0],[0],xrange=xrange,yrange=yrange,xtitle='MJD',ytitle='Relative Flux', $
    pos=pos,charsize=charsize,xticklen=ticklen,xtickformat='(i0)'
   ypos=pos[3] - 0.02
   xpos1= pos[2]-0.8
   xyouts, xpos1,pos[3]-0.03,/norm,'RMS:',charsize=1
   for i=0L, n_elements(name) - 1 do begin
      ind = where(strmatch(tagname, strupcase(name[i])+'*') )
      ngood = LC_data.(ind[2])
      LC_arr = LC_data.(ind[3])
      LC_err = LC_data.(ind[4])
      goodflag = LC_data.(ind[5])
      if ngood gt 0 then begin
         ;LC_rms = LC_data.(ind[0])
         ;ind_good = where(LC_arr gt 0 and LC_err gt 0)
         ind_good = where(goodflag eq 1)
         LC_arr = 10.D^LC_arr[ind_good]*vw_scaling[ind_good]
         LC_err = LC_err[ind_good] ; NB, err in dex
         mjd_plot = mjdlist[ind_good]
         mean_flux = mean(LC_arr)
         LC_rms = stddev(LC_arr)/mean_flux
         LC_err = LC_err*alog(10.D)*LC_arr/mean_flux ;
         LC_arr = LC_arr/mean_flux ; flux relative to the mean flux over the period
         oploterror,mjd_plot,LC_arr,LC_err,psym=symcat(psyms[i]),color=colors[i]
         legend,pos=[pos[2]-0.18,ypos],name[i],/norm,psym=symcat(psyms[i]),color=colors[i],$
          textcolor=colors[i],box=0,charsize=1.2
         xpos1=xpos1+0.07
         xyouts, xpos1,pos[3]-0.03,/norm,string(LC_RMS,format='(f0.2)'),charsize=1, $
           color=colors[i]
         ypos=ypos-0.02
      endif
   endfor

   ; plot the CCF
   pos = [0.11,0.08,0.95,0.33]   
   plot,[0],[0],xrange=[-20,100],yrange=[-0.5,1],xtitle='Time Delay (Days)',ytitle='CCF', $
    pos=pos,charsize=charsize,xticklen=ticklen
   ypos=pos[3] - 0.03

   ; generate the time delay array (in observed days) in computing the ccf
   if not keyword_set(ccf_tau) then begin
      dt = 2. ; days
      ccf_tau = -10. + (findgen(150./dt + 5L)*dt)
   endif

   line=['Halpha','Hbeta','MgII','CIII','CIV']
   ; try to pair the line LC with the continuum LC with the shortest wavelength possible
   conti=['L1350','L1700','L3000','L5100']

   ccf_out = dblarr(n_elements(ccf_tau), n_elements(line))
   for i=0L, n_elements(line) - 1L do begin
      ind = where(strmatch(tagname, strupcase(line[i])+'*') )
      ngood = LC_data.(ind[2])
      LC_arr = LC_data.(ind[3])
      LC_err = LC_data.(ind[4])
      goodflag = LC_data.(ind[5])
      if ngood gt 0 then begin
         ;ind_good = where(LC_arr gt 0 and LC_err gt 0)
         ind_good = where(goodflag eq 1)
         LC_arr = 10.D^LC_arr[ind_good]*vw_scaling[ind_good]
         LC_err = LC_err[ind_good] ; NB, err in dex
         LC_err = LC_err*alog(10.D)*LC_arr
         mjd_good = mjdlist[ind_good]
         ; remove noisy measurements?
         ind_goodsn = where(LC_arr/LC_err gt sncut)  ; 20% flux measurements, which is too precise
         cadence_id = mjd_good[ind_goodsn]
         lc_line=LC_arr[ind_goodsn] & err1=LC_err[ind_goodsn]

         ; get continuum LC, starting from the shortest wavelength
         flag0=0 & iconti=0
         while flag0 eq 0 do begin
            ind = where(strmatch(tagname, strupcase(conti[iconti])+'*') )
            ngood = LC_data.(ind[2])
            LC_arr = LC_data.(ind[3])
            LC_err = LC_data.(ind[4])
            goodflag = LC_data.(ind[5])
            if ngood gt 0 then begin
               ;ind_good = where(LC_arr gt 0 and LC_err gt 0)
               ind_good = where(goodflag eq 1)
               LC_arr = 10.D^LC_arr[ind_good]*vw_scaling[ind_good]
               LC_err = LC_err[ind_good] ; NB, err in dex
               LC_err = LC_err*alog(10.D)*LC_arr
               mjd_good = mjdlist[ind_good]
               ; remove noisy measurements?
               ind_goodsn = where(LC_arr/LC_err gt sncut) ; 20% flux measurements, which is too precise
               cadence_id_c = mjd_good[ind_goodsn]
               lc_conti=LC_arr[ind_goodsn] & err2=LC_err[ind_goodsn]
               ; we have found the continuum LC
               flag0=1L
            endif
            iconti=iconti+1
         endwhile 

         ; ccf[*,0] is based on the original data
         ccf = xcorr_interp(cadence_id,lc_line,cadence_id_c,lc_conti,tau=ccf_tau,ndata=ndata $
            , /nogap,err1=err1,err2=err2, /bootstrap,cent_tau=cent_tau, peak_tau=peak_tau $
            , mirror=mirror)

         rm_dcf,cadence_id_c,lc_conti,cadence_id,lc_line,err1=median(err2),err2=median(err1),$
            DCF=DCF,eDCF=eDCF,tau_grid=tau_grid
        
         oplot, ccf_tau, ccf[*,0],color=colors[i] ; psym=symcat(psyms[i])
         ind1 = where(eDCF gt 0)
         ; only plot the Hbeta DCF
         if ind1[0] ne -1 then begin
            oploterror, tau_grid[ind1],dcf[ind1],edcf[ind1],psym=symcat(psyms[i]),color=colors[i]
            ;print, dcf[ind1]
         endif
         tau_mea = [cent_tau[0], quantile_1d(0.16, cent_tau), quantile_1d(0.84, cent_tau) ]
         for jj=0,2 do xyouts, pos[2]-0.18+jj*0.03,ypos, /norm, color=colors[i], $
            string(tau_mea[jj],format='(i0)'),charsize=1
         ypos=ypos-0.02
         ccf_out[*, i] = ccf[*,0]
         ;print, line[i],ccf_tau, ccf[*,0]
         ;pause
         ; 
      endif
   endfor
   xyouts, pos[2]-0.5,pos[3]-0.03,/norm, textoidl('\tau_{exp}=')+string(more_info.tau_obs,format='(i0)'),charsize=1
   xyouts, pos[2]-0.8,pos[3]-0.03,/norm, textoidl('logL5100=')+string(more_info.logL5100,format='(f0.2)'),charsize=1
   xyouts, pos[2]-0.8,pos[3]-0.05,/norm, textoidl('i_{psf}=')+string(more_info.imag,format='(f0.2)'),charsize=1

   !P.multi = 0

   splog, 'Finished object RM_ID=', RM_ID
end

;----------------------
pro rm_plot_LC, range=range, vw_norm=vw_norm, sncut=sncut, epoch_id=epoch_id, outtag=outtag

   forward_function rm_qsofit
   ; 
   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)
   ; keep only RM targets
   target=fibermap[0:848]

   ; set-up the epochs of the LC
   if n_elements(epoch_id) eq 0 then epoch_id=indgen(18)

   ; set up tags for the figfile and outfile
   if not keyword_set(outtag) then outtag='maxep_'+string(max(epoch_id+1),format='(i2.2)')

   ; get all the LC data
   rm_get_rms_prop,result=LC_data_all, /silent,epoch_id=epoch_id

   ind=sort(target.zfinal)
   if keyword_set(range) then ind=ind[range[0]:range[1]]
   nobj=n_elements(ind)

   ; read in more target information
   file='/data3/quasar/yshen/work/composite_lag/target_info.fits'
   info=mrdfits(file,1)

   figfile='/data3/quasar/yshen/work/composite_lag/LC_all_' + outtag + '.ps'
   if keyword_set(vw_norm) then figfile='/data3/quasar/yshen/work/composite_lag/LC_all_vw_norm_' + outtag + '.ps'
   ;message, 'stop'
   begplot,name=figfile,/color,/cmyk ; /encap
   !P.multi = [0,1,3]

   ; setup the structure to store all the CCF results
   dt = 2. ; days
   ccf_tau = -10. + (findgen(150./dt + 5L)*dt)
   ccf_struct = {RM_ID:-1L, ccf_tau:ccf_tau, ccf:dblarr(n_elements(ccf_tau), 5)}
   ccf_struct = replicate(ccf_struct, nobj)

   for i=0L, nobj - 1L do begin
      rm_ID = ind[i]
      tag = ', z='+string(target[rm_ID].zfinal,format='(f5.3)')
      rm_plot_LC_one, rm_ID, LC_data=LC_data_all[rm_ID], vw_norm=vw_norm, $
        ccf_tau=ccf_tau, ccf_out=ccf_out, tag=tag, more_info=info[rm_ID],epoch_id=epoch_id,sncut=sncut
      ccf_struct[i].RM_ID = rm_ID
      ccf_struct[i].ccf = ccf_out
   endfor

   !P.multi=0
   endplot

   outfile='/data3/quasar/yshen/work/composite_lag/ccf_data_' + outtag + '.fits'
   if keyword_set(vw_norm) then outfile='/data3/quasar/yshen/work/composite_lag/ccf_data_vw_norm_' + outtag + '.fits'
   mwrfits, ccf_struct, outfile, /create

end
