;+
; PURPOSE:
;   perform lag measurement with input light curves
;
;

pro rm_lag_diag,t_conti,f_conti,e_conti,t_line,f_line,e_line,dt=dt,tag=tag,$
     outdir=outdir,ccf_xrange=ccf_xrange, fig_append=fig_append, $
     tau_exp=tau_exp,dcf=dcf,ccf_tau=ccf_tau,lag_mea=lag_mea, peak_sig=peak_sig, $
     noplot=noplot, ntrial=ntrial,mirror=mirror,dosub=dosub, out_ccf=out_ccf

   ; resolution of the CCF
   if ~keyword_Set(dt) then dt=2.

   if n_elements(ntrial) eq 0 then ntrial = 500L

   if ~keyword_set(outdir) then outdir='/data3/quasar/yshen/work/lags/prepspec/'
   if ~keyword_set(tag) then tag='lag_diag'

   if (~keyword_Set(fig_append)) and (~keyword_set(noplot)) then begin ; start a new figure
      figfile=outdir+tag+'.ps'
      begplot,name=figfile,/color,/cmyk ; /encap
   endif
   !P.multi = [0,1,3]

   ; plot the lightcurves
   if ~keyword_set(noplot) then begin
     charsize=2.5 & ticklen=0.05
     pos = [0.11, 0.7, 0.95, 0.95]
     xrange=[min(t_conti),max(t_conti)]
     yrange=[min(f_conti),max(f_conti)]
     plot,[0],[0],xrange=xrange,yrange=yrange,xtitle='Time (Days)',ytitle='Relative Flux', $
      title=tag, pos=pos, charsize=charsize,xticklen=ticklen,xtickformat='(i0)'
     oploterror,t_conti,f_conti,e_conti,psym=symcat(6)
     xyouts, pos[2] - 0.2, pos[3] - 0.04, 'continuum',/norm

     pos = [0.11,0.39, 0.95,0.64]
     xrange=[min(t_line),max(t_line)]
     yrange=[min(f_line),max(f_line)]
     plot,[0],[0],xrange=xrange,yrange=yrange,xtitle='Time (Days)',ytitle='Relative Flux', $
      pos=pos,charsize=charsize,xticklen=ticklen,xtickformat='(i0)'
     oploterror,t_line,f_line,e_line,psym=symcat(6)
     xyouts, pos[2] - 0.2, pos[3] - 0.04, 'line',/norm
   endif

   ; compute the CCF
   if ~keyword_set(ccf_xrange) then ccf_xrange=[-20,100]
   yrange=[-0.5,1]

   if ~keyword_set(ccf_tau) then $ 
    ; ccf_tau = -70. + (findgen(140./dt + 5L)*dt)
     ccf_tau = -10. + (findgen(70./dt + 5L)*dt)

   ; run xcorr_interp once to get the cent_tau and peak sig
   ccf0 = xcorr_interp(t_line,f_line,t_conti,f_conti,tau=ccf_tau,ndata=ndata0 $
            , /nogap,cent_tau=cent_tau0, peak_tau=peak_tau0 $
            , mirror=mirror,peak_sig=peak_sig)
   ; write to output CCF
   out_ccf = ccf0

   if keyword_set(dosub) then begin
      if peak_sig gt 0.999 and cent_tau0 gt 0 then doflag=1 else $
         doflag=0
   endif else doflag=1

   if ntrial gt 0 and doflag eq 1 then begin ; calculate MC error in lag
     ccf = xcorr_interp(t_line,f_line,t_conti,f_conti,tau=ccf_tau,ndata=ndata $
            , /nogap,err1=e_line,err2=e_conti, /bootstrap,cent_tau=cent_tau, peak_tau=peak_tau $
            , mirror=mirror,ntrial=500L)
     tau_mea = [cent_tau[0], quantile_1d(0.16, cent_tau), quantile_1d(0.84, cent_tau) ]
   endif else tau_mea = [cent_tau0, -1., -1.]
 
   if ~keyword_set(noplot) then begin
     pos = [0.11,0.08,0.95,0.33]
     plot,[0],[0],xrange=ccf_xrange,yrange=yrange,xtitle='Time Delay (Days)',ytitle='CCF', $
      pos=pos,charsize=charsize,xticklen=ticklen
     ypos=pos[3] - 0.03
     oplot, ccf_tau, ccf[*,0] ; ,color=colors
     oplot, [tau_mea[0],tau_mea[0]],yrange,color=cgcolor('blue')
     oplot, [tau_mea[1],tau_mea[1]],yrange,color=cgcolor('blue'),line=2
     oplot, [tau_mea[2],tau_mea[2]],yrange,color=cgcolor('blue'),line=2

     ; overplot the expected lag
     oplot, [tau_exp,tau_exp],yrange,color=cgcolor('red'),line=1
     legend, textoidl('\tau_{exp}'),pos=[pos[2]-0.2,pos[3]-0.02],$
       line=1,color=cgcolor('red'),/norm,box=0
   endif

   print, tau_mea, '[tau0, peak sig]=', cent_tau0, peak_sig
   lag_mea=tau_mea

   ; compute the discrete CCF
   if keyword_set(dcf) then begin
     rm_dcf,t_conti,f_conti,t_line,f_line,err1=median(e_conti),err2=median(e_line),$
              DCF=DCF,eDCF=eDCF,tau_grid=tau_grid
     ind1 = where(eDCF gt 0)
     if ~keyword_set(noplot) then begin
       if ind1[0] ne -1 then begin
         oploterror, tau_grid[ind1],dcf[ind1],edcf[ind1],psym=symcat(6),color=cgcolor('red'),$
          errcolor=cgcolor('red')
       endif
     endif
   endif
    
   !P.multi = 0


   if (~keyword_set(fig_append)) and (~keyword_set(noplot)) then endplot

end
