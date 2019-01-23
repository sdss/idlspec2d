;+
; NAME:
;   rm_vw_scaling
;
; PURPOSE:
;   Perform relative (to the reference spectrum) flux scaling using 
;   the Van Groningen & Wanders (1992) method; assuming [OIII]5007 flux is constant
;   
;
;----------------------------------------------
function chi2_vw, pp, order=order,diag=diag,wave_range=wave_range
; pp = [fsca, delta_lam, delta_fwhm]
   common vw_diff_spec,wave1_c,wave2_c,flux1_c,err1_c,flux2_c,err2_c

   if not keyword_set(wave_range) then wave_range=[4950, 5050]
   
   wave1=wave1_c & wave2=wave2_c
   flux1=flux1_c & err1=err1_c
   flux2=flux2_c & err2=err2_c
   
   ; shift the wavelength in units of SDSS pixel, and adjust the error array
   dlnlam=2.30259d-4
   wave2=wave2*(1. - pp[1]*dlnlam)
   aa=floor(pp[1]) & ff=(pp[1]-aa)
   err2=sqrt( (shift(err2,aa))^2*(1-ff)^2 + (shift(err2,aa+1))^2*ff^2 )
   
   ; convolve the input spectrum with delta_fwhm [in pixel]
   sig_pix=abs(pp[2])/2.3548
   khalfsz = round (4*sig_pix+1)
   xx= findgen(khalfsz*2+1) - khalfsz
   kernel = exp(-xx^2/(2*sig_pix^2))
   kernel = kernel/total(kernel)
   if pp[2] gt 0 then begin ; broaden spectrum2
      flux2=convol(flux2, kernel, /center, /edge_truncate)
      err2=sqrt(convol(err2^2,kernel^2,/center,/edge_truncate) )
   endif
   if pp[2] lt 0 then begin ; broaden spectrum1
      flux1=convol(flux1, kernel, /center, /edge_truncate)
      err1=sqrt(convol(err1^2,kernel^2,/center,/edge_truncate) )
   endif

   ; map the spectrum2 onto wave1
   flux2=interpol(flux2,wave2,wave1)
   err2=interpol(err2,wave2,wave1)

   ; now limit the wavelength range
   ind=where(wave1 ge wave_range[0] and wave1 le wave_range[1],npix1)
   wave1=wave1[ind] & flux1=flux1[ind] & err1=err1[ind]
   flux2=flux2[ind] & err2=err2[ind]

   ;generate the difference spectrum
   diff_flux = flux2*pp[0] - flux1
   diff_err = sqrt( (err1)^2 + (err2*pp[0])^2)

   ; fit a low-order polynomial to the difference spectrum
   result=poly_fit(wave1,diff_flux,order,measure_errors=diff_err,chisq=chi2, $
      /double,yfit=yfit)

   ;if keyword_set(diag) then begin
   ;  oplot, wave_c, diff_flux + mean(flux1_c)
   ;  pause
   ;endif

   return, chi2
end

pro rm_vw_scaling,wave_ref,flux_ref,err_ref,wave_input,flux_input,err_input, $
      vw_par=vw_par,fsca=fsca,delta_lam=delta_lam,delta_fwhm=delta_fwhm,diag=diag, $
      skymask=skymask

   ; default is not to plot diagnosis
   if n_elements(diag) eq 0 then diag=0

   wave1=wave_ref & flux1=flux_ref & err1=err_ref
   wave2=wave_input & flux2=flux_input & err2=err_input

   ; setup common block
   common vw_diff_spec,wave1_c,wave2_c,flux1_c,err1_c,flux2_c,err2_c
   wave1_c=wave1 & wave2_c=wave2
   flux1_c=flux1 & err1_c=err1
   flux2_c=flux2 & err2_c=err2

   ; do some diagnosis
   if keyword_set(diag) then begin
      plot, wave1_c, flux1_c
      oplot, wave2_c, flux2_c,color=cgcolor('red')
      pause
   endif

   ; find the parameter fsca that minimize the chi2 of a low-order polynomial 
   ; to the differece spectrum
   if not keyword_set(wave_range) then wave_range=[4950, 5050] ; default [OIII] region
   fcnargs = {order:2,diag:diag,wave_range:wave_range}
   start_value=[1.,0.,0.]
   parinfo=replicate({value:0.0, fixed:0L, limited:[0L,0L], limits:[0.,0.]}, 3)

   ; for SDSS, the spectral resolution and wavelength solution are constant, so don't fit them
   parinfo[1:2].fixed=1L
   parinfo.value=start_value

   bestfit=tnmin('chi2_vw',FUNCTARGS=fcnargs,parinfo=parinfo, $
      BESTMIN=bestmin,/autoderivative,/quiet)
   vw_par=bestfit ; [fsca, delta_lam, delta_fwhm]
   fsca=bestfit[0]  ; the scaling factor for the input spectrum

   ;print, vw_par

   if keyword_set(diag) then begin
      plot, wave1_c, flux1_c
      oplot, wave2_c, flux2_c*fsca,color=cgcolor('cyan')
      oplot, wave1_c, flux2_c*fsca-flux1_c+mean(flux1_c),color=cgcolor('blue'),linestyle=2
      pause
   endif

end
;--------------------------------------------
function rm_vw_fluxing_obj,rm_ID,ref_ep=ref_ep,diag=diag,skymask=skymask, $
 epoch_id=epoch_id

   ; The reference epoch; default is to use the first epoch in epoch_id list
   ; make sure this is a good epoch
   if n_elements(ref_ep) eq 0 then ref_ep=0 

   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)
   ; keep only RM targets
   target=fibermap[rm_ID]

   ; setup the epochs of the LC
   if not keyword_set(epoch_id) then epoch_id=indgen(15)
   nepoch=n_elements(epoch_id)

   plate=(target.plate)[epoch_id]
   fiber=(target.fiberid)[epoch_id]
   mjd=(target.mjd)[epoch_id]
   zz=target.zfinal   

   calibdir='recalib/'
   rm_readspec,plate[ref_ep],fiber[ref_ep],mjd=mjd[ref_ep],calibdir=calibdir, $
     wave=wave_ref,flux=flux_ref,flerr=err_ref
   ; shift to restframe
   wave_ref=wave_ref/(1.+zz)
   fsca_all = dblarr(nepoch)+1.

   ; estimate the scaling factor for all epochs
   for i=0L, nepoch-1 do begin

      rm_readspec,plate[i],fiber[i],mjd=mjd[i],calibdir=calibdir, $
        wave=wave_input,flux=flux_input,flerr=err_input
      wave_input=wave_input/(1.+zz)

      rm_vw_scaling,wave_ref,flux_ref,err_ref,wave_input,flux_input,err_input, $
        fsca=fsca,diag=diag,skymask=skymask

      fsca_all[i] = fsca

   endfor

   return, fsca_all
end



