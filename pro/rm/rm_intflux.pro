;+
; NAME:
;   rm_intflux
;
; PURPOSE:
;   measure integrated line flux from the spectrum; fit a straight line to the pixels just
;   blueward and redward of the line to remove the continuum. The line flux error is estimated
;   using MC trials.
;
; CALLING SEQUENCE:
;   result=rm_intflux(wave,flux,ivar,int_range=[4400., 4700.],/qa,err=)
;
; INPUTS:
;   wave    -  observed wavelength
;   flux    -  flux density
;   ivar    -  inverse variance
;   err     -  flux error [optional]
;   int_range - range of the line
;   dlam_conti - span of continuum wave blue- and redward of the line to estimate the continuum
;
; OUTPUTS:
;   integrated lineflux and lineflux_err from MC trials.
;

function rm_intflux,wave1,flux1,ivar1,err1=err1,ntrial=ntrial,int_range=int_range,$
    qa=qa,mask=mask,dlam_conti=dlam_conti

   ; No of MC trials to compute line flux err
   if ~keyword_set(ntrial) then ntrial=50L

   ; range of integrated line flux
   if ~keyword_set(int_range) then int_range=[4900.,5100.]

   wave=wave1 & flux=flux1
   if keyword_set(err1) then begin
      err=err1
      ivar=err*0
      ind=where(err gt 1d-5)
      ivar[ind]=1./err[ind]^2
   endif else begin
      ivar=ivar1
      err=ivar*0
      ind=where(ivar gt 1d-5)
      err[ind]=1./sqrt(ivar[ind])
   endelse

   ; interpolate over bad pixels
   if ~keyword_set(mask) then $
      maskind = where(ivar le 0, complement=goodind) else begin
      maskind = mask ne 1
      goodind = mask ne 0
   endelse
   if maskind[0] ne -1 then begin
      flux[maskind] = interpol(flux[goodind],wave[goodind],wave[maskind])
   endif

   ; construct a local linear continuum just blueward and redward dlam_conti A of the line
   if ~keyword_set(dlam_conti) then dlam_conti = 25. 
   ind1 = where(wave ge int_range[0] - dlam_conti and $
                wave le int_range[1] + dlam_conti, npix)
   wave=wave[ind1] & flux=flux[ind1] & ivar=ivar[ind1] & err=err[ind1]
   ind_line=where(wave ge int_range[0] and wave le int_range[1], complement=ind_conti)
   ; fit a straight line to the continuum pixels
   result=linfit(wave[ind_conti],flux[ind_conti], measure_errors=err[ind_conti])
   conti_spec=flux
   conti_spec[ind_line]=result[0]+result[1]*wave[ind_line]
   line_spec = flux - conti_spec
   lineflux = int_tabulated(wave[ind_line],line_spec[ind_line],/double)

   ; estimate lineflux_err using MC trials
   conti_spec_mc = fltarr(npix, ntrial)
   line_spec_mc = fltarr(npix, ntrial)
   lineflux_mc = dblarr(ntrial)
   for i=0L, ntrial-1 do begin
      flux_new = flux + randomn(seed, npix)*err
      result=linfit(wave[ind_conti],flux_new[ind_conti], measure_errors=err[ind_conti])
      conti_spec_mc[*,i]=flux_new
      conti_spec_mc[ind_line,i]=result[0]+result[1]*wave[ind_line]
      line_spec_mc[*,i]=flux_new - conti_spec_mc[*,i]
      lineflux_mc[i]=int_tabulated(wave[ind_line],line_spec_mc[ind_line,i],/double)
   endfor
   lineflux_err=(quantile_1d(0.84,lineflux_mc) - quantile_1d(0.16,lineflux_mc) ) $
          * 0.5

   ; make a QA plot
   if keyword_set(qa) then begin
      yrange=[min(flux), max(flux) ]
      xrange=[int_range[0] - dlam_conti,int_range[1] + dlam_conti]
      plot, wave1, flux1, yrange=yrange,xrange=xrange
      oplot, wave, flux, color=cgcolor('red')
      oplot, wave, conti_spec, color=cgcolor('blue')
      oplot, [int_range[0], int_range[0]], [-1,1],thick=3,line=2,color=cgcolor('green')
      oplot, [int_range[1], int_range[1]], [-1,1],thick=3,line=2,color=cgcolor('green')
   endif

   return, [lineflux, lineflux_err]

end
