;+
; NAME:
;   lrgmodel_errors
;
; PURPOSE:
;   Make plots of photo-z's where we increase the errors.
;
; CALLING SEQUENCE:
;   lrgmodel_errors
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_oplot
;   djs_plot
;   hogg_scatterplot
;   lrgmodel_photoz()
;
; REVISION HISTORY:
;   22-Dec-2003  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro lrgmodel_errors

   ;----------
   ; Create a fake photo-z catalog

   foo = lrgmodel_photoz(zsplinearr=zarr, synfluxarr=synfluxarr)
   nobj = 5000
   iz = long(70 * findgen(nobj)/nobj)
   ztrue = zarr[iz]
   flux = synfluxarr[*,iz]

   ;----------
   ; Loop through different errors

   dfpsplot, 'photoz-errors.ps', /square, /color

   randarr = randomn(1234, 5, nobj)
   for magerr=0.01, 0.10, 0.01 do begin
      thiserr = magerr * flux / 1.086
      thisflux = flux + thiserr * randarr
      thisivar = 1. / thiserr^2
;thisivar[0,*] = 0 ; Set u-band S/N to zero
;thisivar[4,*] = 0 ; Set z-band S/N to zero

      zfit = lrgmodel_photoz(thisflux, thisivar)

      djs_plot, ztrue, zfit, $
       psym=3, xr=[0,0.8], yr=[0,0.8], $
       xtitle='True Z', ytitle='Photometric Z'
      djs_oplot, !x.crange, !x.crange, color='red'
      djs_xyouts, 0.1, 0.7, $
       'Errors = ' + string(magerr,format='(f4.2)') + ' mag'
   endfor

   dfpsclose

end
