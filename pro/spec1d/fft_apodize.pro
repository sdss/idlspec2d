;------------------------------------------------------------------------------
; Apodize a spectrum with cos bell before doing FFT's.
; Use FLUXERR to select the first and last good points at which to
; start the apodizations.

pro fft_apodize, flux, fluxerr

   npix = N_elements(flux)

   if (keyword_set(fluxerr)) then begin
      igood = where(fluxerr GT 0, ngood)
      if (ngood LE 2) then return
      i1 = igood[0]
      i2 = igood[ngood-1]
      if (i1 GT 0) then flux[0:i1-1] = 0
      if (i2 LT npix-1) then flux[i2+1:npix-1] = 0
   endif else begin
      i1 = 0
      i2 = npix-1
   endelse

   objbell = 0.0 * flux
   objbell[i1:i2] = cosbell(i2-i1+1, 0.2, $
    double=(size(flux, /tname) EQ 'DOUBLE'))

   fluxmean = total(flux*objbell)/total(objbell)
   flux = (flux - fluxmean) * objbell
   if (keyword_set(fluxerr)) then $
    fluxerr = fluxerr * objbell

   return
end
;------------------------------------------------------------------------------
