;called by veldisp

; this routine:
;   1) interpolates accross missing pixels
;   2) normalizes spectrum to unit flux (or variance)
;   3) apodizes spectrum with cos bell in real space
;   4) pads with zeros up to TWICE the next higher value of 2^N
;   5) computes FFT

PRO veldisp_fft, flux, err, npixbig, fluxfft, fluxfilt, fluxvar0, fluxvariancefft, wave=wave, keep=keep, klo_cut=klo_cut, khi_cut=khi_cut

  npixobj = n_elements(flux)
; interpolate over bad regions
  flux = djs_maskinterp(flux, err LE 0.0, /const)
  
; Normalize the object flux to be near unity
  norm = djs_mean(flux) > djsig(flux)

  IF norm GT 0 THEN BEGIN 
      flux = flux / norm
      err = err / norm
  ENDIF 

  IF keyword_set(wave) AND keyword_set(keep) THEN BEGIN 
      kmask = (wave GE min(keep)) AND (wave LE max(keep))
      err = err*kmask

  ENDIF 
 
; apodize
  fft_apodize, flux, err
  
; pad
  flux = [flux, fltarr(npixbig-npixobj)]
  err = [err, fltarr(npixbig-npixobj)]
  
; take FFT
  fluxfft = fft(flux) * npixbig
  fluxvariancefft = fft(err^2)  * npixbig
  fluxvar0 = float(fluxvariancefft[0])
 

; Band-pass filter the object spectrum
  fluxfilt = bandpassfilter(fluxfft, klo_cut=klo_cut, khi_cut=khi_cut)
  IF total(finite(fluxfilt) EQ 0) NE 0 THEN BEGIN 
      message, 'Infinite value in FFT'
  ENDIF 
  

  return
END
