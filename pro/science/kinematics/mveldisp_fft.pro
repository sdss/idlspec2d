;called by mveldisp

; this routine:
;   1) interpolates accross missing pixels
;   2) normalizes spectrum to unit flux (or variance)
;   3) apodizes spectrum with cos bell in real space
;   4) pads with zeros up to TWICE the next higher value of 2^N
;   5) computes FFT

PRO mveldisp_fft, flux_in, err_in, npixbig, fluxfft, fluxfilt, fluxvar0, fluxvariancefft, err, khicut, wave=wave, keep=keep, klo_cut=klo_cut, khi_cut=khi_cut

; make copies to work on
  flux = flux_in-mean(flux_in)
  err = err_in


  npixobj = n_elements(flux)
; interpolate over bad regions
  flux = djs_maskinterp(flux, err LE 0.0, /const)
  

  IF keyword_set(wave) AND keyword_set(keep) THEN BEGIN 
      kmask = (wave GE min(keep)) AND (wave LE max(keep))

;stop
      err = err*kmask

  ENDIF 
 
; apodize
  fft_apodize, flux, err
  
; pad
  npad = npixbig-npixobj
  IF npad NE 0 THEN BEGIN 
      flux = [flux, fltarr(npixbig-npixobj)]
      err = [err, fltarr(npixbig-npixobj)]
  ENDIF ELSE BEGIN 
      print, 'MVELDISP_FFT:  Warning: arrays already padded...'
  ENDELSE 
; take FFT
  fluxfft = fft(flux) * npixbig
  fluxvariancefft = fft(err^2)  * npixbig
  fluxvar0 = float(fluxvariancefft[0])


  w = where(err NE 0)

; Band-pass filter the object spectrum
  fluxfilt = mbandpassfilter(fluxfft, klo_cut=klo_cut, khi_cut=khi_cut)
  norm = djsig((float(fft(fluxfilt, /inv)))[w], sigrej=6)

;   window,0, retain=2		
;   plot,fluxfilt,ps=0
;stop	
	

  fluxfiltmy=fluxfilt	
  fluxfilt = fluxfilt/norm


   ndata = N_elements(fluxfft)
   knums = fft_wavenums(ndata)

   bin1=where(knums GT 0.4 AND knums LT 0.49)
   lfluxfft=alog10((fluxfft/norm)*conj(fluxfft/norm))
   meanlflux1=mean(lfluxfft[bin1])
   rmslflux1=stddev(lfluxfft[bin1])

   i=1	
   diff=0
   khigh=0.4	
   while (diff LT float(3*rmslflux1) AND khigh GT 0) do begin
    bin=where(knums GT 0.4-0.05*i AND knums LT 0.449-0.05*i)
    meanlflux=mean(lfluxfft[bin])
    diff= abs(float(meanlflux) - float(meanlflux1))
    khigh=0.425-0.05*i	
   i=i+1   
   endwhile
   if (khigh GT klo_cut) then khicut=khigh
   if (khigh LT klo_cut) then khicut=klo_cut+0.001		 

;   window,1,retain=2	
;   plot,knums,alog10((fluxfft/norm)*conj(fluxfft/norm)),ps=0	
;   djs_oplot,knums,alog10(fluxfilt*conj(fluxfilt)),ps=0,color='red'
;stop

  IF total(finite(fluxfilt) EQ 0) NE 0 THEN BEGIN 
      message, 'Infinite value in FFT'
  ENDIF 
  

  return
END
