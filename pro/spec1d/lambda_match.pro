;  27-Jun-2000 D. Finkbeiner
; Force spectra on to same wavelength system 

PRO lambda_match, refwave, objwave, arr

  nref  = n_elements(refwave)  ; number of samples in reference array
  nwave = n_elements(objwave)

  IF n_elements(arr) NE nwave THEN BEGIN
      print, 'LAMBDA_MATCH: Objwave and arr must have same dimensions'
      help, objwave, arr
      return
  ENDIF 

  rat = mean(objwave/refwave)
  shf = round(alog10(rat)*10000)

  IF shf GT 0 THEN BEGIN 
      print, 'LAMBDA_MATCH ', shf, ' pixel shift'
      arr = [fltarr(shf), arr]
  ENDIF 

  npad = nref-nwave-shf   ; how many more to add at end

  IF npad GT 0 THEN BEGIN 
      arr = [arr, fltarr(npad)]
  ENDIF 

  return
END
