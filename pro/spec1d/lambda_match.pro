;  27-Jun-2000 D. Finkbeiner
; Force spectra on to same wavelength system 
;  28-Jun-2000 DF vectorized to handle arr being a 2-d array

PRO lambda_match, refwave, objwave, arr

; right now assume objwave is just 1d array

  nref  = n_elements(refwave)  ; number of samples in reference array
  nwave = (size(objwave))[1]
  sarr  = (size(arr))[1]   ; size of spectra
  narr  = n_elements(arr)/sarr   ; number of spectra

  IF sarr NE nwave THEN BEGIN
      print, 'LAMBDA_MATCH: Objwave and arr must have same dimensions'
      help, objwave, arr
      return
  ENDIF 

  rat = mean(objwave[*, 0]/refwave)
  shf = round(alog10(rat)*10000)

  IF shf GT 0 THEN BEGIN 
      print, 'LAMBDA_MATCH ', shf, ' pixel shift'
      arr = [fltarr(shf, narr), arr]
  ENDIF ELSE BEGIN 
      IF shf LT 0 THEN BEGIN 
          arr = arr[(-shf):sarr-1, *]
      ENDIF 
  ENDELSE 

  npad = nref-nwave-shf   ; how many more to add at end

  IF npad GT 0 THEN BEGIN 
      arr = [arr, fltarr(npad, narr)]
  ENDIF 
  IF npad LT 0 THEN BEGIN    ; remove if too many
      arr = arr[0:nref-1, *]
  ENDIF 

  return
END
