; D. Finkbeiner
; Prime factor routine for calculating FFT efficiency. 
; totally bonehead.
PRO primfac, n_in

  n = float(n_in)
  FOR i=2L, n+1 DO BEGIN 
      WHILE n/i-long(n/i) EQ 0.0 DO BEGIN 
          n = n/i
          print, i, format='(I6,$)'
      ENDWHILE 

  ENDFOR 
  print
  return
END

