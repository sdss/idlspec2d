; DPF 29 Mar 2000
; Compute gaussian with periodic boundary conditions
FUNCTION gauss_periodic, x, a, shft=shft

  amp = a[0]
  cen = a[1]
  sig = a[2]

  n = n_elements(x)
  IF (NOT keyword_set(shft)) THEN $  ; assume linear
    shft = ((max(x)-min(x))*n)/(n-1)
  
  g = 0
  FOR i=-3, 3 DO BEGIN 
      z = (x-cen-i*shft)/sig
      ez = exp(-z^2/2.d)
      g = g+amp*ez
  ENDFOR 

  return, g
END
