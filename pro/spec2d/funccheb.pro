;+
; NAME:
;   funccheb
;
; PURPOSE:
;   Evaluate chebyshev polynomial at x positions
;
; CALLING SEQUENCE:
;    y = funccheb( x, d)
;
; INPUTS:
;   x	       - Input vector (1D) of desired x positions
;   d          - chebyshev coefficients
;
; OUTPUTS:
;   y          - Evaluation of chebyshev polynomials
;

function funccheb,x,d

        NL=N_ELEMENTS(D)
        NX=N_ELEMENTS(X)
        b2 = fltarr(NX)
        b0 = fltarr(NX)
        b1 = fltarr(NX)
        twox = 2.0*x
        FOR I=0,NL-1 DO BEGIN
          b2 = b1
          b1 = b0
          ni = NL - i - 1
          b0 = twox[*] * b1[*] - b2 + d[ni]
        ENDFOR
        T = 0.5*(b0 - b2 + d[0])
        RETURN,T
end

	
