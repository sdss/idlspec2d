function fchebyshev,x,m,halfintwo=halfintwo
;+
; NAME:
;        FCHEBYSHEV
; PURPOSE:
;       Compute the first M terms in a CHEBYSHEV polynomial expansion.  
; EXPLANATION:
;       Meant to be used as a supplied function to SVDFIT.
;
; CALLING SEQUENCE:
;       result = FCHEBYSHEV( X, M, /halfintwo)
;
; INPUTS:
;       X - the value of the independent variable, scalar or vector
;       M - number of term of the CHEBYSHEV expansion to compute, integer scalar 
;
; OUTPUTS:
;       result - (N,M) array, where N is the number of elements in X and M
;               is the order.   Contains the value of each CHEBYSHEV term for
;               each value of X
; EXAMPLE:
;       (1) If x = 2.88 and M = 3 then 
;       IDL> print, fchebyshev(x,3)   ==>   [1.00, 2.88, 15.5888]
;
;       (2) Find the coefficients to an M term Chebyshev polynomial that gives
;               the best least-squares fit to a dataset (x,y)
;               IDL> coeff = SVDFIT( x,y,M,func='fchebyshev')
;       
; METHOD:
;
; REVISION HISTORY:
;       04-Aug-1999  Written by Scott Burles by hacking FLEGENDRE code
;                    by Landsman in the Goddard libraries.
;-      
 On_Error,2

 if N_params() LT 2 then begin
        print,'Syntax - result = FCHEBYSHEV( x, m)
        return,0
 endif  

 if m LT 1 then message, $
        'ERROR - Order of CHEBYSHEV polynomial must be at least 1'
 N = N_elements(x)
 size_x = size(x)
 leg = make_array(n, m, type = size_x[size_x[0]+1] > 4)    

 skip = 0
 if keyword_set(halfintwo) then skip = 1

 if skip then begin
     leg[0,1] = replicate( 1., n)
     leg[0,0] = leg[0,1] * (x GE 0.0)
 endif else begin
     leg[0,0] = replicate( 1., n)
 endelse

 if m GE 2+skip then leg[0,1+skip] = x
 for j=2+skip,m-1 do begin
     leg[0,j] = 2.0 * x * leg[*,j-1] - leg[*,j-2] 
 endfor

 return, leg
 end

