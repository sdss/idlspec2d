function fchebyshev,x,m
;+
; NAME:
;        FCHEBYSHEV
; PURPOSE:
;       Compute the first M terms in a CHEBYSHEV polynomial expansion.  
; EXPLANATION:
;       Meant to be used as a supplied function to SVDFIT.
;
;       This procedure became partially obsolete in IDL V5.0 with the 
;       introduction of the /LEGENDRE keyword to SVDFIT and the associated 
;       SVDLEG function.    However, note that, unlike SVDLEG, FCHEBYSHEV works
;       on vector values of X.     
; CALLING SEQUENCE:
;       result = FLEGENDRE( X, M)
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
;       IDL> print, fCHEBYSHEV(x,3)   ==>   [1.00, 2.88, 11.9416]
;
;       This result can be checked by explicity computing the first 3 Legendre
;       terms, 1.0, x, 0.5*( 3*x^2 -1)
;
;       (2) Find the coefficients to an M term Legendre polynomial that gives
;               the best least-squares fit to a dataset (x,y)
;               IDL> coeff = SVDFIT( x,y,M,func='fCHEBYSHEV')
;       
;           The coefficients can then be supplied to the function POLYLEG to 
;               compute the best YFIT values for any X. 
; METHOD:
;       The recurrence relation for the Legendre polynomials is used to compute
;       each term.   Compare with the function FLEG in "Numerical Recipes"
;       by Press et al. (1992), p. 674
;
; REVISION HISTORY:
;       Written     Wayne Landsman    Hughes STX      April 1995                
;       Converted to IDL V5.0   W. Landsman   September 1997
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
 
 leg[0,0] = replicate( 1., n)
 if m GE 2 then leg[0,1] = x
 for j=2,m-1 do begin
     leg[0,j] = 2.0 * x * leg[*,j-1] - leg[*,j-2] 
 endfor

 return, leg
 end

