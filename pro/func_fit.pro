;+
; NAME:
;   func_fit
;
; PURPOSE:
;   Convert from an array of x,y positions to a trace set
;
; CALLING SEQUENCE:
;   res = function func_fit( x, y, ncoeff, [ function_name=function_name, $
;    /halfintwo, yfit=yfit ]
;
; INPUTS:
;   x          - X values (independent variable)
;   y          - Y values (dependent variable)
;   ncoeff     - Number of coefficients to fit
;
; OPTIONAL KEYWORDS:
;   function_name - Function to fit; options are:
;                'legendre'
;                'chebyshev'
;                Default to 'legendre'
;   halfintwo  - Optional keyword for Chebyshev function
;
; OUTPUTS:
;   res        - Fit coefficients
;
; OPTIONAL OUTPUTS:
;   yfit       - Fit evaluated at the points X
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   fchebyshev()
;   flegendre()
;
; REVISION HISTORY:
;   10-Sep-1999  Written by D. Finkbeiner, APO?
;   16-Nov-1999  Modified by D. Schlegel to never fit more coefficients
;                than there are data points.
;-
;------------------------------------------------------------------------------

function func_fit, x, y, ncoeff, function_name=function_name, $
 halfintwo=halfintwo, yfit=yfit

   if (N_params() LT 3) then begin
     print,'function func_fit, x, y, ncoeff, function_name=function_name' 
     print,'function_name can be legendre or chebyshev'
     return, 0
   endif
   
   if (NOT keyword_set(function_name)) then function_name = 'flegendre'
   nx = N_elements(x)
   ny = N_elements(y)
   if (nx NE ny) then message, 'Dimensions of X and Y do not agree'

   res = fltarr(ncoeff)

   if (nx EQ 1) then begin

      res[0] = y[0]
      yfit = y[0]

   endif else begin

      ; Do not fit more coefficients than there are data points
      ncfit = ncoeff < nx

      if (function_name EQ 'flegendre') then $
       legarr = flegendre(x, ncfit)
      if (function_name EQ 'fchebyshev') then $
       legarr = fchebyshev(x, ncfit, halfintwo=halfintwo)
       
      beta = transpose(y # legarr)
      alpha = transpose(legarr) # legarr
      svdc, alpha, w, u, v, /double
      res[0:ncfit-1] = svsol(u, w, v, beta, /double)
      yfit = legarr # res

   endelse

   return, res
end
