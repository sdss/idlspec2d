;+
; NAME:
;   func_fit
;
; PURPOSE:
;   Convert from an array of x,y positions to a trace set
;
; CALLING SEQUENCE:
;   res = function func_fit( x, y, ncoeff, [invvar=invvar, function_name=function_name, $
;    /halfintwo, yfit=yfit, inputans=inputans, ia=ia]
;
; INPUTS:
;   x          - X values (independent variable)
;   y          - Y values (dependent variable)
;   ncoeff     - Number of coefficients to fit
;   invvar     - (INVERSE VARIANCE) Weight values; only supports 0<= for masked, >0 for unmasked
;
; OPTIONAL KEYWORDS:
;   function_name - Function to fit; options are:
;                'legendre'
;                'chebyshev'
;                Default to 'legendre'
;   halfintwo  - Optional keyword for Chebyshev function, FCHEBYSHEV
;   inputans   - if holding parameters fixed, set this array to values desired
;   ia         - array of 1's and 0's specifying free (1's) and fixed (0's) variables
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

function func_fit, x, y, ncoeff, invvar=invvar, function_name=function_name, $
 halfintwo=halfintwo, yfit=yfit, inputans=inputans, ia=ia

   if (N_params() LT 3) then begin
     print,'function func_fit, x, y, ncoeff, function_name=function_name' 
     print,'  halfintwo=halfintwo, yfit=yfit, invvar=invvar' 
     print,'function_name can be legendre or chebyshev'
     return, 0
   endif

   if (NOT keyword_set(function_name)) then function_name = 'flegendre'
   nx = N_elements(x)
   ny = N_elements(y)
   if (nx NE ny) then message, 'Dimensions of X and Y do not agree'

   if (keyword_set(invvar)) then begin
      nw = N_elements(invvar)
      if (nx NE nw) then message, 'Dimensions of X and W do not agree'
   endif else begin
      invvar = intarr(nx) + 1
   endelse

   if (NOT keyword_set(ia)) then ia = bytarr(ncoeff) + 1

    goodia = ia NE 0

   ; Select unmasked points
   igood = where(invvar GT 0, ngood)
   if (ngood EQ 0) then $
    message, 'No good data points in fit'

   res = fltarr(ncoeff)

   fixed = where(goodia EQ 0, nfix)
   if (nfix GT 0) then $
     if (NOT keyword_set(inputans)) then inputans = fltarr(ncoeff)

   if (ngood EQ 1) then begin

      res[0] = y[igood[0]]
      yfit = y[igood[0]]

   endif else begin

      ; Do not fit more coefficients than there are data points
      ncfit = ncoeff < nx

      if (function_name EQ 'flegendre') then $
       legarr = flegendre(x, ncfit)
      if (function_name EQ 'fchebyshev') then $
       legarr = fchebyshev(x, ncfit, halfintwo=halfintwo)

;
;	Subtract fixed terms first
;
     
      yfix = y * 0.0 
      nonfix = where(goodia[0:ncfit-1], nparams)
      fixed = where(goodia[0:ncfit-1] EQ 0, nfix)
      
      if (nfix GT 0) then begin
         yfix = legarr # (inputans * (1 - goodia)) 
         ysub = y - yfix
         finalarr = legarr[*, nonfix]
      endif else begin
         finalarr = legarr
         ysub = y
      endelse

      extra2 = finalarr * ((fltarr(nparams) + 1) ## (invvar > 0))
      beta = transpose( (ysub * (invvar > 0)) # finalarr)
      alpha = transpose(finalarr) # extra2
;
;	Don't send just one parameter to svdfit
;
      if (nparams GT 1) then begin
        svdc, alpha, w, u, v, /double
        res[nonfix] = svsol(u, w, v, beta, /double)
      endif else begin
        res[0] = beta / alpha
      endelse

      if (nfix GT 0) then res[fixed] = inputans[fixed]

      yfit = legarr # res 

   endelse

   return, res
end
