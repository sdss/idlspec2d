;+
; NAME:
;   djs_polyfit
;
; PURPOSE:
;   Version of POLYFITW that uses SVD fitting
;
; CALLING SEQUENCE:
;   res = djs_polyfit(x, y, ndegree, [ variance=, yfit= ] )
;
; INPUTS:
;   x -
;   y -
;   ndegree - Degree of polynomial fit, where NDEGREE+1 terms are fit
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   res - Polynominal coefficients [NDEGREE+1]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Perform polynomial fits using SVD.  If the fit is singular, then
;   redo the fit with fewer terms until it is not singular.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   16-Feb-2011  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
function djs_polyfit, x, y, ndegree, variance=variance1, yfit=yfit

   npix = n_elements(x)
   if (n_elements(y) NE npix) then $
    message, 'Number of elements in X and Y must agree!'
   if (ndegree LT 0) then $
    message, 'NDEGREE must be non-negative'
   if (keyword_set(variance1)) then begin
      if (n_elements(variance1) NE npix) then $
       message, 'Number of elements in X and VARIANCE must agree!'
      variance = double(reform(variance1, npix))
   endif else begin
      variance = 0
   endelse

   nfit = fix(ndegree) + 1
   status = 1
   while (nfit GT 0 AND status NE 0) do begin
      res = svdfit(reform(double(x),npix), reform(double(y),npix), $
       nfit, variance=variance, yfit=yfit, status=status, singular=singular, $
       sing_values=sing_values)
      nfit -= 1
   endwhile

   ; If some degrees of freedom were dropped, then add them back in as zeros
   ndrop = fix(ndegree) + 1 - n_elements(res)
   if (ndrop GT 0) then res = [res, fltarr(ndrop)]

   return, res
end
;------------------------------------------------------------------------------
