;+
; NAME:
;   djs_iter_sfit
;
; PURPOSE:
;   Surface-fitting code to tabulated data (with iterateive rejection).
;
; CALLING SEQUENCE:
;   acoeff = djs_sfit_iter(fval, xval, yval, degreex, degreey, $
;            sqivar=, yfit=, mask=, maxdev=, maxrej=, $
;            upper=, lower=, maxiter=, freeiter=, outmask=)
;
; INPUTS:
;   fval      - Function values at XVAL,YVAL.
;   xval      - X coordinate values
;   yval      - Y coordinate values
;   degreex   - Degree of polynomial fit in X; 1 for linear, 2 for quadratic
;   degreey   - Degree of polynomial fit in Y; 1 for linear, 2 for quadratic
;
; OPTIONAL INPUTS:
;   sqivar    - Inverse sigma, which are the weights
;   mask      - Input mask for fitting and rejection (1=good, 0=bad)
;   maxdev    - Rejcet pts w/ abs(data - model) > maxdev (passed to djs_reject)
;   upper     - Rejcet pts w/ data > model + upper*sigma (passed to djs_reject)
;   lower     - Rejcet pts w/ data > model + lower*sigma (passed to djs_reject)
;   maxiter   - Maximum # of rejection iterations 
;   freeiter  - Rejected points accumulate for this # of iterations,
;               then the mask is reset to the input mask
;
; OUTPUTS:
;   acoeff    - Fit coefficients as [DEGREEX+1,DEGREEY+1] array
;
; OPTIONAL OUTPUTS:
;   yfit      - Fit values
;   outmask   - Masked points used in final iteration (1=used, 0=rejected)
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Doesn't like NaNs
;
; PROCEDURES CALLED:
;   djs_reject
;   splog
;
; REVISION HISTORY:
;   12-Aug-2003  Adapted from djs_sfit by C. Tremonti, Steward Observatory
;-
;------------------------------------------------------------------------------
function djs_sfit_iter, fval, xval, yval, degreex, degreey, $
   sqivar=sqivar, yfit=yfit, mask = mask, maxdev = maxdev, maxrej = maxrej, $
   upper=upper, lower=lower, maxiter = maxiter, freeiter = freeiter, $
   outmask = outmask

   npts = n_elements(fval)
   if (n_params() LT 5) then $
    message, 'Wrong number of parameters'
   if (n_elements(xval) NE npts OR n_elements(yval) NE npts) then $
    message, 'Dimensions of FVAL,XVAL,YVAL must agree'
   if (degreex LT 0 OR degreey LT 0) then $
    message, 'DEGREEX,DEGREEY must be non-negative integers'

   if not keyword_set(mask) then mask = intarr(npts) + 1
   if not keyword_set(maxiter) then maxiter = 10

   ncoeff = (degreex+1) * (degreey+1)
   if (keyword_set(sqivar)) then thisvec = reform(double(sqivar),npts) $
    else thisvec = dblarr(npts)+1.d0
   if (n_elements(thisvec) NE npts) then $
    message, 'Number of elements in SQIVAR must agree with FVAL'
   thisx = reform(double(xval), npts)
   thisy = reform(double(yval), npts)


   ;--------------------------------------------
   ; Solve matrix for good points only, then reject outliers and solve again

   iiter = 0
   inmask = mask
   outmask = mask
   while (NOT keyword_set(qdone) and iiter LE maxiter) do begin
   
     ok = where(inmask ne 0 and finite(fval), nok)

     mmatrix = dblarr(npts, ncoeff, /nozero)
     for xp=0, degreex do begin
        for yp=0, degreey do begin
           tmpvec = thisvec
           if (xp GT 0) then tmpvec = tmpvec * thisx^xp
           if (yp GT 0) then tmpvec = tmpvec * thisy^yp
           mmatrix[*,yp*(degreex+1)+xp] = tmpvec
        endfor
     endfor

     ; Only use good pts to compute coeff
     bvec = reform(double(fval),npts) * thisvec
     if (ncoeff EQ 1) then begin
        mt_m = total(mmatrix[ok,*] * mmatrix[ok,*])
        mt_b = total(mmatrix[ok,*] * bvec[ok,*])
        mmi = 1. / (mt_m + (mt_m EQ 0))
     endif else begin
        mt_m = matrix_multiply(mmatrix[ok,*], mmatrix[ok,*], /atranspose)
        mt_b = matrix_multiply(mmatrix[ok,*], bvec[ok], /atranspose)
        mmi = invert(mt_m, /double)
     endelse
     acoeff = mmi # mt_b
     zfit = acoeff ## mmatrix

     ; Reset masks after every few iterations 
     if keyword_set(freeiter) then begin
       if (iiter mod freeiter eq 0) then begin 
         outmask = mask  
         inmask = mask
       endif
     endif 

     ; Do the rejection
     qdone = djs_reject(bvec, zfit, inmask = inmask, outmask = outmask, $
             maxdev = maxdev, upper = upper, lower=lower, /sticky, $
             maxrej = maxrej)

     bad = where(outmask ne 1, nbad)
     if nbad gt 0 then inmask[bad] = 0
     iiter = iiter + 1
   endwhile
   splog, 'Rejected Points from 2d fit: ', nbad

   if (arg_present(yfit)) then begin
     mmatrix = mmatrix / rebin(thisvec, npts, ncoeff)  ; remove weights
     yfit = reform(acoeff ## mmatrix, size(yval, /dimens))
   endif

   return, reform(acoeff, degreex+1, degreey+1)
end
;------------------------------------------------------------------------------
