;+
; NAME:
;   extract_image
;
; PURPOSE:
;   Extract the fiber profile flux for an entire image
;
; CALLING SEQUENCE:
;   extract_image(fimage, invvar, xcen, sigma, flux, [error,
;              ymodel=ymodel, fscat=fscat, proftype = proftype, 
;              wfixed=wfixed, sigmacor=sigmacor, xcencor=xcencor, mask=mask,
;              nPoly=nPoly, maxIter=maxIter, highrej=highrej, lowrej=lowrej])
;
; INPUTS:
;   fimage     - Image[nCol, nRow]
;   invvar     - Inverse Variance[nCol, nRow]
;   xcen       - Initial guesses for X centers[nRow, nFibers]
;   sigma      - sigma of gaussian profile; default to 1.0 
;                  (scalar or [nFibers] or [nRow, nFibers])
;
; OUTPUTS:
;   flux       - Total extracted flux in each profile [nRow, nFibers]
;
; OPTIONAL OUTPUTS:
;   error      - Estimated total error in each profile [nRow, nFibers]
;   ymodel     - model best fit of row[nCol, nRow]
;   fscat      - scattered light contribution in each fiber[nRow, nFibers]
;
; OPTIONAL KEYWORDS:
;   proftype   - currently, one can only use 1: Gaussian (scalar)
;   wfixed     - array of 1's and zero's which set which parameters are fixed.
;                e.g. [1] just gaussian's with fixed width sigma
;                     [1, 1] fit gaussian + sigma correction
;                     [1, 0, 1] fit gaussian + center correction
;                     [1, 1, 1] fit gaussian + sigma and center corrections.   
; 
;   sigmacor   - new estimates of sigma, must have second element of wfixed set
;   xcencor    - new estimates of xcen, must have third element of wfixed set
;
;   mask       - byte mask: 1 is good and 0 is bad [nCol,nRow] 
;	         mask can be passed as input, but will be modified as
;                pixels are rejected
;   nPoly      - order of chebyshev scattered light background; default to 5
;   maxIter    - maximum number of profile fitting iterations; default to 10
;   highrej    - positive sigma deviation to be rejected (default 5.0)
;   lowrej     - negative sigma deviation to be rejected (default 5.0)
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;    extract_row.pro
;
; REVISION HISTORY:
;    8-Aug-1999  Version 0.0 Scott Burles, Chicago 
;-
;------------------------------------------------------------------------------
subroutine extract_row, fimage, invvar, xcen, sigma, flux, error, 
               ymodel=ymodel, fscat=fscat, proftype = proftype, 
               wfixed=wfixed, sigmacor=sigmacor, xcencor=xcencor, mask=mask,
               nPoly=nPoly, maxIter=maxIter, highrej=highrej, lowrej=lowrej

   ; Need 5 parameters
   if (N_params() LT 5) then begin
      print, 'Syntax - extract_image(fimage, invvar, xcen, sigma, flux, [error,'
      print, ' ymodel=ymodel, fscat=fscat, proftype = proftype, '
      print, ' wfixed=wfixed, sigmacor=sigmacor, xcencor=xcencor, mask=mask, '
      print, ' nPoly=nPoly, maxIter=maxIter, highrej=highrej, lowrej=lowrej])'
      return
   endif

   nTrace = n_elements(xcen)
   nx = n_elements(fimage)

   if (n_elements(sigma) NE nTrace) then begin
    sigma1 = sigma[0]
    sigma = xcen*0.0 + sigma1
   endif 

   if (NOT keyword_set(nPoly)) then nPoly = 5
   if (NOT keyword_set(maxIter)) then maxIter = 10
   if (NOT keyword_set(highrej)) then highrej = 5.0
   if (NOT keyword_set(lowrej)) then lowrej = 5.0
   if (NOT keyword_set(wfixed)) then wfixed = [1]

   xvar = findgen(nx)

   if (NOT keyword_set(mask)) then mask = lonarr(nx) + 1 $
      else if (nx NE n_elements(mask)) then $
         message, 'Number of elements in FIMAGE and MASK must be equal'

   nCoeff = n_elements(wfixed)       ;Number of parameters per fibers
   proftype = 1                      ;Gaussian

   if (nx NE N_elements(invvar)) then $
    message, 'Number of elements in FIMAGE and INVVAR must be equal'

;
;	Check that xcen is sorted in increasing order
;	with separations of at 3 pixels.
;
   check = where(xcen[0:nTrace-1] GE xcen[1:nTrace-2] - 3,count)

   if (count GT 0) then $
      message, 'XCEN is not sorted or not separated by greater than 3 pixels.'


   nPoly = LONG(nPoly)
   ma = nPoly + nTrace*nCoeff
   maxIter = LONG(maxIter)
   proftype = LONG(proftype)
   calcCovar = LONG(calcCovar)
		
   p = fltarr(ma)         ; diagonal errors
   ymodel = fltarr(nx)
   fscat = fltarr(nTrace)


   if (NOT keyword_set(wfixarr)) then begin
      wfixarr = lonarr(ma) + 1 	       ; Fixed parameter array
      for i=0,nCoeff-1 do wfixarr(lindgen(nTrace)*nCoeff+i) = wfixed[i]
      if (keyword_set(bfixarr)) then wfixarr(nTrace*nCoeff:ma-1) = bfixarr
   endif else if (ma NE n_elements(wfixarr)) then $
      message, 'Number of elements in FIMAGE and WFIXARR must be equal'

   ans = fltarr(ma)       ; parameter values

   if (keyword_set(iback)) then begin
      if (nPoly NE n_elements(iback)) then $
         message, 'Number of elements in IBACK is not equal to nPoly'
      ans(nTrace*nCoeff:ma-1) = iback 
   endif

   if (keyword_set(inputans)) then begin
      if (nTrace*nCoeff NE n_elements(inputans)) then $
         message, 'Number of elements in INPUTANS is not equal to nTrace*nCoeff'
      ans(0:nTrace*nCoeff-1) = inputans
   endif

   p = fltarr(ma)         ; diagonal errors
   covar = fltarr(ma,ma)  ; full covariance matrix

   finished = 0
   niter = 0
   totalreject = 0

   while(finished NE 1) do begin 

      workinvvar = invvar*mask

      result = call_external(getenv('IDL_EVIL')+'libspec2d.so','extract_row',$
       nx, float(xvar), float(fimage), float(invvar), float(ymodel), nTrace, $
       nPoly, float(xcen), float(sigma), proftype, calcCovar, nCoeff, ma, ans, $
       long(wfixarr), p, fscat, covar)

       diffs = (fimage - ymodel)*sqrt(workinvvar) 
       chisq = total(diffs*diffs)
       reducedChi = chisq/(total(mask) - ma)
       scaleError = 1.0
       if (reducedChi GT 1.0) then scaleError = sqrt(reducedChi)

       badhigh = where(diffs GT highrej*scaleError, badhighct)
       badlow = where(diffs LT -lowrej*scaleError, badlowct)

       finished = 1
       if(badhighct GT 0) then begin
          mask(badhigh) = 0
          totalreject = totalreject + badhighct
          finished = 0
       endif 
       if(badlowct GT 0) then begin
          mask(badlow) = 0
          totalreject = totalreject + badlowct
          finished = 0
       endif

       niter = niter + 1
       if (niter EQ maxIter) then finished = 1

   endwhile

   oback = ans[ma-nPoly:ma-1]
   return, reform(ans[0:nTrace*nCoeff-1],nCoeff,nTrace)
end
;------------------------------------------------------------------------------
