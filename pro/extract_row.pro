;+
; NAME:
;   extract_row
;
; PURPOSE:
;   Fit the fiber profiles and background in a single row with least squares
;
; CALLING SEQUENCE:
;   ans = extract_row( fimage, invvar, xcen, sigma, [ymodel=ymodel,
;              fscat = fscat, proftype = proftype, wfixed = wfixed,
;              inputans=inputans, iback = iback, oback =oback, bfixarr=bfixarr,
;              xvar = xvar, mask=mask, relative=relative,
;              diagonal=diagonal, fullcovar=fullcovar, wfixarr = wfixarr,
;              nPoly=nPoly, maxIter=maxIter, highrej=highrej, niter=niter,
;              lowrej=lowrej, calcCovar=calcCovar, squashprofile=squashprofile])
;
; INPUTS:
;   fimage     - Image[nCol]
;   invvar     - Inverse Variance[nCol]
;   xcen       - Initial guesses for X centers[nFibers]
;   sigma      - sigma of gaussian profile; default to 1.0 (scalar or [nFibers])
;
; OPTIONAL KEYWORDS:
;   proftype   - currently, one can only use 1: Gaussian (scalar)
;   wfixed     - array of 1's and zero's which set which parameters are fixed.
;   inputans   - 2d array of input answers [nCoeff, nFibers]
;   iback      - 1d array of input background coeff 
;                    (needed if fixed parameters are non-zero)
;   bfixarr    - array of 1's and zero's which set which background 
;                    parameters are fixed.
;   wfixarr    - 1d integer array of 1's and zero's which specify fixed 
;                    profile parameters
;   xvar       - x values of fimage and invvar, default is findgen(nx) 
;   mask       - pixel mask of 1 is good and 0 is bad (nx) 
;   relative   - use reduced chisq to scale rejection threshold
;   nPoly      - order of chebyshev scattered light background; default to 5
;   maxIter    - maximum number of profile fitting iterations; default to 10
;   highrej    - positive sigma deviation to be rejected (default 5.0)
;   lowrej     - negative sigma deviation to be rejected (default 5.0)
;   calcCovar  - Calculate full covariance matrix (which is symmetric)
;                  and fill lower triangle of fullcovar.
;		   Default is 0, and should be left at 0
;                  unless fullcovar is being analyzed, adds a factor of 2 or
;		   3 to CPU time.
;   niter      - number of rejection iterations performed
;
; OUTPUTS:
;   ans        -  Extracted flux in each parameter [nCoeff, nFiber]
;
; OPTIONAL OUTPUTS:
;   ymodel     - model best fit of row[nCol]
;   fscat      - scattered light contribution in each fiber[nFibers]
;   diagonal   - full 1d diagonal of covariance matrix 
;		     (Currently, this is diagonal from cholesky decompostion,
;                     which is 1/error(j) ).
;   fullcovar  - full 2d covariance matrix
;
; COMMENTS:
;
;    Still need to do:
;       limits on chebyshev polynomial are assumed to be 0.0 <--> nx
;	these may need to be optional if only partial rows are being fit
;
;       Error codes need to be returned, currently no such codes are returned
;
;	
;
;	  
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   Dynamic link to extract_row.c
;
; REVISION HISTORY:
;    8-Aug-1999  Version 0.0 Scott Burles, Chicago 
;-
;------------------------------------------------------------------------------
function extract_row, fimage, invvar, xcen, sigma, ymodel=ymodel, $
                   fscat = fscat, proftype = proftype, wfixed = wfixed, $
                   inputans=inputans, iback = iback, oback =oback, $
                   bfixarr = bfixarr, xvar=xvar, mask=mask, $
		   relative=relative, squashprofile=squashprofile, $
                   diagonal=p, fullcovar=covar, wfixarr = wfixarr, $
                   nPoly=nPoly, maxIter=maxIter, highrej=highrej, $
                   lowrej=lowrej, calcCovar=calcCovar, niter=niter

   ; Need 4 parameters
   if (N_params() LT 4) then begin
      print, 'Syntax - ans = extract_row( fimage, invvar, xcen, sigma, [ymodel=ymodel,'
      print, ' fscat = fscat, proftype = proftype, wfixed = wfixed,'
      print, ' inputans=inputans, iback = iback, oback =oback, bfixarr=bfixarr,'
      print, ' xvar=xvar, mask=mask, relative=relative, '
      print, ' squashprofile=squashprofile, '
      print, ' diagonal=diagonal, fullcovar=fullcovar, wfixarr = wfixarr,'
      print, ' nPoly=nPoly, maxIter=maxIter, highrej=highrej, '
      print, ' lowrej=lowrej, calcCovar=calcCovar, niter=niter])'
      return, -1
   endif

   nTrace = n_elements(xcen)
   nx = n_elements(fimage)

   if (n_elements(sigma) NE nTrace) then begin
    sigma1 = sigma[0]
    sigma = xcen*0.0 + sigma1
   endif 

   if (NOT keyword_set(nPoly)) then nPoly = 5
   if (NOT keyword_set(maxIter)) then maxIter = 10
   if (NOT keyword_set(highrej)) then highrej = 15.0
   if (NOT keyword_set(lowrej)) then lowrej = 20.0
   if (NOT keyword_set(wfixed)) then wfixed = [1]
   if (NOT keyword_set(calcCovar)) then calcCovar = 0
   if (NOT keyword_set(proftype)) then proftype = 1
   relative = keyword_set(relative) 
   squashprofile =keyword_set(squashprofile) 

   if (NOT keyword_set(xvar)) then xvar = findgen(nx) $
      else if (nx NE n_elements(xvar)) then $
         message, 'Number of elements in FIMAGE and XVAR must be equal'

   if (NOT keyword_set(mask)) then mask = make_array(nx, /byte, value=1) $
      else if (nx NE n_elements(mask)) then $
         message, 'Number of elements in FIMAGE and MASK must be equal'

   nCoeff = n_elements(wfixed)       ;Number of parameters per fibers
;   proftype = 1                      ;Gaussian

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
   squashprofile = LONG(squashprofile)

		
   ymodel = fltarr(nx)
   fscat = fltarr(nTrace)


   if (NOT keyword_set(wfixarr)) then begin
      wfixarr = lonarr(ma) + 1 	       ; Fixed parameter array
      i=0
      wfixarr(lindgen(nTrace)*nCoeff+i) = wfixed[i] 
      for i=1,nCoeff-1 do $
	wfixarr(lindgen(nTrace)*nCoeff+i) = wfixed[i] * (1 - squashprofile)
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

      workinvvar = float(invvar*mask)

      result = call_external(getenv('IDL_EVIL')+'libspec2d.so','extract_row',$
       nx, float(xvar), float(fimage), workinvvar, float(ymodel), nTrace, $
       nPoly, float(xcen), float(sigma), proftype, calcCovar, squashprofile, $
       nCoeff, ma, ans, long(wfixarr), p, fscat, covar)

       diffs = (fimage - ymodel)*sqrt(workinvvar) 
       chisq = total(diffs*diffs)
	countthese = total(wfixarr)
       reducedChi = chisq/(total(mask) - countthese)
       scaleError = 1.0
       if (relative) then scaleError = sqrt(chisq/total(mask))

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
;       print,format='($,i)', niter
   endwhile

   oback = ans[ma-nPoly:ma-1]
   return, reform(ans[0:nTrace*nCoeff-1],nCoeff,nTrace)
end
;------------------------------------------------------------------------------
