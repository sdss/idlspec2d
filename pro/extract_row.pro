;+
; NAME:
;   extract_row
;
; PURPOSE:
;   Extract the flux with profile weighting at centroid positions.
;
; CALLING SEQUENCE:
;   ans = extract_row( fimage, invvar, xcen, [ymodel=ymodel,
;              fscat = fscat, proftype = proftype, wfixed = wfixed,
;              inputans=inputans, iback = iback, oback =oback, bfixarr=bfixarr,
;              diagonal=diagonal, fullcovar=fullcovar, wfixarr = wfixarr,
;              sigma=sigma, nPoly=nPoly, maxIter=maxIter, highrej=highrej, 
;              lowrej=lowrej])
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
;   nPoly      - order of chebyshev scattered light background; default to 5
;   maxIter    - maximum number of profile fitting iterations; default to 5
;   highrej    - positive sigma deviation to be rejected (default 5.0)
;   lowrej     - negative sigma deviation to be rejected (default 5.0)
;
; OUTPUTS:
;   ans        -  Extracted flux in each parameter [nCoeff, nFiber]
;
; OPTIONAL OUTPUTS:
;   ymodel     - model best fit of row[nCol]
;   fscat      - scattered light contribution in each fiber[nFibers]
;   diagonal   - full 1d diagonal of covariance matrix
;   fullcovar  - full 2d covariance matrix
;   wfixarr    - 1d integer array of 1's and zero's which specify fixed params
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   Dynamic link to extract_profile.c
;
; REVISION HISTORY:
;   24-Mar-1999  David Schlegel, Princeton.
;   24-Jun-1999  Stolen and modified by Scott Burles, Chicago.
;-
;------------------------------------------------------------------------------
function extract_row, fimage, invvar, xcen, sigma, ymodel=ymodel, $
                   fscat = fscat, proftype = proftype, wfixed = wfixed, $
                   inputans=inputans, iback = iback, oback =oback, $
                   bfixarr = bfixarr, $
                   diagonal=diagonal, fullcovar=fullcovar, wfixarr = wfixarr, $
                   nPoly=nPoly, maxIter=maxIter, highrej=highrej, $
                   lowrej=lowrej

   ; Need 4 parameters
   if (N_params() LT 4) then begin
      print, 'Syntax - ans = extract_row( fimage, invvar, xcen, [ymodel=ymodel,'
      print, ' fscat = fscat, proftype = proftype, wfixed = wfixed,'
      print, ' inputans=inputans, iback = iback, oback =oback, bfixarr=bfixarr,'
      print, ' diagonal=diagonal, fullcovar=fullcovar, wfixarr = wfixarr,'
      print, ' sigma=sigma, nPoly=nPoly, maxIter=maxIter, highrej=highrej, '
      print, ' lowrej=lowrej])'
      return, -1
   endif

   nTrace = n_elements(xcen)
   nx = n_elements(fimage)

   if ((size(sigma))[0] NE 1) then begin
    sigma1 = sigma[0]
    sigma = xcen*0.0 + sigma1
   endif


   if (NOT keyword_set(nPoly)) then nPoly = 5
   if (NOT keyword_set(maxIter)) then maxIter = 5
   if (NOT keyword_set(highrej)) then highrej = 5.0
   if (NOT keyword_set(lowrej)) then lowrej = 5.0
   if (NOT keyword_set(wfixed)) then wfixed = [1]

   nCoeff = n_elements(wfixed)       ;Number of parameters per fibers
   proftype = 1                      ;Gaussian

   if (nx NE N_elements(invvar)) then $
    message, 'Number of elements in FIMAGE and INVVAR must be equal'

;
;	Check that xcen is sorted in increasing order
;	with separations of at 3 pixels.
;
   check = where(xcen[0:nTrace-1] GE xcen[1:nTrace-2] - 3,count)

   if(count GT 0) then $
     message, 'XCEN is not sorted or not separated by greater than 3 pixels.'


   nPoly = LONG(nPoly)
   ma = nPoly + nTrace*nCoeff
   maxIter = LONG(maxIter)
		
   ymodel = fltarr(nx)
   x = findgen(nx)*2.0/nx - 1.0;
   fscat = fltarr(nTrace)

   ia = lonarr(ma) + 1 	       ; Fixed parameter array

   for i=0,nCoeff-1 do ia(lindgen(nTrace)*nCoeff+i) = wfixed[i]
   if (keyword_set(bfixarr)) then ia(nTrace*nCoeff:ma-1) = bfixarr

   ans = fltarr(ma)       ; paramter values

   if (keyword_set(iback)) then ans(nTrace*nCoeff:ma-1) = iback 
   if (keyword_set(inputans)) then ans(0:nTrace*nCoeff-1) = inputans

   p = fltarr(ma)         ; diagonal errors
   covar = fltarr(ma,ma)  ; full covariance matrix

   result = call_external(getenv('IDL_EVIL')+'libspec2d.so','extract_row',$
    nx, x, float(fimage), float(invvar), float(ymodel), nTrace, nPoly, $
    float(xcen), float(sigma), nCoeff, ma, ans, ia, p, covar, fscat)
    

   oback = ans[ma-nPoly:ma-1]
   return, rebin(ans[0:nTrace*nCoeff-1],nCoeff,nTrace)
end
;------------------------------------------------------------------------------
