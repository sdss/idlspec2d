;+
; NAME:
;   extract_row
;
; PURPOSE:
;   Fit the fiber profiles and background in a single row with least squares
;
; CALLING SEQUENCE:
;   ans = extract_row( fimage, invvar, xcen, sigma, [ymodel=, fscat=, 
;              proftype=, wfixed=, inputans=, iback=, bfixarr=, xvar=,
;              mask=, relative=, diagonal=, fullcovar=, wfixarr=, nPoly=,
;              maxIter=, highrej=, niter=, lowrej=, calcCovar=, squashprofile=,
;              whopping=, wsigma=, pixelmask=, reject= ])
;
; INPUTS:
;   fimage     - Image[nCol]
;   invvar     - Inverse Variance[nCol]
;   xcen       - Initial guesses for X centers[nFibers]
;   sigma      - sigma of gaussian profile; default to 1.0 (scalar or [nFibers])
;
; OPTIONAL KEYWORDS:
;   proftype   - one can use 1: Gaussian (scalar)
;                            2: exponential ^3
;			     3: exponential ^2.5
;   wfixed     - array of 1's and zero's which set which parameters are fixed.
;   inputans   - 2d array of input answers [nCoeff, nFibers]
;   iback      - 1d array of input background coeff 
;                    (needed if fixed parameters are non-zero)
;   bfixarr    - array of 1's and zero's which set which background 
;                    parameters are fixed.
;   wfixarr    - 1d integer array of 1's and zero's which specify fixed 
;                    profile parameters
;   xvar       - x values of fimage and invvar, default is findgen(nx) 
;   mask       - image mask: 1 is good and 0 is bad (nx) 
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
;   whopping   - traces with extra high flux need extra terms
;   wsigma     - sigma width of whopping profile (exponential, default 25)
;   pixelmask  - bits set due to extraction rejection [nFiber]
;   reject   - two elements array setting partial and full rejection thresholds 
;                  for profiles  (default [0.8, 0.2])  
;                    ---->used to be [0.8,0.4] when it was hardwired
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
   fscat = fscat, proftype = proftype, wfixed = wfixed, inputans=inputans, $
   iback = iback, bfixarr = bfixarr, xvar=xvar, mask=mask, relative=relative, $
   squashprofile=squashprofile, diagonal=p, fullcovar=fullcovar, $
   wfixarr = wfixarr, nPoly=nPoly, maxIter=maxIter, highrej=highrej, $ 
   lowrej=lowrej, calcCovar=calcCovar, niter=niter, reducedChi = reducedChi, $
   whopping=whopping, wsigma=wsigma, pixelmask=pixelmask, reject=reject

   ; Need 4 parameters
   if (N_params() LT 4) then begin
      print, 'Syntax - ans = extract_row( fimage, invvar, xcen, sigma, [ymodel='
      print, ' fscat=, proftype=, wfixed=, inputans=, iback=, bfixarr=,'
      print, ' xvar=, mask=, relative=, squashprofile=, diagonal=diagonal,'
      print, ' fullcovar=, wfixarr=, nPoly=, maxIter=, highrej=, lowrej=,'
      print, ' calcCovar=, niter=, whopping=, wsigma=])'
      return, -1
   endif

   nTrace = n_elements(xcen)
   nx = n_elements(fimage)

   if (n_elements(sigma) NE nTrace) then begin
    sigma1 = sigma[0]
    sigma = xcen*0.0 + sigma1
   endif 

   if (n_elements(nPoly) EQ 0) then nPoly = 5
   if (NOT keyword_set(maxIter)) then maxIter = 10
   if (NOT keyword_set(highrej)) then highrej = 15.0
   if (NOT keyword_set(lowrej)) then lowrej = 20.0
   if (NOT keyword_set(wfixed)) then wfixed = [1]
   if (NOT keyword_set(calcCovar)) then calcCovar = 0
   if (NOT keyword_set(proftype)) then proftype = 1

   if (n_elements(reject) EQ 2) then begin
     checkreject = sort([0.0,reject,1.0])
     if (total(abs(checkreject - [0,2,1,3])) NE 0) then reject = [0.8,0.2] 
   endif else reject = [0.8,0.2]

   ;------------------------------------------------------------------------
   ;   Pixelmask is an int array (I guess)
   ;    
   if (n_elements(pixelmask) NE nTrace OR size(pixelmask,/type) NE 3) then $
         pixelmask = lonarr(nTrace)

   relative = keyword_set(relative) 
   squashprofile =keyword_set(squashprofile) 

   if (NOT keyword_set(whopping)) then begin
      whoppingct = 0
      whopping = -1.0
   endif else if (whopping[0] LT 0.0) then whoppingct = 0 $
   else whoppingct = n_elements(whopping)

   if (NOT keyword_set(wsigma)) then wsigma = 25.0

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
     x1 = 0
     x2 = nTrace-1
     check = where(xcen[x1:x2-1] GE xcen[x1+1:x2] - 3,count)
     ma = nPoly + nTrace*nCoeff + whoppingct

   if (count GT 0) then $
      message, 'XCEN is not sorted or not separated by greater than 3 pixels.'


   whopping = float(whopping)
   whoppingct = LONG(whoppingct)
   nPoly = LONG(nPoly)

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
      if (keyword_set(bfixarr)) then $
              wfixarr(nTrace*nCoeff:nTrace*nCoeff + nPoly - 1) = bfixarr
   endif else if (ma NE n_elements(wfixarr)) then $
      message, 'Number of elements in FIMAGE and WFIXARR must be equal'

   ans = fltarr(ma)       ; parameter values
   p = fltarr(ma)         ; diagonal errors

   nonzerovar = where(invvar GT 0.0, ngood)
   reducedChi = 0.0
   niter = 0
   if (ngood EQ 0) then return, ans

   if keyword_set(iback) then begin
      if (nPoly NE n_elements(iback)) then $
         message, 'Number of elements in IBACK is not equal to nPoly'
      ans(nTrace*nCoeff:nTrace*nCoeff + nPoly-1) = iback 
;      wfixarr(nTrace*nCoeff:nTrace*nCoeff + nPoly - 1) = 0
   endif


   if (NOT arg_present(fullcovar)) then begin
       print, 'watch out'
       fullcovar = fltarr(ma,ma)  ; full covariance matrix
   endif

   finished = 0
   niter = 0
   totalreject = 0

   while(finished NE 1) do begin 

      workinvvar = float(invvar*mask)
      partial = lonarr(nTrace)
      fullreject = lonarr(nTrace)

     if keyword_set(inputans) then begin
       if (ma-nPoly-whoppingct NE n_elements(inputans)) then $
         message, 'Number of elements in INPUTANS is not equal to nTrace*nCoeff'
       ans(0:ma-nPoly-whoppingct-1) = inputans
     endif

      result = call_external(getenv('IDLSPEC2D_DIR')+'/lib/libspec2d.so', $
       'extract_row',$
       nx, float(xvar), float(fimage), workinvvar, float(ymodel), nTrace, $
       nPoly, float(xcen), float(sigma), proftype, float(reject), partial, $
       fullreject, calcCovar, squashprofile, whopping, whoppingct, float(wsigma), $
       nCoeff, ma, ans, long(wfixarr), p, fscat, fullcovar)

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
          mask[badhigh] = 0
          totalreject = totalreject + badhighct
          finished = 0
       endif 
       if(badlowct GT 0) then begin
          mask[badlow] = 0
          totalreject = totalreject + badlowct
          finished = 0
       endif

       diffs = (fimage - ymodel)*sqrt(invvar) 
       if (finished EQ 0) then begin
         good = where(diffs GE -lowrej*scaleError AND diffs LE highrej*scaleError, goodct)
         if (goodct GT 0) then mask[good] = 1
       endif

       niter = niter + 1
       if (niter EQ maxIter) then finished = 1
;       print,format='($,i)', niter
   endwhile

   ;----------------------------------------------------------------------------
   ;  Sort out pixelmask

   pixelmask = pixelmask OR (pixelmask_bits('PARTIALREJECT') * fix(partial))
   pixelmask = pixelmask OR (pixelmask_bits('FULLREJECT') * fix(fullreject))

; HORRIBLE HACK SINCE EXTRACT_ROW RETURNS SOME NaN's!!!???
jj = where(finite(ans) EQ 0)
if (jj[0] NE -1) then ans[jj] = 0
jj = where(finite(ymodel) EQ 0)
if (jj[0] NE -1) then ymodel[jj] = 0
jj = where(finite(fscat) EQ 0)
if (jj[0] NE -1) then fscat[jj] = 0

   return, ans

end
;------------------------------------------------------------------------------
