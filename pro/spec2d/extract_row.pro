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
;              mask=, relative=, diagonal=, fullcovar=, wfixarr=, npoly=,
;              maxiter=, lowrej=, highrej=, niter=, squashprofile=,
;              whopping=, wsigma=, pixelmask=, reject= ])
;
; INPUTS:
;   fimage     - Vector [nCol]
;   invvar     - Inverse variance [nCol]
;   xcen       - Initial guesses for X centers [nFiber]
;   sigma      - Sigma of gaussian profile; (scalar or [nFiber])
;
; OPTIONAL KEYWORDS:
;   proftype   - Select profile type:
;                  1: Gaussian
;                  2: (exponential)^3
;                  3: (exponential)^2.5
;                Default to 1.
;   wfixed     - Array to describe which parameters to fix in the profile;
;                0=fixed, 1=float; default to [1]. ???
;                The number of parameters to fit per fiber is determined
;                this way; e.g. nCoeff = n_elements(wfixed), so the default
;                is to fit only 1 parameter per fiber.
;   inputans   - 2D array of input answers [ncoeff,nFiber]
;   iback      - 1D array of input background coeff 
;                (needed if fixed parameters are non-zero)
;   bfixarr    - array of 1's and zero's which set which background 
;                parameters are fixed. ??? which is 1 or 0?
;   wfixarr    - 1D integer array of 1's and zero's which specify fixed 
;                profile parameters ??? what the fk???
;   xvar       - X values of fimage and invvar; default is findgen(NX).
;   mask       - Image mask: 1=good, 0=bad [NX]
;   relative   - Set to use reduced chi-square to scale rejection threshold
;   squashprofile - ???
;   npoly      - Order of chebyshev scattered light background; default to 5
;   maxiter    - Maximum number of profile fitting iterations; default to 10
;   lowrej     - Negative sigma deviation to be rejected; default to 5
;   highrej    - Positive sigma deviation to be rejected; default to 5
;   whopping   - X locations to center additional "whopping" terms to describe
;                the exponentail tails of flux near bright fibers; default
;                to -1, which means not to use any such terms.
;   wsigma     - Sigma width for exponential whopping profiles; default to 25
;   pixelmask  - Bits set due to extraction rejection [nFiber]
;   reject     - Two-element array setting partial and full rejection
;                thresholds for profiles; default [0.8, 0.2].
;                What does this mean???
;
; OUTPUTS:
;   ans        - Extracted flux in each parameter [ncoeff, nFiber]
;
; OPTIONAL OUTPUTS:
;   ymodel     - Evaluation of best fit [nCol]
;   fscat      - Scattered light contribution in each fiber [nFiber]
;   diagonal   - 1D diagonal of covariance matrix.  Currently, this is
;                the diagonal from the cholesky decompostion, which is
;                1/error[j].
;   fullcovar  - 2D covariance matrix.  This is a symmetric matrix, and we
;                only fill the lower triangle and set the rest to 0.
;                Computing this increases CPU time by a factor of 2 or 3.
;   niter      - Number of rejection iterations performed
;   pixelmask  - (Modified.)
;
; COMMENTS:
;
; BUGS:
;    Still need to do:
;       limits on chebyshev polynomial are assumed to be 0.0 <--> nx
;       these may need to be optional if only partial rows are being fit
;
;       Error codes need to be returned, currently no such codes are returned
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   Dynamic link to extract_row.c
;
; REVISION HISTORY:
;    8-Aug-1999  Written by Scott Burles, Chicago 
;-
;------------------------------------------------------------------------------
function extract_row, fimage, invvar, xcen, sigma, ymodel=ymodel, $
   fscat=fscat, proftype=proftype, wfixed=wfixed, inputans=inputans, $
   iback=iback, bfixarr=bfixarr, xvar=xvar, mask=mask, relative=relative, $
   squashprofile=squashprofile, diagonal=p, fullcovar=fullcovar, $
   wfixarr=wfixarr, npoly=npoly, maxiter=maxiter, $
   lowrej=lowrej, highrej=highrej, niter=niter, reducedChi=reducedChi, $
   whopping=whopping, wsigma=wsigma, pixelmask=pixelmask, reject=reject

   ; Need 4 parameters
   if (N_params() LT 4) then $
    message, 'Wrong number of parameters'

   ntrace = n_elements(xcen)
   nx = n_elements(fimage)

   if (n_elements(sigma) NE ntrace) then begin
      sigma1 = sigma[0]
      sigma = xcen*0.0 + sigma1
   endif 

   if (n_elements(npoly) EQ 0) then npoly = 5
   if (NOT keyword_set(maxiter)) then maxiter = 10
   if (NOT keyword_set(highrej)) then highrej = 15.0
   if (NOT keyword_set(lowrej)) then lowrej = 20.0
   if (NOT keyword_set(wfixed)) then wfixed = [1]
   if (NOT keyword_set(proftype)) then proftype = 1 ; Gaussian
   relative = keyword_set(relative) 
   squashprofile = keyword_set(squashprofile) 
   if (NOT keyword_set(wsigma)) then wsigma = 25.0

   if (n_elements(reject) EQ 2) then begin
     checkreject = sort([0.0,reject,1.0])
     if (total(abs(checkreject - [0,2,1,3])) NE 0) then reject = [0.8,0.2] 
   endif else reject = [0.8,0.2]

   if (n_elements(pixelmask) NE ntrace OR size(pixelmask,/type) NE 3) then $
    pixelmask = lonarr(ntrace)

   if (NOT keyword_set(whopping)) then whopping = -1
   if (whopping[0] EQ -1) then whoppingct = 0 $
    else whoppingct = n_elements(whopping)

   if (NOT keyword_set(xvar)) then xvar = findgen(nx) $
    else if (nx NE n_elements(xvar)) then $
     message, 'Number of elements in FIMAGE and XVAR must be equal'

   if (NOT keyword_set(mask)) then mask = make_array(nx, /byte, value=1) $
    else if (nx NE n_elements(mask)) then $
     message, 'Number of elements in FIMAGE and MASK must be equal'

   ncoeff = n_elements(wfixed)

   if (nx NE n_elements(invvar)) then $
    message, 'Number of elements in FIMAGE and INVVAR must be equal'

   ;----------
   ; Check that XCEN is sorted in increasing order
   ; with separations of at least 3 pixels.

   junk = where(xcen[0:ntrace-2] GE xcen[1:ntrace-1] - 3, ct)
   if (ct GT 0) then $
    message, 'XCEN is not sorted or not separated by greater than 3 pixels.'

   ;----------
   ; Allocate memory for the C subroutine.

   ymodel = fltarr(nx)
   fscat = fltarr(ntrace)
   ma = npoly + ntrace*ncoeff + whoppingct

   if (NOT keyword_set(wfixarr)) then begin
      wfixarr = lonarr(ma) + 1        ; Fixed parameter array
      i = 0
      wfixarr[lindgen(ntrace)*ncoeff+i] = wfixed[i] 
      for i=1, ncoeff-1 do $
       wfixarr[lindgen(ntrace)*ncoeff+i] = wfixed[i] * (1 - squashprofile)
      if (keyword_set(bfixarr)) then $
       wfixarr[ntrace*ncoeff:ntrace*ncoeff + npoly - 1] = bfixarr
   endif else begin
      if (ma NE n_elements(wfixarr)) then $
       message, 'Number of elements in FIMAGE and WFIXARR must be equal'
   endelse

   ans = fltarr(ma)       ; parameter values
   p = fltarr(ma)         ; diagonal errors

   ; Set the following variables before any possible RETURN

   nonzerovar = where(invvar GT 0.0, ngood)
   reducedChi = 0.0
   niter = 0

   if (ngood EQ 0) then return, ans

   if (keyword_set(iback)) then begin
      if (npoly NE n_elements(iback)) then $
       message, 'Number of elements in IBACK is not equal to NPOLY'
      ans[ntrace*ncoeff:ntrace*ncoeff + npoly-1] = iback 
;      wfixarr[ntrace*ncoeff:ntrace*ncoeff + npoly - 1] = 0
   endif

   if (arg_present(fullcovar)) then qcovar = 1L $
    else qcovar = 0L
   fullcovar = fltarr(ma,ma)

   finished = 0
   totalreject = 0

   while(finished NE 1) do begin 

      workinvvar = FLOAT(invvar * mask)
      partial = lonarr(ntrace)
      fullreject = lonarr(ntrace)

      if (keyword_set(inputans)) then begin
         if (ma-npoly-whoppingct NE n_elements(inputans)) then $
          message, 'Number of elements in INPUTANS is not equal to NTRACE*NCOEFF'
         ans[0:ma-npoly-whoppingct-1] = inputans
      endif

      result = call_external(getenv('IDLSPEC2D_DIR')+'/lib/libspec2d.so', $
       'extract_row',$
       nx, FLOAT(xvar), FLOAT(fimage), workinvvar, FLOAT(ymodel), ntrace, $
       LONG(npoly), FLOAT(xcen), FLOAT(sigma), LONG(proftype), $
       FLOAT(reject), partial, fullreject, qcovar, LONG(squashprofile), $
       FLOAT(whopping), whoppingct, FLOAT(wsigma), $
       ncoeff, ma, ans, LONG(wfixarr), p, fscat, fullcovar)

      diffs = (fimage - ymodel) * sqrt(workinvvar) 
      chisq = total(diffs*diffs)
      countthese = total(wfixarr)
      reducedChi = chisq / (total(mask) - countthese)
      errscale = 1.0
      if (relative) then errscale = sqrt(chisq/total(mask))

      badhigh = where(diffs GT highrej*errscale, badhighct)
      badlow = where(diffs LT -lowrej*errscale, badlowct)

      finished = 1
      if (badhighct GT 0) then begin
         mask[badhigh] = 0
         totalreject = totalreject + badhighct
         finished = 0
      endif 
      if (badlowct GT 0) then begin
         mask[badlow] = 0
         totalreject = totalreject + badlowct
         finished = 0
      endif

      diffs = (fimage - ymodel) * sqrt(invvar) 
      if (finished EQ 0) then begin
         igood = where(diffs GE -lowrej*errscale $
          AND diffs LE highrej*errscale, goodct)
         if (goodct GT 0) then mask[igood] = 1
      endif

      niter = niter + 1
      if (niter EQ maxiter) then finished = 1
   endwhile

   ;----------
   ; Add bits to PIXELMASK

   pixelmask = pixelmask OR (pixelmask_bits('PARTIALREJECT') * fix(partial))
   pixelmask = pixelmask OR (pixelmask_bits('FULLREJECT') * fix(fullreject))

; HORRIBLE HACK SINCE EXTRACT_ROW RETURNS SOME NaN's!!!???
jj = where(finite(ans) EQ 0)
if (jj[0] NE -1) then ans[jj] = 0
;if (jj[0] NE -1) then stop
jj = where(finite(ymodel) EQ 0)
if (jj[0] NE -1) then ymodel[jj] = 0
;if (jj[0] NE -1) then stop
jj = where(finite(fscat) EQ 0)
if (jj[0] NE -1) then fscat[jj] = 0
;if (jj[0] NE -1) then stop

   return, ans

end
;------------------------------------------------------------------------------
