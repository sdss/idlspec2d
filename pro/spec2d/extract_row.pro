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
;              whopping=, wsigma=, pixelmask=, reject=, reducedChi= ])
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
;   inputans   - Input fit, excluding background and whopping terms
;                [ncoeff*nFiber]
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
;   reject     - Two-element array setting partial and full rejection
;                thresholds for profiles; default [0.8, 0.2].
;                What does this mean???
;                When this was hardwired, it was [0.8,0.4].
;
; MODIFIED INPUTS (OPTIONAL):
;   wfixed     - Array to describe which parameters to fix in the profile;
;                0=fixed, 1=float; default to [1].
;                The number of parameters to fit per fiber is determined
;                this way; e.g. nCoeff = n_elements(wfixed), so the default
;                is to fit only 1 parameter per fiber.  For the (default)
;                Gaussian profile, this is the height of the Gaussian.
;                Note that WFIXED is used to build the array WFIXARR.
;   iback      - 1D array of input background coeff 
;                (needed if fixed parameters are non-zero)
;   bfixarr    - 1D integer array to specify which terms of the background
;                coefficients to fix; 0=fixed, 1=float.
;   wfixarr    - 1D integer array to specify which parameters in the full fit
;                to fix; 0=fixed, 1=float.
;                The array is sorted as follows:
;                   [ncoeff] values for fiber #0
;                   [ncoeff] values for fiber #1
;                   ...
;                   [ncoeff] values for fiber #(nFiber-1)
;                   [npoly] values for the background polynomial terms
;                   [whoppingct] values for the whopping terms
;   xvar       - X values of fimage and invvar; default is findgen(NX).
;   mask       - Image mask: 1=good, 0=bad [NX]
;   pixelmask  - Bits set for each fiber due to extraction rejection [nFiber]
;
; OUTPUTS:
;   ans        - Output fit [ncoeff*nFiber+npoly+whoppingct]
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
;   reducedChi - Reduced chi ???
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
function extract_row1, fimage, invvar, xcen, sigma, ymodel=ymodel, $
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

   if (n_elements(pixelmask) NE ntrace $
    OR size(pixelmask,/tname) NE 'LONG') then $
    pixelmask = lonarr(ntrace)

   if (NOT keyword_set(whopping)) then whopping = -1
   if (whopping[0] EQ -1) then whoppingct = 0 $
    else whoppingct = n_elements(whopping)

   if (NOT keyword_set(xvar)) then xvar = findgen(nx) $
    else if (nx NE n_elements(xvar)) then $
     message, 'Number of elements in FIMAGE and XVAR must be equal'

   if (NOT keyword_set(mask)) then mask = bytarr(nx) + 1b $
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
   ma = ntrace*ncoeff + npoly + whoppingct

   ;----------
   ; Test which points are good

   qgood = invvar GT 0.0 AND mask NE 0
   igood = where(qgood, ngood)

   ;----------
   ; Set the following variables before any possible RETURN statement.

   reducedChi = 0.0
   niter = 0
   ans = fltarr(ma)       ; parameter values
   p = fltarr(ma)         ; diagonal errors

   if (ngood EQ 0) then return, ans

   ;----------
   ; Build the fixed parameter array if it was not passed.

   if (NOT keyword_set(wfixarr)) then begin
      wfixarr = lonarr(ma) + 1

      ; Set values for the (gaussian) profile terms
      i = 0
      wfixarr[lindgen(ntrace)*ncoeff+i] = wfixed[i] 
      for i=1, ncoeff-1 do $
       wfixarr[lindgen(ntrace)*ncoeff+i] = wfixed[i] * (1 - squashprofile)

      ; Set values for the background polynomial terms
      if (keyword_set(bfixarr)) then $
       wfixarr[ntrace*ncoeff:ntrace*ncoeff + npoly - 1] = bfixarr

      ; Disable fitting to any profiles out of the data range.
      ; If there are more than one XCEN profile centers to the left of
      ; the first good data point, then reject all but the first
      ; out-of-range profile.  Do the same to the right side.

;      ileft = where(xcen LT xvar[igood[0]], nleft)
;      if (nleft GT 1) then $
;       wfixarr[0:(nleft-1)*ncoeff-1] = 0
      ; Instead, look for any XCEN more than 2 pix off
      ileft = where(xcen LT xvar[igood[0]]-2.0, nleft)
      if (nleft GT 0) then $
       wfixarr[0:nleft*ncoeff-1] = 0

;      iright = where(xcen GT xvar[igood[ngood-1]], nright)
;      if (nright GT 1) then $
;       wfixarr[(ntrace-nright+1)*ncoeff : ntrace*ncoeff-1] = 0
      ; Instead, look for any XCEN more than 2 pix off
      iright = where(xcen GT xvar[igood[ngood-1]]+2.0, nright)
      if (nright GT 0) then $
       wfixarr[(ntrace-nright)*ncoeff : ntrace*ncoeff-1] = 0

      ; Don't fit to any centers where there are no good data points
      ; within 2.0 pix
      for i=0, ntrace-1 do begin
         ii = where(abs(xvar[igood] - xcen[i]) LT 2.0, nn)
         if (nn EQ 0) then $
          wfixarr[i*ncoeff : i*ncoeff+ncoeff-1] = 0
      endfor
   endif else begin
      if (ma NE n_elements(wfixarr)) then $
       message, 'Number of elements in FIMAGE and WFIXARR must be equal'
      wfixarr = LONG(wfixarr)
   endelse

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
         if (ntrace*ncoeff NE n_elements(inputans)) then $
          message, 'Number of elements in INPUTANS is not equal to NTRACE*NCOEFF'
         ans[0:ntrace*ncoeff-1] = inputans
      endif

      result = call_external(getenv('IDLSPEC2D_DIR')+'/lib/libspec2d.so', $
       'extract_row',$
       nx, FLOAT(xvar), FLOAT(fimage), workinvvar, ymodel, ntrace, $
       LONG(npoly), FLOAT(xcen), FLOAT(sigma), LONG(proftype), $
       FLOAT(reject), partial, fullreject, qcovar, LONG(squashprofile), $
       FLOAT(whopping), whoppingct, FLOAT(wsigma), $
       ncoeff, ma, ans, wfixarr, p, fscat, fullcovar)

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
         ii = where(diffs GE -lowrej*errscale $
          AND diffs LE highrej*errscale)
         if (ii[0] NE -1) then mask[ii] = 1 ; These points still good
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
;if (jj[0] NE -1) then ans[jj] = 0
if (jj[0] NE -1) then stop
jj = where(finite(ymodel) EQ 0)
;if (jj[0] NE -1) then ymodel[jj] = 0
if (jj[0] NE -1) then stop
jj = where(finite(fscat) EQ 0)
;if (jj[0] NE -1) then fscat[jj] = 0
if (jj[0] NE -1) then stop
;if (total(ymodel) EQ 0) then stop

   return, ans
end

;------------------------------------------------------------------------------
; DJS hack to not pass fibers outside of the range of data to EXTRACT_ROW.
; I have not implemented returning FULLCOVAR out of sheer laziness.

function extract_row, fimage, invvar, xcen, sigma, ymodel=ymodel, $
 fscat=fscat, wfixed=wfixed, inputans=inputans, xvar=xvar, $
 mask=mask, diagonal=p, wfixarr=wfixarr, npoly=npoly, $
 niter=niter, whopping=whopping, pixelmask=pixelmask, $
 reducedChi=reducedChi, _EXTRA=extra

   ntrace = n_elements(xcen)
   nx = n_elements(fimage)

   if (n_elements(npoly) EQ 0) then npoly = 5
   if (NOT keyword_set(wfixed)) then wfixed = [1]

   if (n_elements(pixelmask) NE ntrace $
    OR size(pixelmask,/tname) NE 'LONG') then $
    pixelmask = lonarr(ntrace)

   if (NOT keyword_set(whopping)) then whopping = -1
   if (whopping[0] EQ -1) then whoppingct = 0 $
    else whoppingct = n_elements(whopping)

   if (NOT keyword_set(xvar)) then xvar = findgen(nx) $
    else if (nx NE n_elements(xvar)) then $
     message, 'Number of elements in FIMAGE and XVAR must be equal'

   if (NOT keyword_set(mask)) then mask = bytarr(nx) + 1b $
    else if (nx NE n_elements(mask)) then $
     message, 'Number of elements in FIMAGE and MASK must be equal'

   ncoeff = n_elements(wfixed)

   if (nx NE n_elements(invvar)) then $
    message, 'Number of elements in FIMAGE and INVVAR must be equal'

   ma = ntrace*ncoeff + npoly + whoppingct

   ;----------
   ; Test which points are good

   qgood = invvar GT 0.0 AND mask NE 0
   igood = where(qgood, ngood)

   ;----------
   ; Set the following variables before any possible RETURN statement.

   reducedChi = 0.0
   niter = 0
   ans = fltarr(ma)       ; parameter values
   p = fltarr(ma)         ; diagonal errors

   if (ngood EQ 0) then return, ans

   ;----------
   ; Disable fitting to any profiles out of the data range.

   fibuse = bytarr(ntrace) + 1

   ; Instead, look for any XCEN more than 0 pix off
   ileft = where(xcen LT xvar[igood[0]]-0.0, nleft)
   if (nleft GT 0) then fibuse[ileft] = 0

   ; Instead, look for any XCEN more than 0 pix off
   iright = where(xcen GT xvar[igood[ngood-1]]+0.0, nright)
   if (nright GT 0) then fibuse[iright] = 0

   ; Don't fit to any centers where there are no good data points
   ; within 2.0 pix
   for i=0, ntrace-1 do begin
      ii = where(abs(xvar[igood] - xcen[i]) LT 2.0, nn)
      if (nn EQ 0) then fibuse[i] = 0
   endfor

   iuse = where(fibuse, nuse)
   tmp_ma = nuse*ncoeff + npoly + whoppingct

   if (keyword_set(inputans)) then begin
      tmp_inputans = fltarr(nuse*ncoeff)
      for i=0, ncoeff-1 do $
       tmp_inputans[lindgen(nuse)*ncoeff+i] = inputans[iuse*ncoeff+i]
   endif else begin
      tmp_inputans = 0
   endelse

   if (keyword_set(wfixarr)) then begin
      tmp_wfixarr = lonarr(tmp_ma) + 1
      for i=0, ncoeff-1 do $
       tmp_wfixarr[lindgen(nuse)*ncoeff+i] = wfixarr[iuse*ncoeff+i]
   endif else begin
      tmp_wfixarr = 0
   endelse

ymodel = 0 ; ???
   tmp_ans = extract_row1( fimage, invvar, xcen[iuse], sigma, ymodel=ymodel, $
    fscat=tmp_fscat, wfixed=wfixed, inputans=tmp_inputans, xvar=xvar, $
    mask=mask, diagonal=tmp_p, wfixarr=tmp_wfixarr, npoly=npoly, $
    whopping=whopping, reducedChi=reducedChi, _EXTRA=extra)

   ; Set WFIXARR for unused fibers equal to zero
   wfixarr = lonarr(ma) + 0
   for i=0, ncoeff-1 do $
    wfixarr[iuse*ncoeff+i] = tmp_wfixarr[lindgen(nuse)*ncoeff+i]

   ; Set P for unused fibers equal to zero
   p = fltarr(ma) + 0
   for i=0, ncoeff-1 do $
    p[iuse*ncoeff+i] = tmp_p[lindgen(nuse)*ncoeff+i]

   ; In FSCAT, linearly interpolate over unused fibers
   fscat = fltarr(ntrace)
   fscat[iuse] = tmp_fscat
   fscat = djs_maskinterp(fscat, fibuse EQ 0)

   ; Set ANS for unused fibers equal to zero
   ans = fltarr(ma)
   for i=0, ncoeff-1 do $
    ans[iuse*ncoeff+i] = tmp_ans[lindgen(nuse)*ncoeff+i]

   nextra = npoly + whoppingct
   if (nextra GT 0) then begin
      wfixarr[ma-nextra:ma-1] = tmp_wfixarr[tmp_ma-nextra:tmp_ma-1]
      ans[ma-nextra:ma-1] = tmp_ans[tmp_ma-nextra:tmp_ma-1]
      p[ma-nextra:ma-1] = tmp_p[tmp_ma-nextra:tmp_ma-1]
   endif

   return, ans
end
;------------------------------------------------------------------------------
