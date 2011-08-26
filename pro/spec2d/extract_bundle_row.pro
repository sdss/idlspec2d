;+
; NAME:
;   extract_bundle_row
;
; PURPOSE:
;   Fit the fiber profiles and background in a single row with least
;   squares, chunking it up by bundle, using Gaussians plus polynomial
;   background.
;
;   Based upon extract_row.pro, which was a call to extract_row.c
;
; NOTES ON CONVERSION TO EXTRACT_BUNDLE_ROW (A. Bolton, Utah, 2011Aug):
;   The extraction guts of this program are entirely new,
;   written in pure IDL rather than in IDL-wrapped C.
;   However, I have attempted (1) to preserve most of the same
;   outward functionality in terms of inputs and outputs, at least
;   those that are relevant to the BOSS case, and (2) to retain
;   the same rejection-iteration behavior, so that the same sorts
;   of cases trigger rejection as before.  I did not parse all
;   the details of the rejection logic.  One thing I know is that
;   the rejection in current form would be more efficient if it were
;   done on each bundle separately, rather than on the entire row.
;   But that would require a fine-grained refactoring that would risk
;   breaking the rejection logic, so I didn't attempt it.
;
;   Some deprecated (and even retained) keyword input/output variables
;   may have unexpected behavior is used in a context other than the
;   current BOSS-implementation calls in EXTRACT_IMAGE.
;
;   I have retained a fair bit of the original code in commented-out
;   form, for forensic purposes.  Anything removed during the
;   conversion to bundle form is commented out with a ";#".
;
; CALLING SEQUENCE:
;; First four positional arguments inherited from extract_row:
;   ans = extract_row( fimage, invvar, xcen, sigma,
;; Next keyword arguments new to extract_bundle_row:
;; (skew and kurt not currently enabled)
;              [nperbun=, buffsize=, skew=, kurt=,
;; Next keyword retained but reinterpreted:
;              npoly=, 
;; Next keywords retained in more or less the same role (I think):
;              maxiter=, lowrej=, highrej=, niter=, reducedChi=, xvar=,
;              mask=, relative=, diagonal=, pixelmask=, reject=,, ymodel=
;; Next keyword arguments deprecated:
;              fscat=, proftype=, wfixed=, inputans=, iback=, bfixarr=,
;              fullcovar=, wfixarr=, squashprofile=, whopping=,
;              wsigma=, nBand=  ])
;
; INPUTS:
;   fimage     - Vector [nCol]
;   invvar     - Inverse variance [nCol]
;   xcen       - Initial guesses for X centers [nFiber]
;   sigma      - Sigma of gaussian profile; (scalar or [nFiber])
;
; OPTIONAL KEYWORDS (new):
;   nperbun    - Number of fibers per bundle (default is 20).
;   buffsize   - Pixel buffer size on the very outside of the image
;                (default is 8)
;   npoly      - order of polynomial background in bundle, default=2 (linear).
;
; OPTIONAL KEYWORDS (retained):
;   maxiter    - Maximum number of profile fitting iterations; default to 10
;   lowrej     - Negative sigma deviation to be rejected; default to 5
;   highrej    - Positive sigma deviation to be rejected; default to 5
;   reject     - Three-element array setting partial and full rejection
;                thresholds for profiles; default [0.2, 0.6, 0.6].
;                When there is less than REJECT[2] of the area is left,
;                  then drop fitting of all higher-order terms.
;                When there is less than REJECT[1] of the area is left,
;                  then the pixel is rejected (inverse variance is set to 0).
;                When there is less than REJECT[0] of the area is left,
;                  then assume that there's no fiber there, and don't fit
;                  for that fiber at all.
;   relative   - Set to use reduced chi-square to scale rejection threshold
;
; DEPRECATED KEYWORDS:
;   inputans   - Input fit, excluding background and whopping terms
;                [ncoeff*nFiber]
;                The array is sorted as follows:
;                  [ncoeff] values for fiber #0
;                   ...
;                  [ncoeff] values for fiber #(nFiber-1)
;   squashprofile - ???
;   nband      - Band-width of covariance matrix of fiber profiles: default 1
;   whopping   - X locations to center additional "whopping" terms to describe
;                the exponentail tails of flux near bright fibers; default
;                to -1, which means not to use any such terms.
;   wsigma     - Sigma width for exponential whopping profiles; default to 25
;
; MODIFIED INPUTS (OPTIONAL, retained):
;   xvar       - X values of fimage and invvar; default is findgen(NX).
;   mask       - Image mask: 1=good, 0=bad [NX]
;   pixelmask  - Bits set for each fiber due to extraction rejection
;                [nFiber]
;
; MODIFIED INPUTS (OPTIONAL, deprecated):
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
;                  [ncoeff] values for fiber #0
;                   ...
;                  [ncoeff] values for fiber #(nFiber-1)
;                  [npoly] values for the background polynomial terms
;                  [whoppingct] values for the whopping terms
;
; OUTPUTS (reinterpreted):
;   ans        - Output fit [nFiber]
;                (extracted Gaussian profile amplitudes)
;
; OPTIONAL OUTPUTS (retained):
;   ymodel     - Evaluation of best fit [nCol]
;   diagonal   - 1D diagonal of covariance matrix
;   niter      - Number of rejection iterations performed
;   reducedChi - Reduced chi ???
;
; OPTIONAL OUTPUTS (deprecated):
;   fscat      - Scattered light contribution in each fiber [nFiber]
;   fullcovar  - 2D covariance matrix.  This is a symmetric matrix, and we
;                only fill the lower triangle.  Computing this increases CPU
;                time by a factor of 2 or 3.
;
; REVISION HISTORY:
;    8-Aug-1999  extract_row Written by Scott Burles, Chicago 
;    Sep 2010    converted to extract_bundle_row by Adam S. Bolton, Utah
;    Aug 2011    attempted documentation cleanup, Adam S. Bolton, Utah
;-
;------------------------------------------------------------------------------
function extract_bundle_row, fimage, invvar, xcen, sigma, ymodel=ymodel, $
 fscat=fscat, proftype=proftype, wfixed=wfixed, inputans=inputans, $
 iback=iback, bfixarr=bfixarr, xvar=xvar, mask=mask, relative=relative, $
 squashprofile=squashprofile, diagonal=p, fullcovar=fullcovar, $
 wfixarr=wfixarr, npoly=npoly, maxiter=maxiter, $
 lowrej=lowrej, highrej=highrej, niter=niter, reducedChi=reducedChi, $
 whopping=whopping, wsigma=wsigma, pixelmask=pixelmask, reject=reject, $
 oldreject=oldreject, nband = nband, contribution=contribution, $
 nperbun=nperbun, buffsize=buffsize, skew=skew, kurt=kurt

   ; Need 4 parameters
   if (N_params() LT 4) then $
    message, 'Wrong number of parameters'

   ntrace = n_elements(xcen)
   nx = n_elements(fimage)

   if (n_elements(sigma) NE ntrace) then begin
      sigma1 = sigma[0]
      sigma = xcen*0.0 + sigma1
   endif 

; extract_row keywords:
   if (n_elements(npoly) EQ 0) then npoly = 2L ; order of background is now per bundle!
   if (NOT keyword_set(nband)) then nband = 1L
   if (NOT keyword_set(maxiter)) then maxiter = 10
   if (NOT keyword_set(highrej)) then highrej = 15.0
   if (NOT keyword_set(lowrej)) then lowrej = 20.0
   if (NOT keyword_set(wfixed)) then wfixed = [1]
   if (NOT keyword_set(proftype)) then proftype = 1 ; Gaussian
   relative = keyword_set(relative) 
   squashprofile = keyword_set(squashprofile) 
   if (NOT keyword_set(wsigma)) then wsigma = 25.0
; extract_bundle_row keywords:
   if (NOT keyword_set(nperbun)) then nperbun = 20L
   if (NOT keyword_set(buffsize)) then buffsize = 8L

   ; Here we want a three element array where both are between 0 and 1, and the
   ; first is larger than the second.  
   ; The first threshold sets the minimum area required to perform 
   ;    single parameter profile fitting
   ; The second threshold is the minimum area required not to reject 
   ;    the pixel in the final extracted spectrum.
   ; The third parameter is the area required
   ;    in the profile fit containing good pixels to do a full fit.

   if (n_elements(reject) NE 3) then reject = [0.2, 0.6, 0.6]

   if (n_elements(pixelmask) NE ntrace $
    OR size(pixelmask,/tname) NE 'LONG') then $
    pixelmask = lonarr(ntrace)

   if (NOT keyword_set(whopping)) then whopping = -1
   if (whopping[0] EQ -1) then whoppingct = 0L $
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
   ; Allocate memory for the C subroutine.

   ymodel = fltarr(nx)
;#   fscat = fltarr(ntrace)
;#   ma = ntrace*ncoeff + npoly + whoppingct

   ;----------
   ; Test which points are good

   qgood = invvar GT 0.0 AND mask NE 0
   igood = where(qgood, ngood)

   ;----------
   ; Set the following variables before any possible RETURN statement.

   reducedChi = 0.0
   niter = 0
;#   ans = fltarr(ma)       ; parameter values
;#   p = fltarr(ma)         ; diagonal errors
ans = fltarr(ntrace)
p = fltarr(ntrace)

   if (ngood EQ 0) then return, ans

   ;----------
   ; Check that XCEN is sorted in increasing order
   ; with separations of at least 3 pixels.

   junk = where(xcen[0:ntrace-2] GE xcen[1:ntrace-1] - 3, ct)
   if (ct GT 0) then $
;    message, 'XCEN is not sorted or not separated by greater than 3 pixels.'
; Should definitely reject here!!!???
    splog, 'XCEN is not sorted or not separated by greater than 3 pixels.'

   ;----------
   ; Build the fixed parameter array if it was not passed.

;#   if (NOT keyword_set(wfixarr)) then begin
;#      wfixarr = lonarr(ma) + 1
      wfixarr = lonarr(ntrace) + 1

      ; Set values for the (gaussian) profile terms
;#      i = 0
;#      wfixarr[lindgen(ntrace)*ncoeff+i] = wfixed[i] 
;#      for i=1, ncoeff-1 do $
;#       wfixarr[lindgen(ntrace)*ncoeff+i] = wfixed[i] * (1 - squashprofile)

      ; Set values for the background polynomial terms
;#      if (keyword_set(bfixarr)) then $
;#       wfixarr[ntrace*ncoeff:ntrace*ncoeff + npoly - 1] = bfixarr

      ; Disable fitting to any profiles out of the data range.
      ; If there are more than one XCEN profile centers to the left of
      ; the first good data point, then reject all but the first
      ; out-of-range profile.  Do the same to the right side.

;      ileft = where(xcen LT xvar[igood[0]], nleft)
;      if (nleft GT 1) then $
;       wfixarr[0:(nleft-1)*ncoeff-1] = 0
      ; Instead, look for any XCEN more than 2 pix off
;      ileft = where(xcen LT xvar[igood[0]]-2.0, nleft)
;      if (nleft GT 0) then $
;       wfixarr[0:nleft*ncoeff-1] = 0

;      iright = where(xcen GT xvar[igood[ngood-1]], nright)
;      if (nright GT 1) then $
;       wfixarr[(ntrace-nright+1)*ncoeff : ntrace*ncoeff-1] = 0
      ; Instead, look for any XCEN more than 2 pix off
;      iright = where(xcen GT xvar[igood[ngood-1]]+2.0, nright)
;      if (nright GT 0) then $
;       wfixarr[(ntrace-nright)*ncoeff : ntrace*ncoeff-1] = 0

      ; Don't fit to any centers where there are no good data points
      ; within 2.0 pix
;      for i=0, ntrace-1 do begin
;         ii = where(abs(xvar[igood] - xcen[i]) LT 2.0, nn)
;         if (nn EQ 0) then $
;          wfixarr[i*ncoeff : i*ncoeff+ncoeff-1] = 0
;      endfor
;#   endif else begin
;#      if (ma NE n_elements(wfixarr)) then $
;#       message, 'Number of elements in FIMAGE and WFIXARR must be equal'
;#      wfixarr = LONG(wfixarr)
;#   endelse

; DISABLING "IBACK":
;#   if (keyword_set(iback)) then begin
;#      if (npoly NE n_elements(iback)) then $
;#       message, 'Number of elements in IBACK is not equal to NPOLY'
;#      ans[ntrace*ncoeff:ntrace*ncoeff + npoly-1] = iback 
;      wfixarr[ntrace*ncoeff:ntrace*ncoeff + npoly - 1] = 0
;#   endif

; DISABLING "FULLCOVAR":
;#   if (arg_present(fullcovar)) then qcovar = 1L $
;#    else qcovar = 0L
;#   fullcovar = fltarr(ma,ma)

; PUTTING IN LOOP OVER BUNDLES:
; Find number of bundles:
nbun = ntrace / nperbun
; Find breaks between the bundles
; and set limits of pixels to be associated with each bundle:
t_lo = nperbun * lindgen(nbun)
t_hi = t_lo + nperbun - 1
xc_lo = xcen[t_lo]
xc_hi = xcen[t_hi]
if (nbun gt 1) then begin
   midpoints = 0.5 * (xc_lo[1:nbun-1] + xc_hi[0:nbun-2])
   jmax = floor(midpoints)
   jmin = jmax + 1
   jmin = [(floor(xc_lo[0]) - buffsize) > 0, jmin]
   jmax = [jmax, (ceil(xc_hi[nbun-1]) + buffsize) < (nx - 1)]
endif else begin
   jmin = (floor(xc_lo) - buffsize) > 0
   jmax = (ceil(xc_hi) + buffsize) < (nx - 1)
endelse

; Mask outside the bundle-fitting zone:
mask = mask * (xvar ge min(jmin)) * (xvar le max(jmax))

; The loop over bundles:
   totalreject = 0
   partial = lonarr(ntrace)
   fullreject = lonarr(ntrace)
   finished = 0
   fiberbase = findgen(nperbun)

   while(finished NE 1) do begin 

      ; Apply current mask to invvar:
      workinvvar = (FLOAT(invvar * mask))

      for ibun = 0L, nbun-1 do begin

; Make the extracting basis:
         workxcen = xcen[t_lo[ibun]:t_hi[ibun]]
         worksigma = sigma[t_lo[ibun]:t_hi[ibun]]
         workpar = (transpose([[workxcen], [worksigma]]))[*]
         workx = xvar[jmin[ibun]:jmax[ibun]]
         workxpoly = 2. * (workx - min(workx)) / (max(workx) - min(workx)) - 1.
         workbasis = gausspix(workx, workpar)
         if (npoly gt 0) then workbasis = [[workbasis], [flegendre(workxpoly, npoly)]]
         workimage = fimage[jmin[ibun]:jmax[ibun]]

; Associate pixels with fibers:
         pixelfiber = round(interpol(fiberbase, workxcen, workx))

; Determine the good area fraction for each fiber,
; limiting to pixels "associated with" that fiber:
         bworkinvvar = workinvvar[jmin[ibun]:jmax[ibun]]
         goodarea = 0. * fltarr(nperbun)
         for ijk = 0L, nperbun-1 do begin
            totalarea = total((pixelfiber eq ijk) * (workbasis[*,ijk]))
            totalgood = total((pixelfiber eq ijk) * (bworkinvvar gt 0.) * (workbasis[*,ijk]))
            if (totalarea gt 0.) then goodarea[ijk] = totalgood / totalarea
         endfor

; Populate the rejection masks:
         partial[t_lo[ibun]:t_hi[ibun]] = goodarea lt reject[2]
         fullreject[t_lo[ibun]:t_hi[ibun]] = goodarea lt reject[1]

; Remove fitting for any fiber with area less than reject[0]:
         use_component = where(goodarea ge reject[0], n_use)
; Don't bother fitting at all if there are no unmasked fibers:
if (n_use gt 0) then begin

         if (npoly gt 0) then use_component = [use_component, nperbun + lindgen(npoly)]
         nfit_this = n_elements(use_component)
         workbasis = workbasis[*,use_component]

; The extraction steps:
         itworkbasis = transpose(workbasis * (bworkinvvar # replicate(1., nfit_this)))
         icovar = itworkbasis # workbasis
         beta = itworkbasis # workimage
         workidx = lindgen(nfit_this)
         extractivar = icovar[workidx,workidx]
         la_choldc, icovar, status=cholstat
         if (cholstat eq 0) then begin
            workcoeffs = la_cholsol(icovar, beta)
         endif else begin
            workcoeffs = replicate(0., nfit_this)
            extractivar[*] = 0.
         endelse
         workmodel = workbasis # workcoeffs

; Unpack results, accounting for masked fibers:
         fullcoeffs = replicate(0., nperbun + npoly)
         fullexivar = replicate(0., nperbun + npoly)
         fullcoeffs[use_component] = workcoeffs
         fullexivar[use_component] = extractivar

; Populate the model:
         ymodel[jmin[ibun]:jmax[ibun]] = workmodel
         ans[t_lo[ibun]:t_hi[ibun]] = fullcoeffs[0:nperbun-1]
         p[t_lo[ibun]:t_hi[ibun]] = sqrt(fullexivar[0:nperbun-1])
endif

endfor

; Apply fullreject mask to output flux and inverse-error:
p = p * (fullreject eq 0)
ans = ans * (fullreject eq 0)

; DISABLING "INPUTANS":
;#      if (keyword_set(inputans)) then begin
;#         if (ntrace*ncoeff NE n_elements(inputans)) then $
;#          message, 'Number of elements in INPUTANS is not equal to NTRACE*NCOEFF'
;#         ans[0:ntrace*ncoeff-1] = inputans
;#      endif

; REMOVING REFERENCE TO "NBAND":
;#      if (proftype EQ 10 AND nband LT 3) then $
;#         message, 'Do we have enough band-width for proftype=10?'

; I DON'T THINK WE NEED THIS:
;#      sigmal = 15.0

; REMOVING THE C CALL:
;#      result = call_external(soname, 'extract_row',$
;#       LONG(nx), FLOAT(xvar), FLOAT(fimage), workinvvar, ymodel, LONG(ntrace), $
;#       LONG(npoly), FLOAT(xcen), FLOAT(sigma), LONG(proftype), $
;#       FLOAT(reject), partial, fullreject, qcovar, LONG(squashprofile), $
;#       FLOAT(whopping), LONG(whoppingct), FLOAT(wsigma), $
;#       ncoeff, LONG(nband), ma, ans, wfixarr, p, fscat, fullcovar, $
;#       sigmal, contribution)


      diffs = (fimage - ymodel) * sqrt(workinvvar) 
      chisq = total(diffs*diffs)
      countthese = total(wfixarr)
      reducedChi = chisq / ((total(mask) - countthese) > 1)
      errscale = 1.0
      if (relative) then errscale = sqrt((chisq/total(mask)) > 1.0)

      finished = 1
      ; I'm changing the rejection algorithm
      if (NOT keyword_set(oldreject)) then begin
         
         goodi = where(workinvvar GT 0)
         if (goodi[0] NE -1) then begin 
           indchi = diffs[goodi] / (highrej*errscale)
           neg = where(indchi LT 0)
           if (neg[0] NE -1) then $
             indchi[neg] = -indchi[neg] * highrej / lowrej

           badcandidate = where(indchi GT 1.0, badct)
           if (badct GT 0) then begin
             finished = 0
             ; find groups
             padbad = [-2, badcandidate, nx + 2]
             startgroup = where(padbad[1:badct] - padbad[0:badct-1] GT 1, ngrp)
             endgroup = where(padbad[2:badct+1] - padbad[1:badct] GT 1)
             for i=0, ngrp - 1 do begin
               worstdiff = max(indchi[startgroup[i]:endgroup[i]],worst)
               mask[goodi[badcandidate[worst+startgroup[i]]]] = 0
	     endfor
             totalreject = totalreject + ngrp
           endif
         endif
      endif else begin
        badhigh = where(diffs GT highrej*errscale, badhighct)
        badlow = where(diffs LT -lowrej*errscale, badlowct)
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
      endelse

      niter = niter + 1
      if (niter EQ maxiter) then finished = 1



   endwhile

   ;----------
   ; Add bits to PIXELMASK

   pixelmask = pixelmask OR (pixelmask_bits('PARTIALREJECT') * fix(partial))
   pixelmask = pixelmask OR (pixelmask_bits('FULLREJECT') * fix(fullreject))

   return, ans
end

