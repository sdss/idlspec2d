
;------------------------------------------------------------------------------
; Fit the minimum of YARR with a quadratic or gaussian.
function zfitmin, yarr, xarr, dofarr=dofarr, $
 xguess=xguess, width=width, xerr=xerr, ypeak=ypeak, doplot=doplot

   npts = n_elements(yarr)
   if (NOT keyword_set(xarr)) then xarr = findgen(npts)
   if (keyword_set(dofarr)) then ydof = yarr / (dofarr + (dofarr EQ 0)) $
    else ydof = yarr
   if (NOT keyword_set(xguess)) then begin
      ypeak = max(ydof, imax)
      xguess = xarr[imax]
   endif else begin
      junk = min(abs(xarr - xguess), indx)
      ypeak = ydof[indx[0]]
   endelse
   if (NOT keyword_set(width)) then width = 1

   ; Set return values in the event of a bad fit
   xerr = 0.0

   ; Insist that there be at least 1 point to the left and right of XGUESS.
   junk = where(xarr LT xguess, nleft)
   junk = where(xarr GT xguess, nright)
   if (nleft EQ 0 OR nright EQ 0) then $
    return, xguess

   xleft = xguess - width
   xright = xguess + width
   indx = where(xarr GE xleft AND xarr LE xright, nthis)
   if (nthis LT 3) then $
    return, xguess

   ; Sort by X, which is necessary for the MPFITPEAK routine.
   indx = indx[sort(xarr[indx])]
   thisx = xarr[indx]
   thisy = ydof[indx] * mean(dofarr[indx])

   ;----------
   ; Case of exactly 3 points: Quadractic fit

   if (nthis EQ 3) then begin

      ndegree = 3
      coeff = svdfit(thisx-xguess, thisy, ndegree, $
       yfit=yfit, covar=covar, sigma=corrsig, /double)
      yerror = sqrt(total( (thisy-yfit)^2 / (nthis - ndegree) ))
      xbest = -0.5 * coeff[1] / coeff[2] + xguess

      ; Compute the fit error of the minimum of the quadratic.
      ; We rescale by the apparent errors in Y, which would be equivalent
      ; to the call SVDFIT(WEIGHTS=REPLICATE(1.,N_ELEMENTS(INDX))/YERROR)
      dx0_db = -0.5 / coeff[2]
      dx0_dc = 0.5 * coeff[1] / (coeff[2])^2
      xerr1 = sqrt( dx0_db^2 * covar[1,1] + dx0_dc^2 * covar[2,2] $
       + 2 * dx0_db * dx0_dc * covar[1,2] ) * yerror

      ; Compute where chi^2 increases by 1
      xerr2 = 1 / sqrt(coeff[2])

      ypeak = poly(xbest-xguess, coeff)

      ; Insist that XBEST is a minimum (not a maximum)
      if (coeff[2] LT 0) then $
       return, xguess

   ;----------
   ; Case of more than 3 points: Gaussian fit

   endif else begin

      nterms = 4
      yfit = mpfitpeak(thisx-xguess, thisy, coeff, nterms=nterms, $
       /gaussian, /negative, perror=perror)
      yerror = sqrt(total( (thisy-yfit)^2 / (nthis - nterms) ))

      ; Compute the fit error of the minimum of the quadratic.
      ; We rescale by the apparent errors in Y.
      xerr1 = coeff[1] * yerror

      ; Compute where chi^2 increases by 1.
      ; Insist that the gaussian fit spans a range of at least one
      ; in the Y-axis, such that we can compute the formal error
      ; where chi^2 increases by 1.
      if (coeff[0] LT -1.) then $
       xerr2 = coeff[2] * sqrt(2. * alog(coeff[0]/(coeff[0]+1.))) $
      else $
       return, xguess

      xbest = coeff[1] + xguess
      ypeak = coeff[0]
      for ic=3, nterms-1 do $
       ypeak = ypeak + coeff[ic] * coeff[1]^(ic-3)

      ; Insist that XBEST is a minimum (not a maximum)
      if (coeff[0] LT 0) then $
       return, xguess

   endelse

   xbesterr = sqrt(xerr1^2 + xerr2^2)

   ; Insist that the minimum is in the fitting range, and not out of bounds
   if (xbest LT xleft OR xbest GT xright) then $
    return, xguess

   if (keyword_set(doplot)) then begin
      djs_plot, thisx, thisy, /yno
      djs_oplot, thisx, yfit, color='green'
      djs_oplot, [xbest], [ypeak], psym=4, symsize=2, color='green'
      djs_oplot, [xbest,xbest]-xerr, !y.crange, linestyle=2, color='green'
      djs_oplot, [xbest,xbest]+xerr, !y.crange, linestyle=2, color='green'
   endif

   xerr = xbesterr
   return, xbest
end
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
function find_nminima, yflux, xvec, dofarr=dofarr, nfind=nfind, minsep=minsep, $
 width=width, ypeak=ypeak, xerr=xerr, npeak=npeak

   ndata = n_elements(yflux)
   if (ndata EQ 1) then $
    return, 0

   if (NOT keyword_set(xvec)) then xvec = lindgen(ndata)
   if (NOT keyword_set(dofarr)) then dofarr = 1
   if (NOT keyword_set(nfind)) then nfind = 1
   if (n_elements(minsep) EQ 0) then minsep = 0
   if (NOT keyword_set(width)) then width = 5
   if (xvec[1] GT xvec[0]) then isign = 1 $ ; ascending X
    else isign = -1 ; descending X

   ;----------
   ; Make a copy of YFLUX/DOFARR for finding local minima; this will be modified
   ; each time a peak is found by filling with values of YDONE where we
   ; are no longer allowed to search.

   ycopy = yflux / (dofarr + (dofarr EQ 0))
   yderiv = [ycopy[1:ndata-1] - ycopy[0:ndata-2], 0]
   ydone = max(ycopy)

   ;----------
   ; Find up to NFIND peaks

   for ifind=0, nfind-1 do begin

      ;----------
      ; Locate next minimum

      junk = min(ycopy, imax)

      ;----------
      ; Centroid on this peak (local minimum)

      xpeak1 = zfitmin(yflux, xvec, dofarr=dofarr, xguess=xvec[imax], $
       width=width, xerr=xerr1, ypeak=ypeak1)

      ;----------
      ; Save return values

      if (ifind EQ 0) then begin
         xpeak = xpeak1
         xerr = xerr1
         ypeak = ypeak1
      endif else begin
         xpeak = [xpeak, xpeak1]
         xerr = [xerr, xerr1]
         ypeak = [ypeak, ypeak1]
      endelse

      ;----------
      ; Exclude from future peak-finding all points within MINSEP of this
      ; peak, up until the function is increasing again.

      junk = min(abs(xvec - xvec[imax]), ixc)
      ix1 = (reverse(where(isign*xvec LT (isign*xvec[imax] - minsep) $
       AND shift(yderiv,1) GT 0)))[0]
      if (ix1 EQ -1) then ix1 = 0
      ix2 = (where(isign*xvec GT (isign*xvec[imax] + minsep) AND yderiv LT 0))[0]
      if (ix2 EQ -1) then ix2 = ndata-1

      ycopy[ix1:ix2] = ydone

      ;----------
      ; Test to see if we can find any more peaks

      junk = where(ycopy LT ydone, ct)
      if (ct EQ 0) then ifind = nfind

   endfor

   npeak = n_elements(xpeak)

   return, xpeak
end
;------------------------------------------------------------------------------
