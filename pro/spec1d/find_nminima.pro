
forward_function mpfit, mpfitfun, mpfitpeak, mpfitpeak_gauss, $
  mpfitpeak_lorentz, mpfitpeak_moffat, mpfitpeak_u

;------------------------------------------------------------------------------
; Fit the minimum of YARR with a quadratic or gaussian.
; Return value is the minimum value of chi^2/DOF.

function zfitmin, yarr, xarr, dofarr=dofarr, $
 xguess=xguess, width=width, xerr=xerr, ypeak=ypeak, errcode=errcode, $
 doplot=doplot, _EXTRA=KeywordsForPlot

   npts = n_elements(yarr)
   if (NOT keyword_set(xarr)) then xarr = findgen(npts)
   if (keyword_set(dofarr)) then ydof = yarr / (dofarr + (dofarr EQ 0)) $
    else ydof = yarr
   if (NOT keyword_set(xguess)) then begin
      junk = min(ydof, imin)
      xguess = xarr[imin]
      ypeak = ydof[imin]
   endif else begin
      junk = min(abs(xarr - xguess), imin)
      ypeak = ydof[imin]
   endelse

   ; Set default return values
   errcode = 0L
   xerr = 0.0

   ; Insist that there be at least 1 point to the left and right of XGUESS.
   junk = where(xarr LT xguess, nleft)
   junk = where(xarr GT xguess, nright)
   if (nleft EQ 0 OR nright EQ 0) then begin
      errcode = -1L
   endif

   if (keyword_set(width)) then begin
      xleft = xguess - width
      xright = xguess + width
   endif else begin
      xleft = min(xarr)
      xright = max(xarr)
   endelse
   indx = where(xarr GE xleft AND xarr LE xright, nthis)
   if (nthis LT 3 AND errcode EQ 0) then begin
      errcode = -2L
   endif

   ; Sort by X, which is necessary for the MPFITPEAK routine.
   ; Note that we always expect at least one point, which is
   ; at XGUESS.
   indx = indx[sort(xarr[indx])]
   thisx = xarr[indx]
   meandof = mean(dofarr[indx])
   thisy = ydof[indx] * meandof

   if (keyword_set(doplot)) then begin
      xplot = thisx[0] + findgen(101) * (thisx[nthis-1] - thisx[0]) / 100.
   endif

   ;----------
   ; Case of more than 3 points: Gaussian fit

   if (nthis GT 3 AND errcode EQ 0) then begin

      nterms = 4
      yfit = mpfitpeak(thisx-xguess, thisy, coeff, nterms=nterms, $
       /gaussian, /negative, perror=perror)
      if (nthis LE nterms) then $
       yerror = 0 $
      else $
       yerror = sqrt(total( (thisy-yfit)^2 / (nthis - nterms) ))

      ; Compute the fit error of the minimum of the quadratic.
      ; We rescale by the apparent errors in Y.
      xerr1 = perror[1] * yerror

      ; Compute where chi^2 increases by 1.
      ; Insist that the gaussian fit spans a range of at least one
      ; in the Y-axis, such that we can compute the formal error
      ; where chi^2 increases by 1.  The following also applies the
      ; constraint that XBEST is a minimum (not a maximum).
      if (coeff[0] LT -1.) then begin
         xerr2 = coeff[2] * sqrt(2. * alog(coeff[0]/(coeff[0]+1.)))
      endif else begin
         errcode = -5L
      endelse

      xbest = coeff[1] + xguess
      ybest = coeff[0]
      for ic=3, nterms-1 do $
       ybest = ybest + coeff[ic] * coeff[1]^(ic-3)
      ybest = ybest / meandof

      ; Insist that XBEST is a minimum (not a maximum)
      if (coeff[0] GT 0) then begin
         errcode = -4L
      endif

      if (keyword_set(doplot)) then $
       yplot = mpfitpeak_gauss(xplot - xguess, coeff)

   ;----------
   ; Case of exactly 3 points: Quadractic fit

   endif else if (nthis EQ 3 AND errcode EQ 0) then begin

      ndegree = 3
      coeff = svdfit(thisx-xguess, thisy, ndegree, $
       yfit=yfit, covar=covar, sigma=corrsig, /double)
      if (nthis LE ndegree) then $
       yerror = 0 $
      else $
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

      ybest = poly(xbest-xguess, coeff) / meandof

      ; Insist that XBEST is a minimum (not a maximum)
      if (coeff[2] LT 0) then begin
         errcode = -3L
      endif

      if (keyword_set(doplot)) then begin
         yplot = coeff[0]
         for ic=1, ndegree-1 do yplot = yplot + coeff[ic] * (thisx - xguess)^ic
      endif

   endif

   if (keyword_set(xerr1) AND keyword_set(xerr2)) then $
    xerr = sqrt(xerr1^2 + xerr2^2)

   ; Insist that the minimum is in the fitting range, and not out of bounds
   if (keyword_set(xbest) AND errcode EQ 0) then begin
      if (xbest LT xleft OR xbest GT xright) then begin
         errcode = -6L
      endif
   endif

   if (errcode EQ 0) then begin
      ypeak = ybest
   endif else begin
      xbest = xguess
   endelse

   if (keyword_set(doplot)) then begin
      !x.ticks = 2
      !x.tickv = [thisx[0], xbest, thisx[nthis-1]]
      !x.tickname = [' ', strtrim(string(xbest),2), ' ']
      djs_plot, [thisx], [thisy]/meandof, psym=-4, $
       _EXTRA=KeywordsForPlot
      if (errcode EQ 0) then color = 'green' $
       else color='red'
      if (keyword_set(yplot)) then $
       djs_oplot, [xplot], [yplot]/meandof, color=color
      djs_oplot, [xbest], [ypeak], psym=2, symsize=2, color=color
      djs_oplot, [xbest,xbest]-xerr, !y.crange, linestyle=2, color=color
      djs_oplot, [xbest,xbest]+xerr, !y.crange, linestyle=2, color=color
   endif

   return, xbest
end

;------------------------------------------------------------------------------
function find_nminima, yflux, xvec, dofarr=dofarr, nfind=nfind, minsep=minsep, $
 width=width, ypeak=ypeak, xerr=xerr, errcode=errcode, npeak=npeak, $
 plottitle=plottitle, doplot=doplot

   ndata = n_elements(yflux)
   if (ndata EQ 1) then $
    return, 0

   if (NOT keyword_set(xvec)) then xvec = lindgen(ndata)
   if (NOT keyword_set(dofarr)) then dofarr = 1
   if (NOT keyword_set(nfind)) then nfind = 1
   if (n_elements(minsep) EQ 0) then minsep = 0
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
   ; Set up for plots

   if (keyword_set(doplot)) then begin
      bangp = !p
      bangx = !x
      bangy = !y
      dxplot = 0.85 / nfind
      !p.position = [0.10, 0.10, 0.10+dxplot, 0.45]
;      !x.margin = [0,0]
;      !x.omargin = [10,3]
      !y.range = minmax(ycopy)
      !y.title = textoidl('\chi^2/DOF')
      !p.multi = [0,nfind+1,1]
      !p.charsize = 1.5
      !x.charsize = 1.5
      !y.charsize = 1.5
   endif

   ;----------
   ; Find up to NFIND peaks

   for ifind=0, nfind-1 do begin

      ;----------
      ; Locate next minimum

      junk = min(ycopy, imin)

      ;----------
      ; Centroid on this peak (local minimum)

      if ((keyword_set(doplot)) $
       AND ifind GT 0) then begin
         !y.tickname = replicate(' ',30)
         !y.title = ''
         !p.position[[0,2]] = !p.position[[0,2]] + dxplot
      endif

      xpeak1 = zfitmin(yflux, xvec, dofarr=dofarr, xguess=xvec[imin], $
       width=width, xerr=xerr1, ypeak=ypeak1, errcode=errcode1, doplot=doplot)

      ;----------
      ; Save return values

      if (ifind EQ 0) then begin
         xpeak = xpeak1
         xerr = xerr1
         errcode = errcode1
         ypeak = ypeak1
      endif else begin
         xpeak = [xpeak, xpeak1]
         xerr = [xerr, xerr1]
         errcode = [errcode, errcode1]
         ypeak = [ypeak, ypeak1]
      endelse

      ;----------
      ; Exclude from future peak-finding all points within MINSEP of this
      ; peak, up until the function is increasing again.

      junk = min(abs(xvec - xvec[imin]), ixc)
      ix1 = (reverse(where(isign*xvec LT (isign*xvec[imin] - minsep) $
       AND shift(yderiv,1) GT 0)))[0]
      if (ix1 EQ -1) then ix1 = 0
      ix2 = (where(isign*xvec GT (isign*xvec[imin] + minsep) AND yderiv LT 0))[0]
      if (ix2 EQ -1) then ix2 = ndata-1

      ycopy[ix1:ix2] = ydone

      ;----------
      ; Test to see if we can find any more peaks

      junk = where(ycopy LT ydone, ct)
      if (ct EQ 0) then ifind = nfind

   endfor

   npeak = n_elements(xpeak)

   if (keyword_set(doplot)) then begin
      !p.position = [0.10, 0.55, 0.95, 0.95]
      xcsize = !x.charsize
      ycsize = !y.charsize
      yrange = !y.range
      !x = bangx
      !y = bangy
      yplot = yflux
      if (keyword_set(dofarr)) then yplot = yplot / dofarr
      djs_plot, xvec, yplot, yrange=yrange, $
       xcharsize=xcsize, ycharsize=ycsize, yrange=yrange, $
       title=plottitle, ytitle='\chi^2/DOF'
      for ipeak=0, npeak-1 do begin
         if (errcode[ipeak] EQ 0) then color = 'green' $
          else color = 'red'
         djs_oplot, [xpeak[ipeak]], [ypeak[ipeak]], psym=2, color=color
      endfor
      !p = bangp
      !p.color = djs_icolor('default')
   endif

   return, xpeak
end
;------------------------------------------------------------------------------
