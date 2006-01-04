;+
; NAME:
;   find_nminima
;
; PURPOSE:
;   Find one or several minima in a vector of chi^2 values.
;
; CALLING SEQUENCE:
;   xpeak = find_nminima( yflux, [ xvec, dofarr=, nfind=, minsep=, $
;    width=, ypeak=, xerr=, errcode=, npeak=, plottitle=, /doplot, /debug ]
;
; INPUTS:
;   yflux          - Y values
;
; OPTIONAL INPUTS:
;   xvec           - X values, which must either be sorted in ascending
;                    or descending order; default to 0-indexed integers.
;   dofarr         - If set, then fit to the minima in the function
;                    YFLUX/DOFARR, but avoiding any points where DOFARR
;                    is set to zero.
;   nfind          - Number of minima to find; default to 1.  It is possible
;                    to find fewer than NFIND minima.
;   minsep         - Minimum separation between local minima.  If a peak
;                    is found closer than MINSEP to an existing peak, then
;                    the latter peak is discarded.
;   width          - Width to use when selecting the points used in the fit.
;                    Only use points where XVEC is within WIDTH of the
;                    the lowest-values point (which is used as the initial
;                    guess); default to using all points.
;
; OUTPUTS:
;   ypeak          - Fit value for either chi^2 or chi^2/DOF at the minima.
;                    If the fit value is less than zero, then change it to zero.
;
; OPTIONAL OUTPUTS:
;   xerr           - Formal errors of XPEAK.
;   errcode        - Error codes for each minima; 0 for no errors in the fit.
;   npeak          - The number of peaks found, between [0,NFIND].
;   plottitle      - Title of plot (if /DOPLOT is set).
;   doplot         - If set, then make plots.  Discarded peaks are not plotted.
;   debug          - If set, then wait for keystroke after plot.
;
; COMMENTS:
;   This routine calls SVDFIT for fitting quadratics, or MPFIT for
;   fitting gaussians.
;   
; EXAMPLES:
;   
; BUGS:
;
; PROCEDURES CALLED:
;   djs_icolor()
;   djs_oplot
;   djs_plot
;   mpfitpeak
;   mpfitpeak_gauss
;   textoidl
;
; INTERNAL SUPPORT ROUTINES:
;   zfitmin()
;  
; REVISION HISTORY:
;   22-Aug-2001  Written by D. Schlegel, Princeton 
;   30-Jul-2002  Fixed bug - yrange passed twice to djs_plot
;-  
;------------------------------------------------------------------------------

forward_function mpfit, mpfitfun, mpfitpeak, mpfitpeak_gauss, $
  mpfitpeak_lorentz, mpfitpeak_moffat, mpfitpeak_u

;------------------------------------------------------------------------------
; Fit the minimum of YARR with a quadratic or gaussian.
; Return value is the minimum value of chi^2/DOF.
; Set XPLOTFIT to 1 in order to return values in XPLOTFIT,YPLOTFIT,
; XPLOTVAL,YPLOTVAL for plotting.

function zfitmin, yarr, xarr, dofarr=dofarr, $
 xguess=xguess, width=width, xerr=xerr, ypeak=ypeak, errcode=errcode, $
 xplotfit=xplotfit, yplotfit=yplotfit, xplotval=xplotval, yplotval=yplotval

   if (keyword_set(xplotfit)) then doplot = 1
   npts = n_elements(yarr)
   if (NOT keyword_set(xarr)) then xarr = dindgen(npts)
   if (keyword_set(dofarr)) then $
    ydof = double(yarr) / (dofarr + (dofarr EQ 0)) $
   else $
    ydof = double(yarr)
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
   xerr1 = 0.0
   xerr2 = 0.0
   xplotfit = 0
   yplotfit = 0

   ; Insist that there be at least 1 point to the left and right of XGUESS.
   junk = where(xarr LT xguess, nleft)
   junk = where(xarr GT xguess, nright)
   if (nleft EQ 0 OR nright EQ 0) then begin
      errcode = -1L
   endif

   ; Insist that not all the Y values are the same
   if (min(yarr) EQ max(yarr)) then begin
      errcode = -7L
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
   thisx = double(xarr[indx])
   meandof = mean(dofarr[indx])
   thisy = ydof[indx] * meandof

   if (keyword_set(doplot)) then begin
      xplotfit = thisx[0] + dindgen(101) * (thisx[nthis-1] - thisx[0]) / 100.
      xplotval = thisx
      yplotval = thisy / meandof
   endif

   ;----------
   ; Case of more than 3 points: Gaussian fit

   if (nthis GT 3 AND errcode EQ 0) then begin

      nterms = 4
      yfit = mpfitpeak(thisx-xguess, thisy, coeff, nterms=nterms, $
       /gaussian, /negative, perror=perror, status=status)
      if (nthis LE nterms) then $
       yerror = 0 $
      else $
       yerror = sqrt(total( (thisy-yfit)^2 / (nthis - nterms) ))

      ; Compute the fit error of the minimum of the quadratic.
      ; We rescale by the apparent errors in Y.
      if (status GT 0) then begin
         xerr1 = perror[1] * yerror
      endif else begin
         errcode = -8L
      endelse

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
       yplotfit = mpfitpeak_gauss(xplotfit - xguess, coeff) / meandof

   ;----------
   ; Case of exactly 3 points: Quadractic fit

   endif else if (nthis EQ 3 AND errcode EQ 0) then begin

      mmatrix = [[1d0,1d0,1d0], [thisx-xguess], [([thisx-xguess])^2]]
      coeff = invert(mmatrix) # thisy
      xbest = -coeff[1] / (2. * coeff[2]) + xguess
      ymin = coeff[0] - (coeff[1])^2 / (4. * coeff[2])

      ; Compute the fit error of the minimum of the quadratic.
      ; This is set to zero, since a quadratic goes exactly through the 3 pts.
      xerr1 = 0

      ; Compute where chi^2 increases by 1
      xerr2 = sqrt((coeff[1])^2 - 4.*coeff[2]*(coeff[0]-ymin-1.)) $
       / (2.*coeff[2])

      ybest = ymin / meandof

;      ndegree = 3
;      coeff = svdfit(thisx-xguess, thisy, ndegree, $
;       yfit=yfit, covar=covar, sigma=corrsig, /double)
;      if (nthis LE ndegree) then $
;       yerror = 0 $
;      else $
;       yerror = sqrt(total( (thisy-yfit)^2 / (nthis - ndegree) ))
;      xbest = -0.5 * coeff[1] / coeff[2] + xguess
;
;      ; Compute the fit error of the minimum of the quadratic.
;      ; We rescale by the apparent errors in Y, which would be equivalent
;      ; to the call SVDFIT(WEIGHTS=REPLICATE(1.,N_ELEMENTS(INDX))/YERROR)
;      dx0_db = -0.5 / coeff[2]
;      dx0_dc = 0.5 * coeff[1] / (coeff[2])^2
;      xerr1 = sqrt( dx0_db^2 * covar[1,1] + dx0_dc^2 * covar[2,2] $
;       + 2 * dx0_db * dx0_dc * covar[1,2] ) * yerror
;
;      ; Compute where chi^2 increases by 1
;      xerr2 = 1 / sqrt(coeff[2])
;
;      ybest = poly(xbest-xguess, coeff) / meandof

      ; Insist that XBEST is a minimum (not a maximum)
      if (coeff[2] LT 0) then begin
         errcode = -3L
      endif

      if (keyword_set(doplot)) then begin
         yplotfit = coeff[0]
         for ic=1, n_elements(coeff)-1 do $
          yplotfit = yplotfit + coeff[ic] * (xplotfit - xguess)^ic
         yplotfit = yplotfit / meandof
      endif

   endif

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

   ypeak = ypeak > 0 ; Insist that the Y minimum does not go below zero

   return, xbest
end

;------------------------------------------------------------------------------
function find_nminima, yflux, xvec, dofarr=dofarr, nfind=nfind, minsep=minsep, $
 width=width, ypeak=ypeak, xerr=xerr, errcode=errcode, npeak=npeak, $
 plottitle=plottitle, xtitle=xtitle, doplot=doplot, debug=debug

   ndata = n_elements(yflux)
   if (ndata EQ 1) then $
    return, 0

   if (NOT keyword_set(xvec)) then xvec = dindgen(ndata)
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
      !p.position = [0.15, 0.10, 0.10+dxplot, 0.45]
;      !x.margin = [0,0]
;      !x.omargin = [10,3]
;      !y.range = minmax(ycopy)
      if (keyword_set(chi2arr)) then $
       !y.title = textoidl('\chi^2/DOF') $
      else $
       !y.title = textoidl('\chi^2')
      oldmulti = !p.multi
      !p.multi = [0,nfind+1,1]
      !p.charsize = 1.5
      !x.charsize = 1.5
      !y.charsize = 1.5
   endif

   ;----------
   ; Find up to NFIND peaks

   npeak = 0
   while (npeak LT nfind) do begin

      ;----------
      ; Locate next minimum

      junk = min(ycopy, imin)

      ;----------
      ; Centroid on this peak (local minimum)

      if (keyword_set(doplot)) then xplotfit = 1
      xpeak1 = zfitmin(yflux, xvec, dofarr=dofarr, xguess=xvec[imin], $
       width=width, xerr=xerr1, ypeak=ypeak1, errcode=errcode1, $
       xplotfit=xplotfit, yplotfit=yplotfit, $
       xplotval=xplotval, yplotval=yplotval)

      ;----------
      ; Save return values

      if (npeak EQ 0) then begin
         ; Always retain the first peak
         xpeak = xpeak1
         xerr = xerr1
         errcode = errcode1
         ypeak = ypeak1
      endif else begin
         ; Retain this peak only if it is at least MINSEP from any
         ; peaks already found.
         if (min(abs(xpeak - xpeak1)) GT minsep) then begin
            xpeak = [xpeak, xpeak1]
            xerr = [xerr, xerr1]
            errcode = [errcode, errcode1]
            ypeak = [ypeak, ypeak1]
         endif else begin
            print, 'Discarding peak #', npeak+1, ' of ', nfind
         endelse
      endelse

      ; Increment NPEAK if a peak has been added
      noldpeak = npeak
      npeak = n_elements(xpeak)

      ;----------
      ; Make plot of this peak

      if (keyword_set(doplot) AND npeak GT noldpeak) then begin
         if (npeak GT 1) then begin
            !y.tickname = replicate(' ',30)
            !y.title = ''
            !p.position[[0,2]] = !p.position[[0,2]] + dxplot
         endif

         !x.ticks = 2
         nthis = n_elements(xplotval)
         !x.tickv = [xplotval[0], xpeak1, xplotval[nthis-1]]
         !x.tickname = [' ', strtrim(string(xpeak1),2), ' ']
         djs_plot, [xplotval], [yplotval], psym=4, /ynozero
         if (errcode1 EQ 0) then color = 'green' $
          else color='red'
         if (keyword_set(yplotfit)) then $
          djs_oplot, [xplotfit], [yplotfit], color=color
         djs_oplot, [xpeak1], [ypeak1], psym=2, symsize=2, color=color
         djs_oplot, [xpeak1,xpeak1]-xerr1, !y.crange, linestyle=2, color=color
         djs_oplot, [xpeak1,xpeak1]+xerr1, !y.crange, linestyle=2, color=color
      endif

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
      if (ct EQ 0) then npeak = nfind

   endwhile

   npeak = n_elements(xpeak)

   if (keyword_set(doplot)) then begin
      !p.position = [0.15, 0.55, 0.95, 0.95]
      xcsize = !x.charsize
      ycsize = !y.charsize
      yrange = !y.range
      !x = bangx
      !y = bangy
      yplot = yflux
      if (keyword_set(dofarr)) then yplot = yplot / dofarr
      djs_plot, xvec, yplot, yrange=yrange, psym=-4, /ynozero, $
       xcharsize=xcsize, ycharsize=ycsize, $
       title=plottitle, xtitle=xtitle
      for ipeak=0, npeak-1 do begin
         if (errcode[ipeak] EQ 0) then color = 'green' $
          else color = 'red'
         djs_oplot, [xpeak[ipeak]], [ypeak[ipeak]], psym=2, color=color
      endfor
      !p = bangp
      !p.color = djs_icolor('default')

      ; Wait for a keystroke...
      if (keyword_set(debug)) then begin
         print, 'Press any key...'
         cc = strupcase(get_kbrd(1))
      endif

      !p.multi = oldmulti
   endif

   return, xpeak
end
;------------------------------------------------------------------------------
