;------------------------------------------------------------------------------
; Fit the maximum of YARR
function zfitmax, yarr, xarr, xguess=xguess, width=width, xerr=xerr, ypeak=ypeak

   npts = n_elements(yarr)
   if (NOT keyword_set(xarr)) then xarr = findgen(npts)
   if (NOT keyword_set(xguess)) then begin
      ypeak = max(yarr, imax)
      xguess = xarr[imax]
   endif else begin
      junk = min(abs(xarr - xguess), indx)
      ypeak = yarr[indx[0]]
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
   indx = where(xarr GE xleft AND xarr LE xright, ct)
   if (ct LT 3) then $
    return, xguess

;   if (!version.release LT '5.4') then begin
;      coeff = poly_fit(xarr[indx]-xguess, yarr[indx], 2, $
;       yfit, junk, junk, corrm)
;   endif else begin
;      coeff = poly_fit(xarr[indx]-xguess, yarr[indx], 2, $
;       yfit=yfit, covar=corrm, /double)
;   endelse
;   xbest = -0.5 * coeff[1] / coeff[2] + xguess
;   xerr = - 0.5 * sqrt(corrm[1,1]) / coeff[2] $
;    + 0.5 * coeff[1] * sqrt(corrm[2,2]) / (coeff[2])^2

   ; Replace above with SVDFIT, and it should work with IDL 5.3 or 5.4.
   ; The above code gave wrong XERR from IDL 5.4 due to changes in POLY_FIT().
   ndegree = 3
   coeff = svdfit(xarr[indx]-xguess, yarr[indx], ndegree, $
    yfit=yfit, sigma=corrsig, /double)
   yerror = sqrt(total( (yarr[indx]-yfit)^2 / (n_elements(indx)-ndegree) ))
   xbest = -0.5 * coeff[1] / coeff[2] + xguess
   xerr = - 0.5 * corrsig[1] / coeff[2] $
    + 0.5 * coeff[1] * corrsig[2] / (coeff[2])^2
   xerr = xerr * yerror ; Rescale by the apparent errors in Y
   ; Would be equiv. to SVDFIT(WEIGHTS=REPLICATE(1.,N_ELEMENTS(INDX))/YERROR)

   ; Insist that XBEST is a maximum (not a minimum)
   if (coeff[2] GE 0) then $
    return, xguess

   ; Insist that the maximum is in the fitting range, and not out of bounds
   if (xbest LT xleft OR xbest GT xright) then $
    return, xguess

;   psigma = sqrt( total((xarr[indx] - xbest)^2 * yfit) / total(yfit) ) ; ???
   ypeak = poly(xbest-xguess, coeff)
;xtmp = min(xarr[indx]) + 0.1*findgen(50)
;plot,xarr[indx], yarr[indx],/yno,ps=-1
;djs_oplot, xtmp, poly(xtmp-xguess,coeff), color='red'
;djs_oplot, [xbest], [ypeak], ps=4,color='red'

   return, xbest
end
;------------------------------------------------------------------------------
