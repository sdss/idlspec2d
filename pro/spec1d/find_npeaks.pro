; This used to be used by ZCOMPUTE, but is no longer.
function find_npeaks, yflux, xvec, nfind=nfind, minsep=minsep, $
 width=width, ypeak=ypeak, xerr=xerr, npeak=npeak

   ndata = n_elements(yflux)
   if (ndata EQ 1) then $
    return, 0

   if (NOT keyword_set(xvec)) then xvec = lindgen(ndata)
   if (NOT keyword_set(nfind)) then nfind = 1
   if (n_elements(minsep) EQ 0) then minsep = 0
   if (NOT keyword_set(width)) then width = 3
   if (xvec[1] GT xvec[0]) then isign = 1 $ ; ascending X
    else isign = -1 ; descending X

   ;----------
   ; Make a copy of YFLUX for finding local maxima; this will be modified
   ; each time a peak is found by filling with values of YDONE where we
   ; are no longer allowed to search.

   ycopy = yflux
   yderiv = [ycopy[1:ndata-1] - ycopy[0:ndata-2], 0]
   ydone = min(ycopy)

   ;----------
   ; Find up to NFIND peaks

   for ifind=0, nfind-1 do begin

      ;----------
      ; Locate next maximum

      junk = max(ycopy, imax)

      ;----------
      ; Centroid on this peak

      xpeak1 = zfitmax(yflux, xvec, xguess=xvec[imax], width=width, $
       xerr=xerr1, ypeak=ypeak1)

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
       AND shift(yderiv,1) LT 0)))[0]
      if (ix1 EQ -1) then ix1 = 0
      ix2 = (where(isign*xvec GT (isign*xvec[imax] + minsep) AND yderiv GT 0))[0]
      if (ix2 EQ -1) then ix2 = ndata-1

      ycopy[ix1:ix2] = ydone
;print,xpeak1,ix1,ixc,ix2, xvec[imax]
;splot,yflux
;soplot,ycopy,color='red'

      ;----------
      ; Test to see if we can find any more peaks

      junk = where(ycopy GT ydone, ct)
      if (ct EQ 0) then ifind = nfind

   endfor

   npeak = n_elements(xpeak)

   return, xpeak
end
