;+
; NAME:
;   trace320crude
;
; PURPOSE:
;   Calling script to return 320 full traces using TRACE_CRUDE.
;
; CALLING SEQUENCE:
;   xset = trace320crude( fimage, invvar, [mthresh=, ystart=, nmed=, $
;    xmask=, yset=, maxerr=, maxshifte=, maxshift0=, xerr=, maxdev=, ngrow=, $
;    fibermask=  ] )
;
; INPUTS:
;   fimage     - Image
;
; OPTIONAL INPUTS FOR TRACE320CEN:
;   mthresh    - Threshold for peak-finding in convolved row; default to 0.5
;                times the dispersion (found with djs_iterstat).
;   ystart     - Y position in image to search for initial X centers; default
;                to the central row
;   nmed       - Number of rows to median filter around YSTART; default to 21
;
; OPTIONAL INPUTS FOR TRACE_CRUDE:
;   invvar     - Inverse variance (weight) image
;   radius     - Radius for centroiding; default to 3.0
;   maxerr     - Maximum error in centroid allowed for valid recentering;
;                default to 0.2
;   maxshifte  - Maximum shift in centroid allowed for valid recentering;
;                default to 0.1
;   maxshift0  - Maximum shift in centroid allowed for initial row;
;                default to 0.5
;
; OPTIONAL INPUTS:
;   maxdev     - Maximum deviation of X in pixels; default to rejecting any
;                XPOS positions that deviate by more than 1.0 pixels from
;                a polynomial remapping of the centroids from other rows.
;   ngrow      - For each trace, replace all centroids within NGROW rows
;                of a bad centroid with the predicted centroid locations.
;                Default to 5.
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;
; OUTPUTS:
;   xset       - X centers for all traces
;
; OPTIONAL OUTPUTS:
;   yset       - Y centers for all traces
;   xerr       - Errors for XSET
;   xmask      - Mask set to 1 for good fiber centers, 0 for bad;
;                same dimensions as XSET.
;   fibermask  - (Modified.)
;
; COMMENTS:
;   Without djs_maskinterp, hot columns skew traces 
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_maskinterp()
;   fibermask_bits()
;   trace_crude()
;   trace_fweight()
;   trace320cen()
;
; REVISION HISTORY:
;   13-Sep-1999  Written by David Schlegel, Princeton.
;    8-Jul-2001  Added djs_maskinterp call
;-
;------------------------------------------------------------------------------
function trace320crude, image, invvar, ystart=ystart, nmed=nmed, $
   xmask=xmask, radius=radius, yset=yset, maxerr=maxerr, maxshifte=maxshifte, $
   maxshift0=maxshift0, xerr=xerr, maxdev=maxdev, ngrow=ngrow, $
   fibermask=fibermask

   if (NOT keyword_set(maxdev)) then maxdev = 1.0
   if (NOT keyword_set(ngrow)) then ngrow = 5

   ;----------
   ; If INVVAR is set, then start by interpolating over bad pixels

   if (keyword_set(invvar)) then $
    fimage = djs_maskinterp(image, (invvar LE 0), iaxis=0) $
   else $
    fimage = image

   ;----------
   ; Find the 320 X-centers in the row specified by YSTART

   ; XGOOD=1 for fibers that were actually found, 0 otherwise
   xstart = trace320cen(fimage, mthresh=mthresh, ystart=ystart, nmed=nmed, $
    xgood=xgood)
   ntrace = N_elements(xstart) ; Better be 320
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(ntrace)

   ;----------
   ; Trace those 320

   xset = trace_crude(fimage, invvar, xstart=xstart, ystart=ystart, $
    radius=radius, yset=yset, maxerr=maxerr, maxshifte=maxshifte, $
    maxshift0=maxshift0, xerr=xerr)
   xmask = xerr LT 990  ; =1 for good centers, =0 for bad

   ;----------
   ; Compare the traces in each row to those in row YSTART.
   ; Our assumption is that those centers should be a polynomial mapping
   ; of the centers from row YSTART.  Centers that are deviant from this
   ; mapping are replaced with the position predicted by this mapping.

   ny = (size(fimage, /dimens))[1]
   ndegree = 4 ; Five terms

   ;----------
   ; Loop to find all deviant centroids, and add these to the mask XMASK.
   ; XMASK=1 for good.

   for iy=0, ny-1 do begin
      xcheck = xgood AND xmask[iy,*] ; Test for good fiber & good centroid
      if (total(xcheck) GT ndegree+1) then begin
         if (!version.release LT '5.4') then begin
            coeff = polyfitw(xstart, xset[iy,*], xcheck, ndegree, xfit)
         endif else begin
            coeff = polyfitw(xstart, xset[iy,*], xcheck, ndegree, xfit, /double)
         endelse

        xdiff = xfit - xset[iy,*]
        ibad = where(abs(xdiff) GT maxdev)
        if (ibad[0] NE -1) then xmask[iy,ibad] = 0
      endif else begin
        xmask[iy,*] = 0 ; Too few good centroids in this row; mark all as bad
      endelse
   endfor

   ;----------
   ; Smooth the bad centroids to NGROW adjacent rows (of the same trace)

   for itrace=0, ntrace-1 do begin
      xmask[*,itrace] = smooth( xmask[*,itrace]+0.0, 2*ngrow+1) EQ 1
   endfor

   ;----------
   ; Loop to fix deviant centroids

   for iy=0, ny-1 do begin
      ixbad = where(xmask[iy,*] EQ 0, nbad)
      if (nbad GT 0 AND nbad LT ntrace-ndegree) then begin
         ixgood = where(xmask[iy,*] EQ 1)

         if (!version.release LT '5.4') then begin
            coeff = polyfitw(xstart, xset[iy,*], xmask[iy,*], $
             ndegree, xfit)
         endif else begin
            coeff = polyfitw(xstart, xset[iy,*], xmask[iy,*], $
             ndegree, xfit, /double)
         endelse

         xset[iy,ixbad] = xfit[ixbad]
      endif
   endfor

   ;----------
   ; Perform a second centering iteration on the fibers initially rejected
   ; by TRACE320CEN.  Those fibers might not actually be bad, but might
   ; have just had bad pixels near YSTART.

   indx = where(xgood EQ 0, ct)
   for ii=0, ct-1 do begin
      itrace = indx[ii]

      tmp_xpos = trace_fweight(fimage, xset[*,itrace], yset[*,itrace], $
       radius=radius, xerr=tmp_xerr, invvar=invvar)

      xset[*,itrace] = tmp_xpos
      xerr[*,itrace] = tmp_xerr
      xmask[*,itrace] = tmp_xerr LT 990 ; =1 for good centers, =0 for bad
   endfor

   ;----------
   ; Replace XSET with a smooth trace-set

;   xy2traceset, yset, xset, tset, ncoeff=5, yfit=xnew, invvar=xmask, $
;    maxdev=maxdev, maxrej=1, /sticky
;   xset = xnew

   ;----------
   ; Set FIBERMASK bit for any fiber with more than 20% of its positions
   ; masked, which includes any positions off the CCD.  Do not pay any
   ; attention to XGOOD, since that may indicate that a fiber is only
   ; bad near YSTART.

   ibad = where(total(1-xmask, 1) GT 0.20*ny)
   if (ibad[0] NE -1) then $
    fibermask[ibad] = fibermask[ibad] OR fibermask_bits('BADTRACE')

   return, xset
end
;------------------------------------------------------------------------------
