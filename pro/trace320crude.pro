;+
; NAME:
;   trace320crude
;
; PURPOSE:
;   Calling script to return 320 full traces using TRACE_CRUDE.
;
; CALLING SEQUENCE:
;   xset = trace320crude( fimage, invvar, [mthresh=, ystart=, nmed=, xgood=, $
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
;   fibermask  - Mask [NFIBER]
;
; OUTPUTS:
;   xset       - X centers for all traces
;
; OPTIONAL OUTPUTS:
;   yset       - Y centers for all traces
;   xerr       - Errors for XSET
;   xgood      - Set to 1 for fibers that were actually found, 0 otherwise
;   xmask      - Mask set to 1 for good fiber centers, 0 for bad;
;                same dimensions as XSET.
;   fibermask  - (Modified.)
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   fibermask_bits()
;   trace_crude()
;   trace320cen()
;
; REVISION HISTORY:
;   13-Sep-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function trace320crude, fimage, invvar, ystart=ystart, nmed=nmed, xgood=xgood, $
   xmask=xmask, radius=radius, yset=yset, maxerr=maxerr, maxshifte=maxshift, $
   maxshift0=maxshift0, xerr=xerr, maxdev=maxdev, ngrow=ngrow, $
   fibermask=fibermask

   if (NOT keyword_set(maxdev)) then maxdev = 1.0
   if (NOT keyword_set(ngrow)) then ngrow = 5

   ;----------
   ; Find the 320 X-centers in the row specified by YSTART

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
   ; Re-fit the centroids.  Do this in case one of the XSTART positions
   ; happenened to be bad.

;   xset = trace_fweight(fimage, xset, yset, radius=radius, invvar=invvar, $
;    xerr=xerr)
;;   xset = trace_gweight(fimage, xset, yset, sigma=1.0, invvar=invvar, $
;;    xerr=xerr)
;   xmask = xerr LT 990  ; =1 for good centers, =0 for bad

   ;----------
   ; Replace XSET with a smooth trace-set

;   xy2traceset, yset, xset, tset, ncoeff=5, yfit=xnew, invvar=xmask, $
;    maxdev=maxdev, /singlerej
;   xset = xnew

   ;----------
   ; Identify bad traces as those with more than 1/4 of their centroids bad

;   nx = (size(xset, /dimens))[0]
;   quarterbad = (total(xmask,1) LT 0.75 * nx)
;   ixgood = where(xgood AND NOT quarterbad)
;   xstart[ixgood] = xset[ystart,ixgood]

   ;----------
   ; Compare the traces in each row to those in row YSTART.
   ; Our assumption is that those centers should be a polynomial mapping
   ; of the centers from row YSTART.  Centers that are deviant from this
   ; mapping are replaced with the position predicted by this mapping.

   ny = (size(fimage, /dimens))[1]
   ndegree = 4 ; Five terms

   ;----------
   ; Loop to find all deviant centroids, and add these to the mask XMASK.

   for iy=0, ny-1 do begin
      xcheck = xgood AND xmask[iy,*] ; Test for good fiber & good centroid
      if (total(xcheck) GT ndegree+1) then begin
        coeff = polyfitw(xstart, xset[iy,*], xcheck, ndegree, xfit)  
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
         coeff = polyfitw(xstart, xset[iy,*], xmask[iy,*], $
          ndegree, xfit)
         xset[iy,ixbad] = xfit[ixbad]
      endif
   endfor

   ;----------
   ; Set FIBERMASK bit for any fiber with more than 20% of its positions
   ; masked, which includes any positions off the CCD

   ibad = where(total(1-xmask, 1) GT 0.20*ny)
   if (ibad[0] NE -1) then $
    fibermask[ibad] = fibermask[ibad] OR fibermask_bits('BADTRACE')

   return, xset
end
;------------------------------------------------------------------------------
