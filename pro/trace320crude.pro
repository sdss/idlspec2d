;+
; NAME:
;   trace320crude
;
; PURPOSE:
;   Calling script to return 320 full traces using TRACE_CRUDE.
;
; CALLING SEQUENCE:
;   xset = trace320crude( fimage, invvar, [mthresh=, ystart=, nmed=, xgood=, $
;    yset=, maxerr=, maxshifte=, maxshift0=, xerr=, minsep=, ngrow=  ] )
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
; OPTIONAL INPUTS FOR TRACE_FIX:
;   minsep     - Minimum separation between adjacent traces.  Smaller
;                separations are regarded as bad traces.  Default to 5.5.
;   ngrow      - Replace all pixels within MINSEP of its adjacent trace,
;                plus NGROW of its neighboring pixels.  Default to 20.
;
; OUTPUTS:
;   xset       - X centers for all traces
;
; OPTIONAL OUTPUTS:
;   yset       - Y centers for all traces
;   xerr       - Errors for XSET
;   xgood        - Set to 1 for fibers that were actually found, 0 otherwise
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   trace_crude
;   trace320cen
;
; REVISION HISTORY:
;   13-Sep-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function trace320crude, fimage, invvar, ystart=ystart, nmed=nmed, xgood=xgood, $
 radius=radius, yset=yset, maxerr=maxerr, maxshifte=maxshift, $
 maxshift0=maxshift0, xerr=xerr, minsep=minsep, ngrow=ngrow

   if (NOT keyword_set(minsep)) then minsep = 5.5
   if (NOT keyword_set(ngrow)) then ngrow = 20

   ; Find the 320 X-centers in the row specified by YSTART
   xstart = trace320cen(fimage, mthresh=mthresh, ystart=ystart, nmed=nmed, $
    xgood=xgood)
   ntrace = N_elements(xstart) ; Better be 320

   ; Trace those 320
   xset = trace_crude( fimage, invvar, xstart=xstart, ystart=ystart, $
    radius=radius, yset=yset, maxerr=maxerr, maxshifte=maxshifte, $
    maxshift0=maxshift0, xerr=xerr )

   ; Construct an image of bad X positions.
   ; Start by setting all bad initial X centers.
   xmask = byte( 0*xset )
   i = where(xgood EQ 0)
   if (i[0] NE -1) then xmask(i,*) = 1

   xnew = xset
   ny = (size(xnew,/dimens))[0]

   ; Decide where neighboring traces are too close to one another.
   ; Do this in every row by looking at distances between neighboring centers.
   xdiff = abs( xnew[*,1:ntrace-1] - xnew[*,0:ntrace-2] )
   xbad = xdiff LT minsep


; ???

   ; Now look for traces that are sometimes too close to their neighbors,
   ; and shift them so that they are not any more.
   for itrace=0, ntrace-2 do begin
      ibad = where(xbad[*,itrace], nbad)
      if (nbad GT 0) then begin
         igood = where(xbad[*,itrace] EQ 0, ngood)

         ; Identify which trace number has gone bad
         if (xgood[itrace] EQ 1 AND xgood[itrace+1] EQ 0) then begin
            ifix = itrace+1
            ihelp = itrace
         endif else if (xgood[itrace] EQ 0 AND xgood[itrace+1] EQ 1) then begin
            ifix = itrace
            ihelp = itrace+1
         endif else begin
            if (itrace MOD 20 EQ 0) then begin
               space1 = median( xnew[igood,itrace+1] - xnew[igood,itrace] )
               space2 = median( xnew[ibad,itrace+1] - xnew[ibad,itrace] )
               if (space2 GT space1+1.5) then begin
                  ifix = itrace 
                  ihelp = itrace + 1
               endif else begin
                  ifix = itrace + 1
                  ihelp = itrace
               endelse
            endif else begin
               space1 = median( xnew[igood,itrace] - xnew[igood,itrace-1] )
               space2 = median( xnew[ibad,itrace] - xnew[ibad,itrace-1] )
               if (space2 GT space1+1.5) then begin
                  ifix = itrace - 1 
                  ihelp = itrace 
               endif else begin
                  ifix = itrace 
                  ihelp = itrace - 1
               endelse
            endelse
         endelse

print,'fix ', ifix
         ; Fix trace number IFIX
         if (ifix GE 0 AND ihelp GE 0) then begin
            ; Grow the bad pixels to their neighbors
            ibad = where( smooth(float(xbad[*,itrace]), 1+2*ngrow) GT 0)
            xadd = median( xnew[igood,ifix] - xnew[igood,ihelp] )
            xnew[ibad,ifix] = xnew[ibad,ihelp] + xadd 
         endif

      endif
   endfor

   return, xnew
end
;------------------------------------------------------------------------------
