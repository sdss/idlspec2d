;+
; NAME:
;   trace_fix
;
; PURPOSE:
;   Fix a set of trace centers by replacing traces that converge.
;
; CALLING SEQUENCE:
;   xnew = trace_fix( xcen, [minsep= , ngrow=, ycen=, xerr=])
;
; INPUTS:
;   xcen       - X centers for all traces [ny,nTrace]
;
; OPTIONAL INPUTS:
;   minsep     - Minimum separation between adjacent traces.  Smaller
;                separations are regarded as bad traces.  Default to 5.5.
;   ngrow      - Replace all pixels within MINSEP of its adjacent trace,
;                plus NGROW of its neighboring pixels.  Default to 20.
;   ycen       - Y centers corresponding to XCEN.
;   xerr       - X errors corresponding to XCEN.
;
; OUTPUTS:
;   xnew       - Modified XCEN; traces may be removed or shifted.
;   ycen       - Modified YCEN; columns may be removed.
;   xerr       - Modified XERR; columns may be removed.
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; INTERNAL PROCEDURES:
;   remove_column
;
; REVISION HISTORY:
;   13-Aug-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
; Remove column number COLNUM from an array, ARR.

pro remove_column, arr, colnum

   if (keyword_set(arr)) then begin
      ncol = (size(arr,/dim))[1]

      if (colnum EQ 0) then arr = arr[*,1:ncol-1] $
       else if (colnum EQ ncol-1) then arr = arr[*,0:ncol-2] $
       else arr = [ [arr[*,0:colnum-1]], [arr[*,colnum+1:ncol-1]] ]
   endif

   return
end
;------------------------------------------------------------------------------
function trace_fix, xcen, minsep=minsep, ngrow=ngrow, ycen=ycen, xerr=xerr

   if (NOT keyword_set(minsep)) then minsep = 5.5
   if (NOT keyword_set(ngrow)) then ngrow = 20

   xnew = xcen
   ny = (size(xnew,/dimens))[0]
   ntrace = (size(xnew,/dim))[1]

   ; Decide where neighboring traces are too close to one another.
   ; Do this in every row by looking at distances between neighboring centers.
   xdiff = abs( xnew[*,1:ntrace-1] - xnew[*,0:ntrace-2] )
   xbad = xdiff LT minsep

   ; First, look for any traces that are ALWAYS too close to their neighbor.
   itrace = 0
   while (itrace LT ntrace-1) do begin
      if (total(xbad[*,itrace]) EQ ny) then begin
         ; Either trace number ITRACE or ITRACE+1 should be discarded.
         ; Get rid of the one with the most curvature.
         disp1 = stddev(xnew[0:ny-2,itrace] - xnew[1:ny-1,itrace], /double)
         disp2 = stddev(xnew[0:ny-2,itrace+1] - xnew[1:ny-1,itrace+1], /double)
         if (disp1 GT disp2) then badnum = itrace $
          else badnum = itrace + 1
         remove_column, xnew, badnum
         remove_column, ycen, badnum
         remove_column, xerr, badnum

         ; With one trace removed, recompute NTRACE, XDIFF, and XBAD
         ntrace = ntrace-1
         xdiff = abs( xnew[*,1:ntrace-1] - xnew[*,0:ntrace-2] )
         xbad = xdiff LT minsep
      endif else begin
         itrace = itrace+1
      endelse
   endwhile

   ; Now look for traces that are sometimes too close to their neighbors,
   ; and shift them so that they are not any more.
   for itrace=0, ntrace-2 do begin
      ibad = where(xbad[*,itrace], nbad)
      if (nbad GT 0) then begin
         igood = where(NOT xbad[*,itrace], ngood)

         ; Identify which trace number has gone bad
         if (itrace EQ 0) then begin
            space1 = median( xnew[igood,itrace+1] - xnew[igood,itrace] )
            space2 = median( xnew[ibad,itrace+1] - xnew[ibad,itrace] )
            if (space2 GT space1+1.5) then ifix = itrace+1 $
             else ifix = itrace
         endif else begin
            space1 = median( xnew[igood,itrace] - xnew[igood,itrace-1] )
            space2 = median( xnew[ibad,itrace] - xnew[ibad,itrace-1] )
            if (space2 GT space1+1.5) then ifix = itrace $
             else ifix = itrace+1
         endelse

         ; Fix trace number IFIX
         if (ifix GT 1) then begin
            ; Grow the bad pixels to their neighbors
            ibad = where( smooth(float(xbad[*,itrace]), 1+2*ngrow) GT 0)
            xadd = median( xnew[igood,ifix] - xnew[igood,ifix-1] )
            xnew[ibad,ifix] = xnew[ibad,ifix-1] + xadd
         endif

      endif
   endfor

   return, xnew
end
;------------------------------------------------------------------------------
