;+
; NAME:
;   trace_fix
;
; PURPOSE:
;   Fix a set of trace centers by replacing traces that converge.
;
; CALLING SEQUENCE:
;   xnew = trace_fix( xcen, [minsep= ])
;
; INPUTS:
;   xcen       - X centers for all traces [ny,nTrace]
;
; OPTIONAL INPUTS:
;   minsep     - Minimum separation between adjacent traces.  Smaller
;                separations are regarded as bad traces.  Default to 1.5.
;   ngrow      - Replace all pixels within MINSEP of its adjacent trace,
;                plus NGROW of its neighboring pixels.  Default to 20.
;
; OUTPUTS:
;   xnew       - Modified XCEN [ny,nTrace]
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   13-Aug-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function trace_fix, xcen, minsep=minsep, ngrow=ngrow

   if (NOT keyword_set(minsep)) then minsep = 1.5
   if (NOT keyword_set(ngrow)) then ngrow = 20

   xnew = xcen

   ; Decide where neighboring traces are too close to one another
   ntrace = (size(xcen,/dim))[1]
   xdiff = abs( xcen[*,1:ntrace-1] - xcen[*,0:ntrace-2] )
   xbad = xdiff LT minsep

   for itrace=0, ntrace-2 do begin
      ibad = where(xbad[*,itrace], nbad)
      if (nbad GT 0) then begin
         igood = where(NOT xbad[*,itrace], ngood)

         ; Identify which trace number has gone bad
         if (itrace EQ 0) then begin
            space1 = median( xcen[igood,itrace+1] - xcen[igood,itrace] )
            space2 = median( xcen[ibad,itrace+1] - xcen[ibad,itrace] )
            if (space2 GT space1+1.5) then ifix = itrace+1 $
             else ifix = itrace
         endif else begin
            space1 = median( xcen[igood,itrace] - xcen[igood,itrace-1] )
            space2 = median( xcen[ibad,itrace] - xcen[ibad,itrace-1] )
            if (space2 GT space1+1.5) then ifix = itrace $
             else ifix = itrace+1
         endelse

         ; Fix trace number IFIX
         if (ifix GT 1) then begin
            ; Grow the bad pixels to their neighbors
            ibad = where( smooth(float(xbad[*,itrace]), 1+2*ngrow) GT 0)
            xadd = median( xcen[igood,ifix] - xcen[igood,ifix-1] )
            xnew[ibad,ifix] = xnew[ibad,ifix-1] + xadd
         endif

      endif
   endfor

   return, xnew
end
