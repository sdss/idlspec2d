;+
; NAME:
;   trace_fweight
;
; PURPOSE:
;   Recenter a trace using flux-weighting centers.
;
; CALLING SEQUENCE:
;   xnew = trace_fweight( fimage, xcen, ycen, [radius=radius, xerr=xerr] )
;
; INPUTS:
;   fimage     - Image
;   xcen       - Initial guesses for X centers
;   ycen       - Y positions corresponding to "xcen" (expected as integers)
;
; OPTIONAL KEYWORDS:
;   radius     - Radius for centroiding; default to 3.0
;
; OUTPUTS:
;   xnew       - New X centers
;
; OPTIONAL OUTPUTS:
;   xerr       - Formal errors for XNEW
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   Dynamic link to trace_fweight.c
;
; REVISION HISTORY:
;   24-Mar-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function trace_fweight, fimage, xcen, ycen, radius=radius, xerr=xerr

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - xnew = trace_fweight( fimage, xcen, ycen, [radius=radius, $'
      print, ' xerr=xerr] )'
      return, -1
   endif
   if (NOT keyword_set(radius)) then radius = 3.0

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   ncen = N_elements(xcen)
   xnew = float(xcen)
   xerr = fltarr(ncen)
   invtemp = 1.0 / (abs(fimage) + 20.0)

   result = call_external(getenv('IDL_EVIL')+'libspec2d.so', 'trace_fweight', $
    nx, ny, float(fimage), float(invtemp), $
         float(radius), ncen, xnew, long(ycen), xerr)

   return, xnew
end
;------------------------------------------------------------------------------
