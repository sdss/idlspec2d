;+
; NAME:
;   trace_gweight
;
; PURPOSE:
;   Recenter a trace using gaussian-weighted centers.
;
; CALLING SEQUENCE:
;   xnew = trace_fweight( fimage, xcen, ycen, [sigma=sigma, xerr=xerr, 
;    invvar=invvar] )
;
; INPUTS:
;   fimage     - Image
;   xcen       - Initial guesses for X centers
;   ycen       - Y positions corresponding to "xcen" (expected as integers)
;
; OPTIONAL KEYWORDS:
;   radius     - Sigma in pixels; default to 1.0
;   invvar     - Inverse variance of image used only in computing errors XERR.
;                If not set, then INVVAR=1 is used.
;
; OUTPUTS:
;   xnew       - New X centers
;
; OPTIONAL OUTPUTS:
;   xerr       - Formal errors for XNEW; set equal to 999.0 if there are any
;                masked pixels in a centroiding region (e.g., if INVVAR=0)
;
; COMMENTS:
;
; EXAMPLES:
;
; REVISION HISTORY:
;   17-Jan-2000  Written by Scott Burles, Chicago
;-
;------------------------------------------------------------------------------
function trace_gweight, fimage, xcen, ycen, sigma=sigma, xerr=xerr, $
 invvar=invvar

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - xnew = trace_gweight( fimage, xcen, ycen, $'
      print, '       [sigma=sigma, xerr=xerr, invvar=invvar] )'
      return, -1
   endif
   if (NOT keyword_set(sigma)) then sigma = 1.0

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   ncen = N_elements(xcen)
   xnew = float(xcen)
   xerr = 0.0 * xnew ; Allocate memory

   if (NOT keyword_set(invvar)) then invvar = 0.0 * fimage + 1.0

;   result = call_external(getenv('IDLSPEC2D_DIR')+'/lib/libspec2d.so', $
;    'trace_gweight', $
;    nx, ny, float(fimage), float(invvar), $
;    float(sigma), ncen, xnew, long(ycen), xerr)

   soname = filepath('libspec2d.so', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='lib')
   result = call_external(soname, 'trace_gweight', $
    nx, ny, float(fimage), float(invvar), $
    float(sigma), ncen, xnew, long(ycen), xerr)

   return, xnew
end
;------------------------------------------------------------------------------
