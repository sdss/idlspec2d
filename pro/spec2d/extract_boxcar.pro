;+
; NAME:
;   extract_boxcar
;
; PURPOSE:
;   Extract the total flux within a boxcar window at many positions.
;
; CALLING SEQUENCE:
;   fextract = extract_boxcar( fimage, xcen, ycen, [radius=radius] )
;
; INPUTS:
;   fimage     - Image
;   xcen       - Initial guesses for X centers
;   ycen       - Y positions corresponding to "xcen" (expected as integers)
;
; OPTIONAL KEYWORDS:
;   radius     - Radius of extraction; default to 3.0
;
; OUTPUTS:
;   fextract   - Extracted flux at each position specified by (xcen, ycen)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   Dynamic link to extract_boxcar.c
;
; REVISION HISTORY:
;   24-Mar-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function extract_boxcar, fimage, xcen, ycen, radius=radius

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - fextract = extract_boxcar( fimage, xcen, ycen, [radius=radius] )'
      return, -1
   endif
   if (NOT keyword_set(radius)) then radius = 3.0
   if (N_elements(xcen) NE N_elements(ycen)) then $
    message, 'Number of elements in XCEN and YCEN must be equal'

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   ncen = N_elements(xcen)

;   if (min(ycen) LT 0 OR max(ycen) GT ny-y) then $
;    message, 'YCEN contains values out of range'

   fextract = float(0 * xcen)
   result = call_external(getenv('IDL_EVIL')+'libspec2d.so', 'extract_boxcar', $
    nx, ny, float(fimage), float(radius), ncen, float(xcen), long(ycen), $
    fextract)

   return, fextract
end
;------------------------------------------------------------------------------
