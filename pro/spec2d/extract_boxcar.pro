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
function extract_boxcar, fimage, xcen, ycen, radius=radius, idl=idl

   ; Need 2 parameters
   if (N_params() LT 2) then begin
      print, 'Syntax - fextract = extract_boxcar( fimage, xcen, [ycen, radius=radius] )'
      return, -1
   endif
   if (NOT keyword_set(radius)) then radius = 3.0

   if (NOT keyword_set(ycen)) then begin
	ycen=xcen
        ndim = (size(xcen))[0]
        if (ndim EQ 1) then $
          ycen = findgen(N_elements(xcen)) $
        else if (ndim EQ 2) then begin
          npix = (size(xcen))[1]
          nTrace = (size(xcen))[2]
	  for i=0,nTrace-1 do ycen[*,i] = findgen(npix)
        endif else message, 'xcen is not 1 or 2 dimensional'
   endif
         
   if (N_elements(xcen) NE N_elements(ycen)) then $
    message, 'Number of elements in XCEN and YCEN must be equal'

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   ncen = N_elements(xcen)

;   if (min(ycen) LT 0 OR max(ycen) GT ny-y) then $
;    message, 'YCEN contains values out of range'

   fextract = float(0 * xcen)

   if keyword_set(idl) then begin
      left = xcen - radius
      right = xcen + radius
      fextract = extract_asymbox2(fimage, left, right, ycen)
   endif else begin
     soname = filepath('libspec2d.'+idlutils_so_ext(), $
      root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='lib')
     result = call_external(soname, 'extract_boxcar', $
      nx, ny, float(fimage), float(radius), ncen, float(xcen), long(ycen), $
      fextract)
   endelse
   return, fextract
end
;------------------------------------------------------------------------------
