;+
; NAME:
;   trace_crude
;
; PURPOSE:
;   Create a crude trace set given one position (eg, a center) in each trace.
;
; CALLING SEQUENCE:
;   xset = trace_crude( fimage, [xstart, ystart, radius=radius, yset=yset]  )
;
; INPUTS:
;   fimage     - Image
;   xstart     - Initial guesses for X centers (one for each trace).
;                If not set, then this code searches for all peaks at YSTART.
;   ystart     - Y positions corresponding to "xstart" (expected as integers).
;                There are three options for this parameter:
;                (1) One element of YSTART for each value of XSTART,
;                (2) A scalar value that is used for every XSTART, or
;                (3) Not set, in which case the center row is used.
;
; OPTIONAL KEYWORDS:
;   radius     - Radius for centroiding; default to 3.0
;
; OUTPUTS:
;   xset       - X centers for all traces
;
; OPTIONAL OUTPUTS:
;   yset       - Y centers for all traces
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_laxisgen()
;   Dynamic link to trace_crude.c
;
; REVISION HISTORY:
;   14-May-1999  Written by David Schlegel, Princeton.
;   12-Jul-1999  Added optional output YSET (DJS).
;-
;------------------------------------------------------------------------------
function trace_crude, fimage, xstart, ystart=ystart, radius=radius, yset=yset, $
	nave = nave, nmed=nmed

   ; Need 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - xset = trace_crude( fimage, [xstart, ystart=ystart, radius=radius] )'
      return, -1
   endif

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   if (NOT keyword_set(ystart)) then ystart = ny/2
   if (NOT keyword_set(radius)) then radius = 3.0
   if (NOT keyword_set(nmed)) then nmed = 1
   if (NOT keyword_set(nave)) then nave = 5

   ; Make a copy of the image
   ftemp = float(fimage)

   ; Median filter the entire image along columns by NMED rows
   if (nmed GT 1) then $   
   for ix=1, nx-1 do ftemp[ix,*] = median(transpose(ftemp[ix,*]), nmed)

   ; Boxcar-smooth the entire image along columns by NAVE rows
   for ix=1, nx-1 do ftemp[ix,*] = smooth(transpose(ftemp[ix,*]), nave)

   if (NOT keyword_set(xstart)) then begin
      ; Automatically find peaks for XSTART

      ; Extract NAVE rows from the image at YSTART
      nave = 1
      imrow = $
       ftemp[*,long(ystart[0]-0.5*(nave-1)):long(ystart[0]+0.5*(nave-1))]
      imrow = rebin(imrow, nx, 1)

      ; Boxcar smooth in both X and Y
;      imrow = smooth(imrow,3)

      ; Find all local peaks that are also above 1.25 times the median
      medval = median(imrow)
      xstart = where( imrow[1:nx-2] GT imrow[0:nx-3] $
                  AND imrow[1:nx-2] GT imrow[2:nx-1] $
                  AND imrow[1:nx-2] GT 1.25*medval) + 1
   endif

   if (N_elements(ystart) EQ 1) then $
    ypass = replicate(long(ystart), N_elements(xstart)) $
    else ypass = long(ystart)

   if (N_elements(xstart) NE N_elements(ypass)) then $
    message, 'Wrong number of elements for YSTART'

   ntrace = N_elements(xstart)
   xset = fltarr(ny, ntrace)
   result = call_external(getenv('IDL_EVIL')+'libspec2d.so', 'trace_crude', $
    nx, ny, ftemp, float(radius), ntrace, float(xstart), ypass, $
    xset)

   yset = djs_laxisgen([ny,nTrace], iaxis=0)

   return, xset
end
;------------------------------------------------------------------------------
