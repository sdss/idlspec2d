;+
; NAME:
;   trace_crude
;
; PURPOSE:
;   Create a crude trace set given one position (eg, a center) in each trace.
;
; CALLING SEQUENCE:
;   xset = trace_crude( fimage, invvar, [xstart=, ystart=, radius=, yset=, $
;    nave=, nmed=, thresh=, maxerr=, maxshifte=, maxshift0=, xerr= ] )
;
; INPUTS:
;   fimage     - Image
;
; OPTIONAL INPUTS:
;   invvar     - Inverse variance (weight) image
;   xstart     - Initial guesses for X centers (one for each trace).
;                If not set, then this code searches for all peaks at YSTART.
;   ystart     - Y positions corresponding to "xstart" (expected as integers).
;                There are three options for this parameter:
;                (1) One element of YSTART for each value of XSTART,
;                (2) A scalar value that is used for every XSTART, or
;                (3) Not set, in which case the center row is used.
;   radius     - Radius for centroiding; default to 3.0
;   nmed       - Median filtering size down columns before performing trace;
;                default to 1
;   nave       - Averaging size down columns before performing trace;
;                default to 5
;   thresh     - Threshold for initial peak finding; if not set, then use
;                1.0 times the median of the row(s) used for the initial peaks.
;   maxerr     - Maximum error in centroid allowed for valid recentering;
;                default to 0.2
;   maxshifte  - Maximum shift in centroid allowed for valid recentering;
;                default to 0.1
;   maxshift0  - Maximum shift in centroid allowed for initial row;
;                default to 0.5
;
; OUTPUTS:
;   xset       - X centers for all traces
;
; OPTIONAL OUTPUTS:
;   yset       - Y centers for all traces
;   xerr       - Errors for XSET
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_laxisgen()
;
;   Dynamic link to trace_crude.c
;
; REVISION HISTORY:
;   14-May-1999  Written by David Schlegel, Princeton.
;   12-Jul-1999  Added optional output YSET (DJS).
;   06-Aug-1999  Added optional outpust XERR (DJS).
;-
;------------------------------------------------------------------------------
function trace_crude, fimage, invvar, xstart=xstart, ystart=ystart, $
 radius=radius, yset=yset, nave=nave, nmed=nmed, thresh=thresh, $
 maxerr=maxerr, maxshifte=maxshift, maxshift0=maxshift0, xerr=xerr

   ; Need 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - xset = trace_crude( fimage, [ invvar, xstart=, ystart=, $'
      print, ' radius=, nave=, nmed=, maxerr=, maxshifte=, maxshift0=, xerr= ] )'
      return, -1
   endif

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   if (NOT keyword_set(ystart)) then ystart = ny/2
   if (NOT keyword_set(radius)) then radius = 3.0
   if (NOT keyword_set(nmed)) then nmed = 1
   if (NOT keyword_set(nave)) then nave = 5
   if (NOT keyword_set(maxerr)) then maxerr = 0.2
   if (NOT keyword_set(maxshift)) then maxshift = 0.1
   if (NOT keyword_set(maxshift0)) then maxshift0 = 0.5

   ; Make a copy of the image and error map
   imgtemp = float(fimage)
   if (keyword_set(invvar)) then begin
      invtemp = float(invvar)
   endif else begin
      invtemp = 1.0 / (fimage > 1)
   endelse

   ; Median filter the entire image along columns by NMED rows
   if (nmed GT 1) then $   
    for ix=1, nx-1 do imgtemp[ix,*] = median(transpose(imgtemp[ix,*]), nmed)

   ; Boxcar-sum the entire image along columns by NAVE rows
   if (nave GT 1) then begin
      kernal = transpose(intarr(nave) + 1.0)
      imgconv = convol(imgtemp*invtemp, kernal, /edge_truncate)

      ; Add the variances
      invtemp = convol(invtemp, kernal, /edge_truncate)

      ; Look for pixels with infinite errors - replace with original values
      ibad = where(invtemp EQ 0, nbad)
      if (nbad GT 0) then begin
         invtemp[ibad] = 1
         imgconv[ibad] = imgtemp[ibad]
      endif

      ; Renormalize the summed image by the weights
      imgtemp = imgconv / invtemp
   endif

   if (NOT keyword_set(xstart)) then begin
      ; Automatically find peaks for XSTART

      ; Extract NSUM rows from the image at YSTART
      nsum = 1
      imrow = $
       imgtemp[*,long(ystart[0]-0.5*(nsum-1)):long(ystart[0]+0.5*(nsum-1))]
      imrow = rebin(imrow, nx, 1)

      if (keyword_set(thresh)) then mthresh = thresh $
       else mthresh = median(imrow)

      ; Boxcar smooth along X
      imrow = smooth(imrow,3)

      ; Find all local peaks that are also above THESH times the median
;      xstart = where( imrow[1:nx-2] GT imrow[0:nx-3] $
;                  AND imrow[1:nx-2] GT imrow[2:nx-1] $
;                  AND imrow[1:nx-2] GT 1.0*medval) + 1
      rderiv = imrow[1:nx-1] - imrow[0:nx-2]
      izero = where( rderiv[0:nx-3] GT 0 AND rderiv[1:nx-2] LE 0 $
       AND imrow[1:nx-2] GT mthresh)

      if izero[0] EQ -1 then message, 'No peaks found'

      xstart = izero + 0.5 + rderiv[izero] / (rderiv[izero] - rderiv[izero+1])
   endif

   if (N_elements(ystart) EQ 1) then $
    ypass = replicate(long(ystart), N_elements(xstart)) $
    else ypass = long(ystart)

   if (N_elements(xstart) NE N_elements(ypass)) then $
    message, 'Wrong number of elements for YSTART'

   ntrace = N_elements(xstart)
   xset = fltarr(ny, ntrace)
   xerr = fltarr(ny, ntrace)

   soname = filepath('libspec2d.so', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='lib')
   result = call_external(soname, 'trace_crude', $
    nx, ny, imgtemp, invtemp, float(radius), ntrace, float(xstart), ypass, $
    xset, xerr, float(maxerr), float(maxshift), float(maxshift0))

   yset = djs_laxisgen([ny,nTrace], iaxis=0)

   return, xset
end
;------------------------------------------------------------------------------
