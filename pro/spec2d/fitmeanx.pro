;+
; NAME:
;   fitmeanx
;
; PURPOSE:
;   Return the position where the wavelength should fall (MX) plus
;   a fitted-function to the deviations from those positions.
;
; CALLING SEQUENCE:
;   xdiff = fitmeanx(wset, lambda, xpos, aveinvvar, $
;    nord=, maxdev=, mx=mx)
;
; INPUTS:
;   wset     - Initial wavelength solution
;   lambda   - Air wavelengths corresponding to XPOS (in log-10 Angstroms)
;   xpos     - Centroid positoins of sky lines in object image
;
; OPTIONAL KEYWORDS:
;   nord     - Order of fit to delta x positions; default to 4
;   maxdev   - Max abs difference (in pix) to reject outlying sky line
;              positions; default to 0.4
;
; OUTPUTS:
;   xdiff    - Smooth fit to difference between measured sky positions
;              and those predicted from arc wavelength solution
;
; OPTIONAL OUTPUTS:
;   aveinvvar- Weights in xpos, set to zero for rejected positions
;   mx       - Sky line positions predicted from arc line solution
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   19-Oct-1999  Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
function fitmeanx, wset, lambda, xpos, aveinvvar, $
 nord=nord, maxdev=maxdev, mx=mx

   if (NOT keyword_set(nord)) then nord = 4
   if (NOT keyword_set(maxdev)) then maxdev = 0.4

   dims = size(xpos, /dim)
   nfiber = dims[0]
   nlambda = dims[1]

   ; Evaluate the trace set to get the fit pixel position at each arc lambda
   mx = transpose(traceset2pix(wset, lambda))

   xaxis = findgen(nfiber) / float(nfiber) ; From [0,1)
   xfit = xpos
   aveinvvar = fltarr(nlambda)
   xshift = xpos

   ; Loop through each arc line...

   for i=0, nlambda-1 do begin
      rawdiff = xpos[*,i] - mx[*,i] ; Measured position minus predicted

      djs_iterstat, rawdiff, median=mm

      ; Reject pixels from initial fit that are very deviant from the median
      qgood = abs(rawdiff-mm) LT maxdev

      ; Fit to RAWDIFF as a function of fiber number
      junk = polyfitw(xaxis, rawdiff, qgood, nord, yfit)
      res1 = yfit - rawdiff

      ; Find which points are most deviant from the fit
      djs_iterstat, res1, sigrej=4.0, maxiter=3, mask=mask

      ; Re-fit to RAWDIFF as a function of fiber number (rejecting outliers)
      junk = polyfitw(xaxis, rawdiff, 1-mask, nord, yfit)
      xfit[*,i] = yfit

      igood = where(mask EQ 1, ngood)
      if (ngood LE 1) then  sdev = 0.0 $
       else sdev = stddev((yfit-rawdiff)[igood])
      if (sdev EQ 0.0) then aveinvvar[i] = 0 $
       else  aveinvvar[i]  = 1.0 / sdev^2

      splog, 'In skyline number' ,i,' std-dev is ', sdev, ' pix'

   endfor

   return, xfit
end
;------------------------------------------------------------------------------
