;+
; NAME:
;   fitmeanx
;
; PURPOSE:
;   Use fitans to tweak trace and sigma
;   This just perturbs the input values of xcen and sigma
;
; CALLING SEQUENCE:
;   xdiff = fitmeanx(wset, lambda, xpos, xinvvar, $
;                      nord=nord, maxdev=maxdev, mx=mx)
;
; INPUTS:
;   wset     - Initial wavelength solution
;   lambda   - air wavelengths corresponding to xpos
;   xpos     - centroids of sky lines in object image
;
; OPTIONAL KEYWORDS:
;   nord     - order of fit to delta x positions
;   maxdev   - max abs difference to reject outlying sky line positions
;
; OUTPUTS:
;   xdiff    - smooth fit to difference between measured sky positions
;                and those predicted from arc wavelength solution
;
; OPTIONAL OUTPUTS:
;   xinvvar  - weights in xpos (mask)
;   mx       - sky line positions predicted from arc line solution
;
; COMMENTS:
;    LAMBDA = log10-wavelength
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   19-Oct-1999  Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
function fitmeanx, wset, lambda, xpos, xinvvar, nord=nord, maxdev=maxdev, mx=mx

   if (NOT keyword_set(nord)) then nord = 4
   if (NOT keyword_set(maxdev)) then maxdev = 0.4
   dims = size(xpos, /dim)
   nfiber = dims[0]
   nlambda = dims[1]

   ; Evaluate the trace set to get the fit pixel position at each arc lambda
   mx = transpose(traceset2pix(wset, lambda))

   x = findgen(nfiber) / float(nfiber) ; From [0,1)
   xnew = xpos
   xtemp = xpos
   xfit = xpos
   xinvvar = xpos
   xshift = xpos

   ; Loop through each arc line...

   for i=0, nlambda-1 do begin
      dif = xpos[*,i] - mx[*,i] ; Measured position minus MX
      xtemp[*,i] = dif

      djs_iterstat, dif, median = mm
;
;	Throw out pixels with that are really bad
;
      qgood = abs(dif-mm) LT maxdev

      ; Fit to DIF as a function of fiber number
      junk = polyfitw(x, dif, qgood, nord, yfit)
      res1 = yfit - dif

      ; Find which points are most deviant (4-sigma cut; 2 iterations)
      qgood = abs(res1) LT 4.*stddev(res1)
      qgood = abs(res1) LT 4.*stddev(res1*qgood)
      xfit[*,i] = yfit

      ; Re-fit to DIF as a function of fiber number (rejecting outliers)
      junk = polyfitw(x, dif, qgood, nord, yfit)

      xinvvar[*,i] = 0.0
      thisdev = stddev((yfit-dif)*qgood)
      if (thisdev GT 0.0) then xinvvar[*,i]  = 1.0/thisdev^2
      splog, 'In trace' ,i,' Stan. dev is ', thisdev*1000.0, ' mPix'

      ; Return the position where the wavelength should fall (MX) plus
      ; a fitted-function to the deviations from those positions.
      xnew[*,i] = mx[*,i] + yfit

   endfor
   return, xfit
end
;------------------------------------------------------------------------------
