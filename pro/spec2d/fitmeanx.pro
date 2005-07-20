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
;    nord=, maxdev=, minsdev=, inmask=, mx=)
;
; INPUTS:
;   wset     - Initial wavelength solution for all fibers
;   lambda   - Air wavelengths corresponding to XPOS (in log-10 Angstroms)
;              [NLINE]
;   xpos     - Centroid positoins of sky lines in object image [NFIBER,NLINE]
;
; OPTIONAL KEYWORDS:
;   nord     - Order of fit to delta x positions; default to 4
;   maxdev   - Max abs difference (in pix) to reject outlying sky line
;              positions; default to 0.4
;   minsdev  - Minimum standard deviation for computing AVEINVVAR;
;              default to 0.01 pix.
;   inmask   - Input mask, set to 1 for good points, 0 for bad points
;              [NFIBER,NLINE]
;
; OUTPUTS:
;   xdiff    - Smooth fit to difference between measured sky positions
;              and those predicted from arc wavelength solution [NFIBER,NLINE]
;
; OPTIONAL OUTPUTS:
;   aveinvvar- Weights in xpos, set to zero for rejected positions [NLINE]
;   mx       - Sky line positions predicted from arc line solution
;              [NFIBER,NLINE]
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_iterstat
;
; REVISION HISTORY:
;   19-Oct-1999  Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
function fitmeanx, wset, lambda, xpos, aveinvvar, $
 nord=nord, maxdev=maxdev, minsdev=minsdev, inmask=inmask, mx=mx

   if (NOT keyword_set(nord)) then nord = 4
   if (NOT keyword_set(maxdev)) then maxdev = 0.4
   if (NOT keyword_set(minsdev)) then minsdev = 0.01

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
      if (!version.release LT '5.4') then begin
         junk = polyfitw(xaxis, rawdiff, qgood, nord, yfit)
      endif else begin
         junk = polyfitw(xaxis, rawdiff, qgood, nord, yfit, /double)
      endelse

      res1 = yfit - rawdiff

      ; Find which points are most deviant from the fit
      ; Set MASK=1 for bad points
      djs_iterstat, res1, sigrej=4.0, maxiter=3, mask=mask

      if (keyword_set(inmask)) then mask = mask AND inmask[*,i]

      ; Re-fit to RAWDIFF as a function of fiber number (rejecting outliers)
      igood = where(mask NE 0, ngood)
      if (ngood GT nord) then begin
         if (!version.release LT '5.4') then begin
            junk = polyfitw(xaxis, rawdiff, mask, nord, yfit)
         endif else begin
            junk = polyfitw(xaxis, rawdiff, mask, nord, yfit, /double)
         endelse

         xfit[*,i] = yfit
         sdev = stddev((yfit-rawdiff)[igood], /double)
      endif else begin
         ; Case with too few points to fit
         sdev = 0.0
      endelse

      if (sdev EQ 0.0) then aveinvvar[i] = 0 $
       else aveinvvar[i] = 1.0 / (sdev > minsdev)^2
      splog, 'In skyline number' ,i,' std-dev is ', sdev, ' pix'
   endfor

   return, xfit
end
;------------------------------------------------------------------------------
