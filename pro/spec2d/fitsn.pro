;+
; NAME:
;   fitsn
;
; PURPOSE:
;   Perform a simple parameter fit to log S/N vs magnitude
;
; CALLING SEQUENCE:
;   coeffs = fitsn(mag, snvec, [ sigrej=, maxiter=, $
;    fitmag=, sigma= ] )
;
; INPUTS:
;   mag        - Fiber magnitudes
;   snvec      - S/N vector for fibers
;
; OPTIONAL KEYWORDS:
;   sigrej     - Sigma rejection threshold; default to 3
;   maxiter    - Maximum number of rejection iterations; default to 5
;   fitmag     - Magnitude range over which to fit (S/N) as function of mag;
;                default to [14,22]
;
; OUTPUTS:
;   coeffs     - Coefficients from fit; return 0 if fit failed
;
; OPTIONAL OUTPUTS:
;   sigma      - Standard deviation of residuals
;
; COMMENTS:
;   If there are fewer than 3 points, then return COEFFS=0.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_iterstat
;
; REVISION HISTORY:
;   15-Apr-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function fitsn, mag, snvec, sigrej=sigrej, maxiter=maxiter, $
 fitmag=fitmag, sigma=sigma

    if (NOT keyword_set(sigrej)) then sigrej = 3
    if (NOT keyword_set(maxiter)) then maxiter = 5
    if (NOT keyword_set(fitmag)) then fitmag = [14, 22]

    sigma = 0
    nspec = n_elements(snvec)
    mask = (snvec GT 0 AND mag GT fitmag[0] AND mag LT fitmag[1])

    igood = where(mask, ngood)
    if (ngood LE 2) then return, 0

    logsn = snvec*0.0 - 1.0 ; Arbitrarily set bad values to -1, though these
                         ; values are masked from the fit anyway
    logsn[igood] = alog10(snvec[igood])

    for i=0, maxiter-1 do begin
       igood = where(mask, ngood)
       if (ngood LE 2) then return, 0
       if (!version.release LT '5.4') then begin
          coeffs = polyfitw(mag, logsn, mask, 1, yfit)
       endif else begin
          coeffs = polyfitw(mag, logsn, mask, 1, yfit, /double)
       endelse
      
       diff = logsn - yfit
       djs_iterstat, diff[igood], sigrej=sigrej, sigma=sigma, mask=smask
       treject = total(1-smask)
       if (treject EQ 0) then return, coeffs

       mask[igood] = mask[igood] * smask
       if (total(mask) LE 2) then return, 0
    endfor

    return, coeffs
end
;------------------------------------------------------------------------------
