;+
; NAME:
;   fitsn
;
; PURPOSE:
;   Perform a simple parameter to fit to log S/N vs magnitude
;
; CALLING SEQUENCE:
;   coeffs = fitsn(mag, sn, /physical, sigrej=sigrej, maxiter=maxiter, $
;                  maxmag=maxmag, minmag=minmag, sigma=sigma)
;
; INPUTS:
;   mag        - S/N vector for fibers
;   sn         - Fiber magnitudes
;
; OPTIONAL KEYWORDS:
;   physical   - fit with a model of background and throughput
;   sigrej     - sigma rejection threshold
;   sigma      - standard deviation of residuals
;   maxiter    - Maximum number of iterations
;   minmag     - minimum magnitude limit to fit
;   maxmag     - maximum magnitude limit to fit
;
; OUTPUTS:
;   coeffs     - coefficients from fit
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   physical model with flux and sky background is not implemented
;
; PROCEDURES CALLED:
;   poly_fit
;
; REVISION HISTORY:
;   15-Apr-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function fitsn, mag, sn, physical=physical, sigrej=sigrej, maxiter=maxiter, $
    maxmag=maxmag, minmag=minmag, sigma=sigma

    if (NOT keyword_set(sigrej)) then sigrej = 5
    if (NOT keyword_set(maxiter)) then maxiter=5
    if (NOT keyword_set(maxmag)) then maxmag=22
    if (NOT keyword_set(minmag)) then minmag=14

    nspec = n_elements(sn)
    mask = (sn GT 0 AND mag GT minmag AND mag LT maxmag)

    good = where(mask)
    if (good[0] EQ -1) then return, -1

    logsn = sn*0.0-1.0
    logsn[good] = alog10(sn[good])

    for i=0, maxiter-1 do begin
      good = where(mask)
      if (good[0] EQ -1) then return, -1
      a = polyfitw(mag, logsn, mask, 1, yfit)
      
      diff = logsn - yfit
      djs_iterstat, diff[good], sigrej=sigrej, sigma=sigma, mask=smask
      treject = total(smask)
      if (treject EQ 0) then return,a $
      else begin
         mask[good] = mask[good] * (1-smask)
         if (total(mask) LT 3) then return,-1
      endelse
    endfor

    return,-1

end


    
