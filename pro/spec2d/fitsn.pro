;+
; NAME:
;   fitsn
;
; PURPOSE:
;   Perform a simple parameter fit to log S/N vs magnitude
;
; CALLING SEQUENCE:
;   coeffs = fitsn(mag, snvec, [ sigrej=, maxiter=, $
;    fitmag=, colorband=, sigma= ] )
;
; INPUTS:
;   mag        - Fiber magnitudes
;   snvec      - S/N vector for fibers
;
; OPTIONAL KEYWORDS:
;   sigrej     - Sigma rejection threshold; default to 3
;   maxiter    - Maximum number of rejection iterations; default to 5
;   fitmag     - Magnitude range over which to fit (S/N) as function of mag;
;                no defaults unless COLOR is set.
;   colorband  - If set, then override FITMAG with the following fitting ranges:
;                  'B' : [18.20, 19.70]  <-- Used for Son-of-Spectro
;                  'R' : [17.90, 19.40]  <-- Used for Son-of-Spectro
;                  'G' : [18.20, 19.70]
;                  'R' : [18.25, 19.75]
;                  'I' : [17.90, 19.40]
;                The above ranges will be extended to [0,23] if fewer than
;                20 points are found in the above ranges, but at least 3
;                good points at any magnitude.
;
; OUTPUTS:
;   coeffs     - Coefficients from fit; return 0 if fit failed
;
; OPTIONAL OUTPUTS:
;   sigma      - Standard deviation of residuals
;   fitmag     - (Modified if COLOR is set)
;
; COMMENTS:
;   If there are fewer than 3 points, then return COEFFS=0.
;
;   This function is called by the following routes in Son-of-Spectro:
;     APO_PLOTSN -> PLOTSN -> FITSN
;     QUICKEXTRACT -> FITSN
;
;   This function is called by the following routes in Spectro-2D:
;     PLATESN -> PLOTSN -> FITSN
;     EXTRACT_OBJECT -> FITSN
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
 fitmag=fitmag, colorband=colorband, sigma=sigma

   if (NOT keyword_set(sigrej)) then sigrej = 3
   if (NOT keyword_set(maxiter)) then maxiter = 5

   if (keyword_set(colorband)) then begin
      case strupcase(colorband) of
         'B' : fitmag = [18.20, 19.70]
         'R' : fitmag = [17.90, 19.40]
         'G' : fitmag = [18.20, 19.70]
         'R' : fitmag = [18.25, 19.75]
         'I' : fitmag = [17.90, 19.40]
         else: message, 'Invalid COLOR keyword value'
      endcase
      if (total(snvec GT 0 AND mag GT fitmag[0] AND mag LT fitmag[1]) $
       LT 20 AND total(snvec GT 0) GT 2) then fitmag = [0,23]
   endif

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
