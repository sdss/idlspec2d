;+
; NAME:
;   fitsn
;
; PURPOSE:
;   Perform a simple parameter fit to log S/N vs magnitude
;
; CALLING SEQUENCE:
;   coeffs = fitsn(mag, snvec, [ sigrej=, maxiter=, sncode=, $
;    filter=, sigma=, specsnlimit=, sn2= ] )
;
; INPUTS:
;   mag        - Fiber magnitudes
;   snvec      - S/N vector for fibers
;   sncode     - Pipeline code for determining fit range
;   filter     - Filter band for determining fit range
;
; OPTIONAL KEYWORDS:
;   sigrej     - Sigma rejection threshold; default to 3
;                (only used in computing SIGMA)
;   maxiter    - Maximum number of rejection iterations; default to 5
;                (no longer used with median-fit line)
;   plugmap    - Plugmap structure (required for sn2ext output)
;
; OUTPUTS:
;   coeffs     - Coefficients from fit; return 0 if fit failed;
;                log10(S/N) = coeffs[0] - coeffs[1] * mag
;
; OPTIONAL OUTPUTS:
;   sigma      - Standard deviation of residuals
;   specsnlimit- Structure with FITMAG,SNMAG,FIDUCIAL_COEFF
;   sn2        - Fit value of (S/N)^2 at fiducial magnitude
;   dered_sn2  - Extinction corrected (S/N)^2 at fiducial magnitude
;                (requires plugmap as input)
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
function fitsn, mag, snvec, sigrej=sigrej, maxiter=maxiter, plugmap=plugmap, $
 sncode=sncode, filter=filter, sigma=sigma, specsnlimit=specsnlimit1, $
 sn2=sn2, dered_sn2=dered_sn2

   common com_fitsn, specsnlimit

   if (NOT keyword_set(sigrej)) then sigrej = 3
   if (NOT keyword_set(maxiter)) then maxiter = 5

   if (NOT keyword_set(specsnlimit)) then begin
      specsnlimit = yanny_readone(djs_filepath('opSNlimits.par', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdir='opfiles'), 'SPECSNLIMIT')
      if (NOT keyword_set(specsnlimit)) then $
       message, 'opSNlimits.par file not found!'
   endif

   i = where(specsnlimit.sncode EQ sncode $
    AND specsnlimit.filter EQ filter, ct)
   if (ct EQ 0) then $
    message, 'No limits found for specified CODE,FILTER!'
   specsnlimit1 = specsnlimit[i]
   fitmag = specsnlimit1.fitmag
   snmag = specsnlimit1.snmag
   slope = specsnlimit1.slope

   sigma = 0
   sn2 = 0
   nspec = n_elements(snvec)
   mask = (snvec GT 0 AND mag GT fitmag[0] AND mag LT fitmag[1])

   igood = where(mask, ngood)
   if (ngood LE 2) then return, 0

   logsn = snvec*0.0 - 1.0 ; Arbitrarily set bad values to -1, though these
                        ; values are masked from the fit anyway
   logsn[igood] = alog10(snvec[igood])

   if (ngood GE 3) then begin
      coeffs = [median(logsn[igood] - slope * mag[igood]), slope]
      yfit = coeffs[0] + coeffs[1] * mag
      sigma = djsig(logsn[igood] - yfit[igood], sigrej=sigrej)
      sn2 = 10^(2.0*poly(snmag,coeffs))
   endif else begin
      coeffs = 0
   endelse

;   for i=0, maxiter-1 do begin
;      igood = where(mask, ngood)
;      if (ngood LT 3) then return, 0
;      coeffs = poly_fit(mag[igood], logsn[igood], 1, /double)
;      yfit = poly(mag, coeffs)
;      sn2 = 10^(2.0*poly(snmag,coeffs))
; 
;      diff = logsn - yfit
;      djs_iterstat, diff[igood], sigrej=sigrej, sigma=sigma, mask=smask
;      treject = total(1-smask)
;      if (treject EQ 0) then return, coeffs
;
;      mask[igood] = mask[igood] * smask
;      if (total(mask) LE 2) then return, 0
;   endfor

   ;----------
   ; Correct sn2 for galactic dust reddening to match quick SOS reductions
   ; Recalculate reddening since information isn't propagated through plugmap.
   ; NOTE: These constants are also hardwired in design_plate.pro
   ;       If you change them here, also change them there
   if keyword_set(plugmap) then begin
      indx = where(strtrim(plugmap.holetype,2) EQ 'OBJECT', nobj)
      euler, plugmap[indx].ra, plugmap[indx].dec, ll, bb, 1
      reddenvec = [5.155, 3.793, 2.751, 2.086, 1.479] $
       * median(dust_getval(ll, bb, /interp))
      case filter of
         'u' : dered_sn2 = sn2 + coeffs[1]*reddenvec[0]
         'g' : dered_sn2 = sn2 + coeffs[1]*reddenvec[1]
         'r' : dered_sn2 = sn2 + coeffs[1]*reddenvec[2]
         'i' : dered_sn2 = sn2 + coeffs[1]*reddenvec[3]
         'z' : dered_sn2 = sn2 + coeffs[1]*reddenvec[4]
      endcase
   endif

   return, coeffs
end
;------------------------------------------------------------------------------
