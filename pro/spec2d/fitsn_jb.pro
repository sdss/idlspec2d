;+
; NAME:
;   fitsn_jb
;
; PURPOSE:
;   Perform a two parameter fit to S/N vs magnitude
;
; CALLING SEQUENCE:
;   coeffs = fitsn_jb(mag, snvec, [ sigrej=, maxiter=, redden=, sncode=, $
;    filter=, sigma=, specsnlimit=, sn2=, dered_sn2= ] )
;
; INPUTS:
;   mag        - Fiber magnitudes
;   snvec      - S/N vector for fibers
;   sncode     - Pipeline code for determining fit range
;   filter     - Filter band for determining fit range; only works for
;                ugriz filters
;
; OPTIONAL KEYWORDS:
;   sigrej     - Sigma rejection threshold; default to 3
;                (only used in computing SIGMA)
;   maxiter    - Maximum number of rejection iterations; default to 5
;                (no longer used with median-fit line)
;   redden     - 5-element array with median extinction in all bands
;                for use in computing DERED_SN2
;
; OUTPUTS:
;   coeffs     - Coefficients from fit; return 0 if fit failed;
;                flux = 10^(0.4*(22.5-mag))
;                S/N = coeffs[0] * flux / sqrt( flux + coeffs[1])
;
; OPTIONAL OUTPUTS:
;   sigma      - Standard deviation of residuals
;   specsnlimit- Structure with FITMAG,SNMAG,FIDUCIAL_COEFF
;   sn2        - Fit value of (S/N)^2 at fiducial magnitude
;   dered_sn2  - Extinction corrected (S/N)^2 at fiducial magnitude
;                (requires REDDEN as input)
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
;   general_sn
;
; REVISION HISTORY:
;   15-Apr-2000  Written by S. Burles, FNAL
;   12-Jun-2015  Realistic model by J. Bautista
;-
;------------------------------------------------------------------------------
function fitsn_jb, mag, snvec, sigrej=sigrej, maxiter=maxiter, redden=redden, $
 sncode=sncode, filter=filter, sigma=sigma, specsnlimit=specsnlimit1, $
 sn2=sn2, dered_sn2=dered_sn2

   if (NOT keyword_set(sigrej)) then sigrej = 3

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
   dered_sn2 = 0
   nspec = n_elements(snvec)
   ; JEB - avoiding too low SN2 in fit
   mask = (snvec GT 0. AND mag GT fitmag[0] AND mag LT fitmag[1])

   igood = where(mask, ngood)
   splog, 'Default fit range contains ', ngood, ' values'

   ;-- JEB: using full range of mags, there should be more than 3 points

   if (ngood LE 2) then return, 0

   if (ngood GE 3) then begin
      ;-- JEB: converting back to flux
      x = 10^(0.4*(22.5-mag[igood]))
      y = snvec[igood]
      ;-- JEB: uniform weighting is fine
      erry = replicate(1., ngood)
      guess = [ 3.7, 10.]
      ;-- JEB: limiting sky term to positive values
      parsinfo = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},n_elements(guess))
      parsinfo(1).limited(0) = 1 
      parsinfo(1).limits(0) = 1.
      coeffs = MPFITFUN('general_sn', x, y, yerr, guess, status=status, parinfo=parsinfo) 
   endif 

   if status GT 0 then begin
      splog, 'Best fit parameters', coeffs
      yfit = general_sn(x, coeffs) 
      sigma = djsig(alog10(y) - alog10(yfit), sigrej=sigrej)
      sn2 = general_sn( 10^(0.4*(22.5-snmag)),coeffs)^2
   endif else begin
      splog, 'FITSN_JB FAILED'
      coeffs = 0
   endelse

   ;----------
   ; Correct sn2 for galactic dust reddening to match quick SOS reductions
   ; NOTE: These constants are also hardwired in design_plate.pro
   ;       If you change them here, also change them there
   if keyword_set(redden) then begin
      splog, 'Computing with REDDEN=', redden, format='(a,5f7.3)'
      case filter of
         'u' : dered_sn2 = sn2 ;+ coeffs[1]*redden[0]    ; JEB fit already done in 
         'g' : dered_sn2 = sn2 ;+ coeffs[1]*redden[1]    ; deredden mags
         'r' : dered_sn2 = sn2 ;+ coeffs[1]*redden[2]
         'i' : dered_sn2 = sn2 ;+ coeffs[1]*redden[3]
         'z' : dered_sn2 = sn2 ;+ coeffs[1]*redden[4]
      endcase
   endif

   return, coeffs
end
;------------------------------------------------------------------------------
