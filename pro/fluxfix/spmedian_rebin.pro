;+
; NAME:
;   spmedian_rebin
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   smearcorr, bsmearfile, rsmearfile, bcalibfile, rcalibfile, $
;              bscifile, rscifile, corrfile, noplot = noplot
;
; INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:  Files containing the polynomial coefficients of the fit to 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Oct-2000  Formerly fluxcorr_new -- written by S. Burles
;-
;------------------------------------------------------------------------------
function spmedian_rebin, loglam, flux, ivar, color, mask=mask, sigrej=sigrej, $
         outwave = outwave, sn = sn, quality = quality

   if NOT keyword_set(sigrej) then sigrej = 20.0

   ;---------
   ; Define new wavelength sampling

   if color eq 'r' then begin
     w1 = findgen(54)*4.0e-3 + 3.756
     w2 = findgen(54)*4.0e-3 + w1[1]
   endif 
   if color eq 'b' then begin
     w1 = findgen(60)*4.0e-3 + 3.568
     w2 = findgen(60)*4.0e-3 + w1[1]
   endif
   if color eq 'full' then begin
     w1 = findgen(96)*4.0e-3 + 3.568
     w2 = findgen(96)*4.0e-3 + w1[1]
   endif
   range = transpose([[w1],[w2]])
   outwave = djs_median(range,1)

   ;--------------------------
   ; Set up vectors for output
   nr = (size(range))[2]
   ntrace = (size(flux))[2]
   fit = fltarr(nr, ntrace)
   mask = fltarr(nr, ntrace)

   ;--------------------------
   ; Do iterative rejection on each bin of each fiber
   for itrace=0,ntrace-1 do begin
     for irange = 0, nr -1 do begin
        inside = where(loglam[*,itrace] GE range[0,irange] $
                 AND loglam[*,itrace] LT range[1,irange], ninside)
        if ninside GT 0 then begin
           good = where(ivar[inside,itrace] GT 0, ngood)
       
           if ngood GT 1 then begin 
              djs_iterstat, flux[inside[good],itrace], median=md, sigma=sig
              fit[irange, itrace] = md
              if ngood GT 0.5 * ninside then $
                   mask[irange, itrace]  = total(ivar[inside[good],itrace])
              sn = sig * sqrt(mask[irange, itrace])
              if sn GT sigrej OR sn LE 0 then mask[irange, itrace] = 0.0
           endif
        endif
     endfor
   endfor

   sn = djs_median(flux * sqrt(ivar),1)

   badval = ( fibermask_bits('NOPLUG')   OR fibermask_bits('BADTRACE') $
       OR fibermask_bits('BADFLAT')  OR fibermask_bits('BADARC')   $
       OR fibermask_bits('NEARWHOPPER') )

   quality = (mask[0,*] AND badval)

   return, fit
end

;------------------------------------------------------------------------------
