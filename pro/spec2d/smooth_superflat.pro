;+
; NAME:
;   smooth_superflat
;
; PURPOSE:
;   Take the superflat fit and target wavelengths and filter superflat
;    to be sure to remove spurious features
;
; CALLING SEQUENCE:
;   smooth_fit = smooth_superflat( superflatset, airset )
;
; INPUTS:
;   superflatset - Superflat bsplineset returned from "superflat"
;   airset       - Wavelength solution (preferably shifted to match skylines)
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   smooth_fit   - Filtered superflat fit smoothed over about 4 pixels
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
;   bspline_iterfit()
;   bspline_valu()
;
; REVISION HISTORY:
;   27-Jul-2001  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function smooth_superflat, superflatset, airset

   if N_PARAMS() LT 2 then return, 0

   traceset2xy, airset, pixnorm, loglam

   ; sparse sample loglam
   npix    = (size(loglam))[1]
   nfibers = (size(loglam))[2]
   ntotal  = n_elements(loglam)

   ; evaluate at every half-pixel...

   nsparse = npix * 2 + 1
   sparselam = (loglam[sort(loglam)])[lindgen(nsparse)*nfibers/2 < (ntotal-1)]
   model = bspline_valu(sparselam, superflatset)

   ; fit again with fewer breakpoints, one every 4 pixels (8 half pixels)
   ;  and reject with impunity

   smoothset = bspline_iterfit(sparselam, model, invvar=model*0 + 5.0e5, $
                  everyn = 8, lower=3, upper=2, yfit=yfit, outmask=outmask)
 
   bad = where(outmask EQ 0, nbad)

   ; if pixels were rejected (especially emission lines), then grow one pixel
   ; and refit.

   if nbad GE 1 then begin 
     inmask =  outmask
     inmask[bad - 1 > 0] = 0
     inmask[bad + 1 < (nsparse - 1)] = 0
     smoothset = bspline_iterfit(sparselam, model, invvar=inmask, $
                        everyn = 8, yfit=yfit)
   endif

   fullfit = bspline_valu(loglam, smoothset)

   return, fullfit
end


   

