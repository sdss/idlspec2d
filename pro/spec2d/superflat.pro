;+
; NAME:
;   superflat
;
; PURPOSE:
;   Create a "superflat" from an extracted flat-field image
;
; CALLING SEQUENCE:
;   superflat, flux, fluxivar, wset, fullbkpt, coeff, $
;    [ fibermask=, minval=, lower=, upper=, medval= ]
;
; INPUTS:
;   flux       - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   fluxivar   - Inverse variance map for FLUX.
;   wset       - Wavelength solution
;
; OPTIONAL KEYWORDS:
;   fibermask  - Mask of 0 for good fibers [NFIBER]
;   minval     - Minimum value to use in fits to flat-field vectors;
;                default to 0.
;
; PARAMETERS FOR SLATEC_SPLINEFIT:
;   lower      -
;   upper      -
;
; OUTPUTS:
;   fullbkpt   - Breakpoints describing spline fit to superflat
;   coeff      - Coefficients describing spline fit to superflat
;
; OPTIONAL OUTPUTS:
;   medval     - Median value of each fiber [NFIBER]
;   fibermask  - (Modified)
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_maskinterp()
;   djs_mean()
;   djs_oplot
;   djs_plot
;   slatec_bvalu()
;   slatec_splinefit()
;   traceset2xy
;
; REVISION HISTORY:
;   02-Jan-2000  Excised code from SPFLATTEN2 (DJS).
;-
;------------------------------------------------------------------------------

pro superflat, flux, fluxivar, wset, fullbkpt, coeff, $
 fibermask=fibermask, minval=minval, lower=lower, upper=upper, medval=medval

   if (NOT keyword_set(minval)) then minval = 0.0

   dims = size(flux, /dimens)
   ny = dims[0]
   ntrace = dims[1]

   if (N_elements(fibermask) NE ntrace) then fibermask = bytarr(ntrace) 
   igood = where(fibermask EQ 0, ngood)

   ;------
   ; Determine LOGLAM from the wavelength solution

   traceset2xy, wset, xx, loglam

   ;------
   ; Determine the range of wavelengths, [LOGMIN,LOGMAX] in common w/all fibers

   if (loglam[1,0] GT loglam[0,0]) then begin ; Ascending wavelengths
      logmin = max(loglam[0,igood])
      logmax = min(loglam[ny-1,igood])
   endif else begin ; Descending wavelengths
      logmin = max(loglam[ny-1,igood])
      logmax = min(loglam[0,igood])
   endelse

   ;------
   ; Find the approximate scalings between all fibers
   ; Do this with a straight mean value for all wavelengths in common,
   ; interpolating over bad pixels.
   ;   FRACPTS = Fraction of unmasked pixels in each flat vector
   ;   MEDVAL = Mean value for each flat vector, after a median-filter
   ;            which hopefully removes cosmic rays and the like

   filtsz = 11
   qq = loglam GE logmin AND loglam LE logmax
   medval = fltarr(ntrace)
   fracpts = fltarr(ntrace)
   for i=0, ntrace-1 do begin
      indx = where(qq[*,i], ntmp)
      tmpmask = fluxivar[indx,i] EQ 0
      tmpflux = djs_maskinterp( flux[indx,i], tmpmask )
      fracpts[i] = 1.0 - total(tmpmask)/N_elements(tmpmask)
      if (ntmp GT filtsz) then $
       medval[i] = djs_mean( $
        median( tmpflux[ (filtsz-1)/2 : ntmp-(filtsz-1)/2 ], filtsz ) ) $
      else $
       medval[i] = median(tmpflux)
   endfor

   ;------
   ; Limit the superflat to use only good fibers, and those fibers that
   ; have at least 95% good wavelength range as compared to the best fiber
   ; and whose counts are within 30% of the median (good) fiber throughput.

   globalmed = median(medval[igood])
   if (globalmed LT 0) then $
    message, 'Median flat-field vector is negative!'
   igood = where(fibermask EQ 0 AND fracpts GE 0.95*max(fracpts) $
    AND abs(medval-globalmed)/globalmed LT 0.30, ngood)
; ??? Should we set a bit in FIBERMASK for fibers unused for the superflat ???

   ;-----
   ; Prevent divide-by-zeros below

   izero = where(medval LE 0)
   if (izero[0] NE -1) then medval[izero] = 1.0

   ;------
   ; Create a version of flux (and fluxivar) that has all fibers
   ; approximately scaled to have a median value of 1

   scalef = fltarr(ny,ntrace)
   scalefivar = fltarr(ny,ntrace)
   for i=0, ntrace-1 do $
    scalef[*,i] = flux[*,i] / medval[i]
   for i=0, ntrace-1 do $
    scalefivar[*,i] = fluxivar[*,i] * (medval[i])^2

   ;------
   ; Create a "superflat" spectrum, analogous to the "supersky"

   splog, 'Creating superflat from ', ngood, ' fibers'
   isort = sort(loglam[*,igood])
   allwave = (loglam[*,igood])[isort]
   allflux = (scalef[*,igood])[isort]
   allivar = (scalefivar[*,igood])[isort]
   indx = where(flux[*,igood] GT minval)
   if (indx[0] EQ -1) then $
    message, 'No points above MINVAL'
   fullbkpt = slatec_splinefit(allwave[indx], allflux[indx], coeff, $
    maxiter=maxiter, upper=upper, lower=lower, $
    invvar=allivar[indx], nord=4, nbkpts=ny, mask=mask)

; Should move this plotting elsewhere ???
   ;------
   ; QA plot of superflat   ; Plot sampled every 1 Ang

   wmin = fix(10^min(allwave))
   wmax = ceil(10^max(allwave))
   plot_lam = wmin + lindgen(wmax-wmin+1)
   plot_fit  = slatec_bvalu(alog10(plot_lam), fullbkpt, coeff)
   djs_plot, plot_lam, plot_fit, xrange=[wmin,wmax], xstyle=1, $
    xtitle='\lambda [A]', ytitle='Normalized flux', $
    title='Superflat for ???'

   ; Overplot pixels masked from the fit
   ii = where(mask EQ 0)
   if (ii[0] EQ -1) then splog, 'No pixels masked from fit, YIPPEE!!???' $
   else djs_oplot, 10^allwave[indx[ii]], allflux[indx[ii]], ps=3, color='red'

   return
end

;------------------------------------------------------------------------------
