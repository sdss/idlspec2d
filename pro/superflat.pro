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
;   fibermask  - Mask of 0 for bad fibers and 1 for good fibers [NFIBER]
;   minval     - Minimum value to use in fits to flat-field vectors;
;                default to 0.03.
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
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
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

   if (N_elements(fibermask) NE ntrace) then fibermask = bytarr(ntrace) + 1
   igood = where(fibermask NE 0, ngood)

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
   ; Do this with a straight median value for all wavelengths in common

   qq = loglam GE logmin AND loglam LE logmax
   medval = fltarr(ntrace)
   for i=0, ntrace-1 do $
    medval[i] = median( flux[where(qq[*,i]),i] )
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
   indx = where(allflux GT minval)
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
   djs_oplot, 10^allwave[indx[ii]], allflux[indx[ii]], ps=3, color='red'

   return
end

;------------------------------------------------------------------------------
