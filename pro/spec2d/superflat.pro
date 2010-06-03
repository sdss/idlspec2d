;+
; NAME:
;   superflat
;
; PURPOSE:
;   Create a "superflat" from an extracted flat-field image
;
; CALLING SEQUENCE:
;   sset = superflat(flux, fluxivar, wset, fullbkpt, coeff, $
;    [ fibermask=, minval=, medval=, title=, $
;    x2=, nord=, npoly=, upper=, lower= ])
;
; INPUTS:
;   flux       - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   fluxivar   - Inverse variance map for FLUX.
;   wset       - Wavelength solution
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   minval     - Minimum value to use in fits to flat-field vectors;
;                default to 0.
;   title      - TITLE of plot; if not set, then do not make this plot.
;
; OPTIONAL PARAMETERS FOR BSPLINE_ITERFIT:
;   x2         - Orthogonal dependent variable for B-spline fit;
;                this will typically be the X position on the CCD.
;   nord       - Order of b-splines; default of 4 for cubic
;   npoly      - Order of X2 polynomial fit; default to 0 for none
;   lower      -
;   upper      -
;
; OUTPUTS:
;   sset       - Output structure describing spline fit to superflat.
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
;   bspline_valu()
;   bspline_iterfit()
;   traceset2xy
;
; REVISION HISTORY:
;   02-Jan-2000  Excised code from SPFLATTEN2 (DJS).
;-
;------------------------------------------------------------------------------
function superflat, flux, fluxivar, wset, x2=x2, $
 fibermask=fibermask, minval=minval, medval=medval, title=title, $
 nord=nord, npoly=npoly, upper=upper, lower=lower

   if (NOT keyword_set(minval)) then minval = 0.0
   if (NOT keyword_set(nord)) then nord = 4

   dims = size(flux, /dimens)
   ny = dims[0]
   ntrace = dims[1]

   if (N_elements(fibermask) NE ntrace) then fibermask = bytarr(ntrace) 
   igood = where(fibermask EQ 0, ngood)

   ;------
   ; Determine LOGLAM from the wavelength solution

   traceset2xy, wset, xx, loglam

   ;------
   ; Determine the range of wavelengths, [LOGMIN,LOGMAX]
   ; There must be at least one unmasked pixel at each included wavelengths,
   ; and it must at least fall on the CCD for all good fibers.

   goodpix = where((fluxivar GT 0) $
    AND transpose(rebin(fibermask,ntrace,ny)) EQ 0)
   logmin = min(loglam[goodpix], max=logmax)
   if (loglam[1,0] GT loglam[0,0]) then begin ; Ascending wavelengths
      logmin = logmin > max(loglam[0,igood])
      logmax = logmax < min(loglam[ny-1,igood])
   endif else begin ; Descending wavelengths
      logmin = logmin > max(loglam[ny-1,igood])
      logmax = logmax < min(loglam[0,igood])
   endelse
   if (logmin GE logmax) then begin
      splog, 'WARNING: Wavelength mapping makes no sense!'
      return, 0
   endif
   splog, 'Wavelength range = ', 10.^logmin, 10.^logmax

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
       medval[i] = $
        median([ median(tmpflux[ (filtsz-1)/2 : ntmp-(filtsz-1)/2 ], filtsz) ]) $
      else $
       medval[i] = median([tmpflux])
   endfor

   ;------
   ; Limit the superflat to use only good fibers, and those fibers that
   ; have at least 95% good wavelength range as compared to the best fiber
   ; and whose counts are within 30% of the median (good) fiber throughput.

   globalmed = median([medval[igood]])
   if (globalmed LT 0) then $
    message, 'Median flat-field vector is negative!'
   igood = where(fibermask EQ 0 $
    AND fracpts GE 0.95*max(fracpts[where(fibermask EQ 0)]) $
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
   if (keyword_set(x2)) then allx2 = (x2[*,igood])[isort]
   indx = where(flux[*,igood] GT minval)
   if (keyword_set(x2)) then thisx2 = allx2[indx]
   if (indx[0] EQ -1) then $
    message, 'No points above MINVAL'

   sset = bspline_iterfit(allwave[indx], allflux[indx], $
    invvar=allivar[indx], x2=thisx2, nord=nord, npoly=npoly, everyn=ngood, $
    maxiter=maxiter, upper=upper, lower=lower, outmask=mask, requiren=2)

;   generate model fit for full frame
;   yy = bspline_valu(loglam, sset, x2=x2)

; Should move this plotting elsewhere ???
   ;------
   ; QA plot of superflat   ; Plot sampled every 1 Ang

   if (keyword_set(title)) then begin
      ; Interpolate all the spectra to the wavelength mapping
      ; of the central fiber
      wplot = loglam[*,ntrace/2]
      fplot = fltarr(ny,ntrace)
      for i=0L, ntrace-1L do begin
         ii = where(fluxivar[*,i] GT 0, ct)
         if (ct GT 0) then $
          fplot[*,i] = interpol(flux[*,i], loglam[*,i], wplot)
      endfor

      if (keyword_set(allx2)) then begin
        plot_x2 = x2[*,ntrace/2]
        plot_fit  = bspline_valu(wplot, sset, x2=plot_x2)
      endif else plot_fit  = bspline_valu(wplot, sset)

      percentiles = [0.023,0.158,0.50,0.842,0.977]
      fmed = fltarr(ny, n_elements(percentiles))
      for j=0L, ny-1L do begin
         isort = sort(fplot[j,*])
         fmed[j,*] = fplot[j,isort[percentiles*ntrace]]
      endfor

      wrange = minmax( $
       sset.fullbkpt[sset.nord-1:n_elements(sset.fullbkpt)-sset.nord+1] )
      ii = where(wplot GT wrange[0] AND wplot LT wrange[1])
      fmed *= mean(plot_fit[ii]) / mean(fmed[ii,(n_elements(percentiles)-1)/2])
      yrange = [0, 1.25 * weighted_quantile(plot_fit[ii], quant=0.90)]

      djs_plot, 10.^wrange, yrange, /nodata, xrange=10.^wrange, xstyle=1, $
       yrange=yrange, ystyle=1, thick=2, $
       xtitle='\lambda [A]', ytitle='Normalized flux', $
       title=title
      for k=0,n_elements(percentiles)-1 do $
       djs_oplot, 10.^wplot, fmed[*,k]
      djs_oplot, 10.^wplot, plot_fit, color='red'

      djs_xyouts, total([0.95,0.05]*!x.crange), total([0.05,0.95]*!y.crange), $
       'RED = Superflat fit to central fiber', color='red'
      djs_xyouts, total([0.95,0.05]*!x.crange), total([0.10,0.90]*!y.crange), $
       'BLACK =' + string(percentiles*100, $
       format='('+string(n_elements(percentiles))+'f6.1)')+' percentiles'
   endif

   return, sset
end
;------------------------------------------------------------------------------
