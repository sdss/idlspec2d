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
   ; Determine the range of wavelengths, [LOGMIN,LOGMAX] in common w/all fibers

   if (loglam[1,0] GT loglam[0,0]) then begin ; Ascending wavelengths
      logmin = max(loglam[0,igood])
      logmax = min(loglam[ny-1,igood])
   endif else begin ; Descending wavelengths
      logmin = max(loglam[ny-1,igood])
      logmax = min(loglam[0,igood])
   endelse
   if (logmin GE logmax) then begin
      splog, 'WARNING: Wavelength mapping makes no sense!'
      return, 0
   endif

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
        median([ tmpflux[ (filtsz-1)/2 : ntmp-(filtsz-1)/2 ], filtsz ]) ) $
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
   if (keyword_set(x2)) then allx2 = (x2[*,igood])[isort]
   indx = where(flux[*,igood] GT minval)
   if (keyword_set(x2)) then thisx2 = allx2[indx]
   if (indx[0] EQ -1) then $
    message, 'No points above MINVAL'

; THE BELOW FAILS WITH NORD=4, SO MUST SET NORD=3!!!???
; ALSO, ONLY WORKS WITH NPOLY=1,2 OR 4!!!???
;   sset = bspline_iterfit(allwave[indx], allflux[indx], $
;    invvar=allivar[indx], x2=thisx2, nord=nord, npoly=npoly, nbkpts=ny, $
;    maxiter=maxiter, upper=upper, lower=lower, mask=mask)

;
;  David:  It's better to use "everyn=" instead or "nbkpts="
;          Due to the sparse sampling near the ends of wavelength range
;           you have effectively a higher bkpt density and this can give
;           too much freedom to very few data points.  
;          everyn places a breakpoint at every Nth good data point.
;          For the example below everyn=ngood gives about 2040 fullbkpts
;           instead of 2052 with nbkpts=ny
;          We use everyn= in skysubtraction as well, just because there
;           is always one or two fibers extending much further than the rest.
;
   sset = bspline_iterfit(allwave[indx], allflux[indx], $
    invvar=allivar[indx], x2=thisx2, nord=nord, npoly=npoly, everyn=ngood, $
    maxiter=maxiter, upper=upper, lower=lower, outmask=mask, requiren=2)

;   generate model fit for full frame
;   yy = bspline_valu(loglam, sset, x2=x2)

;;--------------------
; De-bugging tests...
;stop
;sset2 = sset ; 2D-fit
;sset1 = bspline_iterfit(allwave[indx], allflux[indx], $
; invvar=allivar[indx], nord=4, nbkpts=ny, $
; maxiter=maxiter, upper=upper, lower=lower, mask=mask2)
;ymodel1 = bspline_valu(allwave, sset1)
;ymodel2 = bspline_valu(allwave, sset2, x2=allx2)
;jj=indx[long(randomu(123,10000)*640000)] ; random sampling of points
;jj=jj[sort(allwave[jj])]
;rmap1 = allflux[jj] - ymodel1[jj]
;rmap2 = allflux[jj] - ymodel2[jj]
;splot,allwave[jj],rmap1,ps=3
;splot,allwave[jj],rmap2,ps=3
;;--------------------
;jj=indx[long(randomu(123,10000)*640000)] ; random sampling of points
;jj=jj[sort(allwave[jj])]
;sset1 = bspline_iterfit(allwave[jj], allflux[jj], $
; invvar=allivar[jj], nord=3, nbkpts=ny, $
; maxiter=maxiter, upper=upper, lower=lower, mask=mask)
;sset2 = bspline_iterfit(allwave[jj], allflux[jj], $
; invvar=allivar[jj], x2=allx2[jj], nord=3, npoly=2, nbkpts=ny, $
; maxiter=maxiter, upper=upper, lower=lower, mask=mask)
;sset3 = bspline_iterfit(allwave[jj], allflux[jj], $
; invvar=allivar[jj], x2=allx2[jj], nord=3, npoly=3, nbkpts=ny, $
; maxiter=maxiter, upper=upper, lower=lower, mask=mask)
;ymodel1 = bspline_valu(allwave, sset1)
;ymodel2 = bspline_valu(allwave, sset2, x2=allx2)
;ymodel3 = bspline_valu(allwave, sset2, x2=allx2)
;rmap1 = allflux[jj] - ymodel1[jj]
;rmap2 = allflux[jj] - ymodel2[jj]
;rmap3 = allflux[jj] - ymodel3[jj]
;splot,allwave[jj],rmap1,ps=3,yr=[-1,1]/10.
;splot,allwave[jj],rmap2,ps=3,yr=[-1,1]/10.
;splot,allwave[jj],rmap3,ps=3,yr=[-1,1]/10.
;print,djsig(rmap1),djsig(rmap2),djsig(rmap3)


; Should move this plotting elsewhere ???
   ;------
   ; QA plot of superflat   ; Plot sampled every 1 Ang

   if (keyword_set(title)) then begin
      wmin = fix(10^min(allwave))
      wmax = ceil(10^max(allwave))
      plot_lam = wmin + lindgen(wmax-wmin+1)
      if keyword_set(allx2) then begin
        plot_x2 = 0*plot_lam + median(allx2)
        plot_fit  = bspline_valu(alog10(plot_lam), sset, x2=plot_x2)
      endif else plot_fit  = bspline_valu(alog10(plot_lam), sset)

      djs_plot, plot_lam, plot_fit, xrange=[wmin,wmax], xstyle=1, $
       xtitle='\lambda [A]', ytitle='Normalized flux', $
       title=title

      ; Overplot pixels masked from the fit
      ii = where(mask EQ 0)
      if (ii[0] NE -1) then $
       djs_oplot, 10^allwave[indx[ii]], allflux[indx[ii]], ps=3, color='red'
   endif

   return, sset
end
;------------------------------------------------------------------------------
