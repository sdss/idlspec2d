;+
; NAME:
;   qaplot_skyline
;
; PURPOSE:
;   Generate QA plots for flux in particular skylines
;
; CALLING SEQUENCE:
;   qaplot_skyline, lwave, obj, objivar, objsub, objsubivar, plugsort, wset, $
;    [fibermask=, dwave=, title= ]
;
; INPUTS:
;   lwave      -
;   obj        - Image
;   objivar    - Inverse variance for OBJ
;   objsub     - Image after sky-subtraction
;   objsubivar - Inverse variance for image after sky-subtraction
;   plugsort   - Plugmap structure trimmed to one element per fiber
;   wset       - Wavelength solution
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   dwave      - Half-width about LWAVE for fitting sky line;
;                default to 5.0 Ang
;   titel      - File name to use for TITLE of plot
;
; OUTPUTS:
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
;   djs_iterstat
;   djs_median
;   djs_oplot
;   splog
;   traceset2xy
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   01-Jan-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro qaplot_skyline, lwave, obj, objivar, objsub, objsubivar, plugsort, wset, $
   iskies, fibermask=fibermask, dwave=dwave, title=title

   if (NOT keyword_set(title)) then title = ''
   if (NOT keyword_set(lwave)) then return
   if (NOT keyword_set(dwave)) then dwave = 5.0

   dims = size(objsub, /dimens)
   ncol = dims[0]
   nrow = dims[1]

   if (n_elements(fibermask) NE nrow) then fibermask = bytarr(nrow) 

   ;----------
   ; Solve for wavelength of each pixel 

   traceset2xy, wset, pixnorm, wave

   ;----------
   ; Find sky fibers

   nskies = n_elements(iskies)
   if (nskies EQ 0) then begin
      splog, 'No sky fibers!'
      return
   endif

   ;---------------------------------------------------------------------------
   ; Spline all the spectra near the specified wavelength

   lflux = fltarr(nrow)
   for i=0, nrow-1 do begin
      ipix = where( abs(10^wave[*,i]-lwave) LT dwave $
              AND objivar[*,i] GT 0, npix)
      if (npix GT 6 AND (fibermask[i] EQ 0)) then begin

         estimates = [16000.0, alog10(lwave), 1.0e-4, 500.0]

         bb = gaussfit(wave[ipix,i], obj[ipix,i], aa, nterms=4, $
          estimates=estimates)
         lflux[i] = total(bb - aa[3])

      endif
   endfor

;jflux = fltarr(nrow)
;for i=0, nrow-1 do begin
;   wd = abs(10^wave[*,i]-lwave)
;   ipix = where(wd LT dwave AND objivar[*,i] GT 0, npix)
;   iback = where(wd GE dwave AND wd LT 2*dwave AND objivar[*,i] GT 0, nback)
;   if (npix GT 6 AND nback GT 6 AND fibermask[i]) then begin
;      backflux = median(obj[iback,i])
;      jflux[i] = total(obj[ipix,i] - backflux)
;   endif
;endfor

   ; Compute the mean flux for this line in the sky fibers
   ii = where(lflux[iskies] GT 0)
   if (ii[0] EQ -1) then begin
      splog, 'WARNING: No good fits to this sky feature at ', lwave
      return
   endif
   djs_iterstat, lflux[iskies[ii]], mean=lmean, sigma=lsig
   splog, format='(a,f7.1,a,f6.2,a)', 'Flux dispersion for skyline ', $
    lwave, ' = ', 100*lsig/lmean, ' %'

   ;---------------------------------------------------------------------------

   ; Set multi-plot format
   !p.multi = [0,1,2]

   fibernum = indgen(nrow) + 1
   yrange = lmean+[-7,7]*lsig

   plot, fibernum, lflux, psym=1, $
    xrange=[0,nrow], xstyle=1, yrange=yrange, $
    xtitle='Fiber number', $
    ytitle=string('Flux at ', lwave, format="(a,f7.1)"), $
    title=title
   djs_oplot, fibernum[iskies], lflux[iskies], psym=2, color='red'
   djs_oplot, [0,nrow], [lmean,lmean], color='red'
   xyouts, 20, 0.80*!y.crange[0]+0.20*!y.crange[1], 'RED = sky fiber'
   xyouts, 20, 0.90*!y.crange[0]+0.10*!y.crange[1], $
    string(format='(a,f6.2,a)', 'Dispersion = ', 100*lsig/lmean, ' %')

   radius = sqrt(plugsort.xfocal^2 + plugsort.yfocal^2)
   plot, radius, lflux, psym=1, $
    yrange=yrange, $
    xtitle='Focal Distance [mm]', $
    ytitle=string('Flux at ', lwave, format="(a,f7.1)")
   djs_oplot, fibernum[iskies], lflux[iskies], psym=2, color='red'
   djs_oplot, !x.crange, [lmean,lmean], color='red'

   !p.multi = 0

   ;---------------------------------------------------------------------------

   return
end
;------------------------------------------------------------------------------
