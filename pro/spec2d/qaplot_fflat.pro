;+
; NAME:
;   qaplot_fflat
;
; PURPOSE:
;   Generate QA plot for fiber-flats
;
; CALLING SEQUENCE:
;   qaplot_fflat, fflat, wset, [ fibermask=, plotsig=, dx=, filename= ]
;
; INPUTS:
;   fflat      - Array of flat-field flat-field vectors for each fiber
;                that remove relative flat-field variations as a function
;                of wavelength between fibers [Nrow, Ntrace]
;   wset       - Wavelength solution
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   plotsig    - Plot error bars out to this many standard deviations of
;                the fiber flats at each wavelength; only plot individual
;                outliers that are beyond this many deviations; default to 2.0
;   dx         - Bin log-lambda by this number; default to 1.e-3 (about 10 pix)
;   filename   - File name to use for TITLE of plot
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
;   djs_oplot
;   djs_plot
;
; REVISION HISTORY:
;   23-Nov-1999  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro qaplot_fflat, fflat, wset, fibermask=fibermask, $
 plotsig=plotsig, dx=dx, filename=filename

   if (NOT keyword_set(plotsig)) then plotsig = 2.0
   if (NOT keyword_set(dx)) then dx = 1.e-3
   if (NOT keyword_set(filename)) then filename = ''

   ; Compute the wavelengths for all flat vectors from the trace set
   traceset2xy, wset, xx, loglam

   dims = size(fflat, /dimens)
   ny = dims[0]
   nfiber = dims[1]
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber) + 1

   ; Strip any trailing ".fit" from the file name
   i = rstrpos(filename, '.fit')
   if (i GT 0) then trimname = strmid(filename, 0, i) $
    else trimname = filename

   xrange = [min(loglam), max(loglam)]

   ; Divide the plotting range up into bins in LOGLAM separated by DX
   nbin = fix( (xrange[1]-xrange[0]) / dx ) + 1
   meanval = fltarr(nbin)
   sigval = fltarr(nbin)
   ioutlier = [0]
   for ibin=0, nbin-1 do begin
      j = where(loglam GE xrange[0]+dx*ibin $
            AND loglam LT xrange[0]+dx*(ibin+1), ct)
      if (ct GT 0) then begin
         djs_iterstat, fflat[j], mean=mn, sigma=sg
         meanval[ibin] = mn
         sigval[ibin] = plotsig * sg ; Make the cut at this many sigma
         iout = where(fflat[j] LT meanval[ibin]-sigval[ibin] $
                   OR fflat[j] GT meanval[ibin]+sigval[ibin])
         if (iout[0] NE -1) then ioutlier = [ioutlier, j[iout]]
      endif

      ; Burles counter of bin number...
      print, format='($, ".",i4.4,a5)', ibin, string([8b,8b,8b,8b,8b])
   endfor

   ; Set up the plot
   xaxis = xrange[0] + (findgen(nbin)+0.5) * dx
   djs_plot, xaxis, meanval, /nodata, xrange=xrange, yrange=[-0.1,1.9], $
    xstyle=1, ystyle=1, $
    xtitle='log \lambda [A]', ytitle='Fiber-Flat', $
    title='Fiber-Flats for  '+trimname

   ; Overplot individual outlier points in green
   noutlier = N_elements(ioutlier) - 1
   if (noutlier GT 1) then begin ; The first outlier is a dummy set to [0]
      ioutlier = ioutlier[1:noutlier-1]
      djs_oplot, loglam[ioutlier], fflat[ioutlier], ps=3, color='green'
   endif

   ; Overplot the mean flat-field vector, and the dispersion
   errplot, xaxis, (meanval-sigval), (meanval+sigval)

   ; Overplot masked fibers from FIBERMASK in red
   ibfiber = where(fibermask EQ 0, nbfiber)
   if (nbfiber GT 0) then begin
      djs_oplot, loglam[*,ibfiber], fflat[*,ifiber], ps=3, color='red'
   endif

   xyouts, xaxis[fix(nbin/20)], 1.8, $
    'ENVELOPE = ' + string(plotsig,format='(f4.1)') + ' * stdev'
   xyouts, xaxis[fix(nbin/20)], 1.7, $
    'GREEN = Outlier points'
   xyouts, xaxis[fix(nbin/20)], 1.6, $
    'RED = Masked fibers'
   xyouts, xaxis[fix(nbin/20)], 1.5, $
    'Min value =' + string(min(fflat),format='(f9.3)')
   xyouts, xaxis[fix(nbin/20)], 1.4, $
    'Max value =' + string(max(fflat),format='(f9.3)')
   xyouts, 0.95, 0., systime(), /normal, align=1 ; , charsize=0.5

   return
end
;------------------------------------------------------------------------------
