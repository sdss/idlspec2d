;+
; NAME:
;   qaplot_fflat
;
; PURPOSE:
;   Generate QA plot for fiber-flats
;
; CALLING SEQUENCE:
;   qaplot_fflat, fflat, wset, [ fibermask=, dx=, filename= ]
;
; INPUTS:
;   fflat      - Array of flat-field flat-field vectors for each fiber
;                that remove relative flat-field variations as a function
;                of wavelength between fibers [Nrow, Ntrace]
;   wset       - Wavelength solution
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   dx         - Bin log-lambda by this number; default to 1.e-3 (about 10 pix)
;   filename   - File name to use for TITLE of plot
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Do not plot bad fibers or masked pixels in good fibers.
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

pro qaplot_fflat, fflat, wset, fibermask=fibermask, dx=dx, filename=filename

   if (NOT keyword_set(dx)) then dx = 1.e-3
   if (NOT keyword_set(filename)) then filename = ''

   ;----------
   ; Compute the wavelengths for all flat vectors from the trace set.
   ; Set these equal to zero for bad fibers, and ignore in these plots.

   traceset2xy, wset, xx, loglam
   xrange = [min(loglam), max(loglam)]

   if (keyword_set(fibermask)) then begin
      indx = where(fibermask AND fibermask_bits('BADFLAT'))
      if (indx[0] NE -1) then loglam[*,indx] = 0
   endif

   dims = size(fflat, /dimens)
   ny = dims[0]
   nfiber = dims[1]

   ;----------
   ; Strip any trailing ".fit" from the file name

   i = rstrpos(filename, '.fit')
   if (i GT 0) then trimname = strmid(filename, 0, i) $
    else trimname = filename

   ; Divide the plotting range up into bins in LOGLAM separated by DX
   nbin = fix( (xrange[1]-xrange[0]) / dx ) + 1
   meanval = fltarr(nbin)
   botval = fltarr(nbin)
   topval = fltarr(nbin)
   ioutlier = [0]
   for ibin=0, nbin-1 do begin
      j = where(loglam GE xrange[0]+dx*ibin $
            AND loglam LT xrange[0]+dx*(ibin+1) $
            AND fflat NE 0, ct)
      if (ct GT 0) then begin
         fsort = fflat[j]
         fsort = fsort[sort(fsort)]

         botval[ibin] = fsort[floor(0.01*ct)]
         topval[ibin] = fsort[ceil(0.99*ct)]

         meanval[ibin] = median([fflat[j]])

         iout = where(fflat[j] LT botval[ibin] OR fflat[j] GT topval[ibin])
         if (iout[0] NE -1) then ioutlier = [ioutlier, j[iout]]
      endif

      ; Burles counter of bin number...
      print, format='($, ".",i4.4,a5)', ibin, string([8b,8b,8b,8b,8b])
   endfor

   ; Set up the plot
   xaxis = xrange[0] + (findgen(nbin)+0.5) * dx
   djs_plot, xaxis, meanval, xrange=xrange, yrange=[-0.1,1.9], $
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
   errplot, xaxis, botval, topval

   xyouts, xaxis[fix(nbin/20)], 1.8, $
    'BLACK = Median value'
   xyouts, xaxis[fix(nbin/20)], 1.7, $
    'ENVELOPE = 1% to 99% range'
   xyouts, xaxis[fix(nbin/20)], 1.6, $
    'GREEN = Outlier points'
   xyouts, xaxis[fix(nbin/20)], 1.5, $
    'Min value =' + string(min(fflat[where(fflat NE 0)]),format='(f9.3)')
   xyouts, xaxis[fix(nbin/20)], 1.4, $
    'Max value =' + string(max(fflat[where(fflat NE 0)]),format='(f9.3)')
   xyouts, 0.95, 0., systime(), /normal, align=1 ; , charsize=0.5

   return
end
;------------------------------------------------------------------------------
