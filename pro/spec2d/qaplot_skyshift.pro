;+
; NAME:
;   qaplot_skyshift
;
; PURPOSE:
;   Generate QA plot: scatter in sky line positions
;
; CALLING SEQUENCE:
;   qaplot_skyshift, wset, xsky, skywaves, skyshift, [title= ]
;
; INPUTS:
;   wset       - Wavelength solution from arc line spectra
;   xsky       - Pixel position of sky lines [NFIBER,NLINE]
;   skywaves   - Wavelengths of sky lines used for computing shifts (Ang)
;   skyshift   - Shifts to apply to sky lines in pix [NROW,NTRACE]
;
; OPTIONAL KEYWORDS:
;   title      - TITLE of plot
;
; OUTPUTS:
;   Output plots only
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
;   djs_laxisgen()
;   djs_oplot
;   djs_plot
;   traceset2pix
;
; REVISION HISTORY:
;   11-Apr-2000 Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro qaplot_skyshift, wset, xsky, skywaves, skyshift, title=title

   dims = size(xsky, /dimens)
   nfiber = dims[0]

   nskyline = n_elements(skywaves)
   ysky = dindgen(nfiber) # (dblarr(nskyline)+1)
   xpredict = transpose( traceset2pix(wset, alog10(skywaves)) )

   scalefac = 1.0 ; Scale factor in pixel shift per tic mark on the X axis

   plot, [0], [0], /nodata, xrange=[-2.5,nskyline+0.5], yrange=[0,nfiber], $
    xtitle='Sky Line Wavelength', ytitle='Fiber Number', $
    xticks=nskyline-1, xtickv=indgen(nskyline), $
    xtickname=string(fix(skywaves),format='(i4)'), $
    title=title, charsize=0.8, /xstyle, /ystyle

   ;----------
   ; Plot zero lines for each sky wavelength

   xplot = djs_laxisgen([nfiber,nskyline],iaxis=1)
   yplot = djs_laxisgen([nfiber,nskyline],iaxis=0)
   for i=0, nskyline-1 do $
    djs_oplot, xplot[*,i], yplot[*,i]

   ;----------
   ; Plot sky line positions before shift

   xplot2 = xplot + (xsky - xpredict) / scalefac
   for i=0, nskyline-1 do $
    djs_oplot, xplot2[*,i], yplot[*,i], color='red'

   ;----------
   ; Plot sky line positions after shift

   xplot2 = xplot + (xsky - xpredict - skyshift) / scalefac
   for i=0, nskyline-1 do $
    djs_oplot, xplot2[*,i], yplot[*,i], color='green'

   ;----------
   ; Label plot

   djs_oplot, [-2.0,-1.0], [1,1]*nfiber*0.90
   xyouts, -2.0, nfiber*0.85, string(format='(f4.2, " pix")', scalefac)

   djs_oplot, [-2.0,-2.0], nfiber*[0.50,0.62], color='red'
   xyouts, -2.0, nfiber*0.65, 'Before shift', orientation=90

   djs_oplot, [-1.0,-1.0], nfiber*[0.50,0.62], color='green'
   xyouts, -1.0, nfiber*0.65, 'After shift', orientation=90

   return
end
