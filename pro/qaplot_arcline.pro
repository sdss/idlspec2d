;+
; NAME:
;   qaplot_arcline
;
; PURPOSE:
;   Generate QA plot for arcline fits
;
; CALLING SEQUENCE:
;   qaplot_arcline, xdif, lambda, filename=filename, color=color
;
; INPUTS:
;   xdif       - Deviations in arc line fits in pixels, array [NFIBER, NLAMP]
;   lambda     - Wavelength for each lamp in Angstroms, vector [NLAMP]
;
; OPTIONAL KEYWORDS:
;   filename   - File name to use for TITLE of plot
;   color      - Specify 'blue' or 'red' to fix plotting limits
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
;   djs_oplot
;   djs_plot
;
; REVISION HISTORY:
;   15-Oct-1999  Written by D. Finkbeiner, APO
;   23-Nov-1999  Modified by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro qaplot_arcline, xdif, lambda, filename=filename, color=color

   if (NOT keyword_set(filename)) then filename = ''

   ; Set multi-plot format
   pmulti = !p.multi
   !p.multi = [0,1,2]

   dims = size(xdif, /dimens)
   nfiber = dims[0]
   nlamp = dims[1]

   ; Strip any trailing ".fit" from the file name
   i = rstrpos(filename, '.fit')
   if (i GT 0) then trimname = strmid(filename, 0, i) $
    else trimname = filename

   ; Determine the plot limits in wavelength
   if (color EQ 'blue') then begin
      xrange = [3800, 6400] 
   endif else if (color EQ 'red') then begin
      xrange = [5600, 9400]
   endif else begin
      xrange = [min(lambda), max(lambda)]
      ; Pad wavelength range by an additional 10% on either end
      xrange = [1.1*xrange[0]-0.1*xrange[1], 1.1*xrange[1]-0.1*xrange[0]]
   endelse

   ; Plot panel of where fits to line centers fall in wavelength

   djs_plot, [0], [0], /nodata, xrange=xrange, yrange=[-10,nfiber+10], $
    xstyle=1, ystyle=1, $
    xtitle='\lambda [A] + 500 * Deviation', ytitle='Fiber Number', $
    title='Arcline Fit for  '+trimname
   fibernum = findgen(nfiber)+1
   for k=1, nlamp-1 do $
    djs_oplot, 500*xdif[*,k]+lambda[k], fibernum

   ; Compute offset + stddev for each line center
   meandif = fltarr(nlamp)
   sig = fltarr(nlamp)
   for k=0,nlamp-1 do begin 
      djs_iterstat, xdif[*,k], mean=mn
      meandif[k] = mn
      sig[k] = djsig(xdif[*,k])
   endfor

   ; Make plot of deviations
   ; Multiply by 1000 to convert to milli-pixels

   djs_plot, lambda, meandif*1000, xrange=xrange, yrange=[-100,100], $
    xstyle=1, ystyle=1, $
    xtitle='\lambda [A]', ytitle='Deviation [mPix]'
   errplot, lambda, (meandif-sig)*1000, (meandif+sig)*1000

   xyouts, 0.95, 0., systime(), /normal, align=1 ; , charsize=0.5

   ; End multi-plot format
   !p.multi = pmulti

   return
end

;------------------------------------------------------------------------------
