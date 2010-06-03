;+
; NAME:
;   qaplot_arcline
;
; PURPOSE:
;   Generate QA plot for arcline fits
;
; CALLING SEQUENCE:
;   qaplot_arcline, xdif, lambda, [rejline=, color=, fibermask=, title=]
;
; INPUTS:
;   xdif       - Deviations in arc line fits in pixels, array [NFIBER, NLAMP]
;   lambda     - Wavelength for each lamp in Angstroms, vector [NLAMP]
;   rejline    - String set to nonzero for any rejected arc line (plot in red)
;
; OPTIONAL KEYWORDS:
;   color      - Specify 'blue' or 'red' to fix plotting limits
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   title      - TITLE of plot
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Skip over bad fibers in plots.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_iterstat
;   djs_oplot
;   djs_plot
;   errplot
;   fibermask_bits()
;   traceset2xy
;   splog
;
; REVISION HISTORY:
;   15-Oct-1999  Written by D. Finkbeiner, APO
;   23-Nov-1999  Modified by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro qaplot_arcline, xdif, wset, lambda, rejline=rejline, $
 color=color, fibermask=fibermask, title=title

   if (NOT keyword_set(title)) then title = ''

   dims = size(xdif, /dimens)
   nfiber = dims[0]
   nlamp = dims[1]

   if (keyword_set(fibermask)) then $
    igfib = where((fibermask AND fibermask_bits('BADARC')) EQ 0) $
   else $
    igfib = lindgen(nfiber)

   ;----------
   ; Compute offset + stddev for each line center

   mnarr = fltarr(nlamp)
   sgarr = fltarr(nlamp)
   for k=0, nlamp-1 do begin 
      djs_iterstat, xdif[igfib,k], mean=mn, sigma=sg
      mnarr[k] = mn
      sgarr[k] = sg
   endfor

   ;----------
   ; Determine the starting and ending wavelengths in each fiber
   ; (Note that MINWAVE is larger than MAXWAVE if the wavelengths
   ; run in reverse direction.)

   traceset2xy, wset, transpose(replicate(wset.xmin,nfiber)), minwave
   traceset2xy, wset, transpose(replicate(wset.xmax,nfiber)), maxwave

   ;---------------------------------------------------------------------------
   ; Print residuals in arc fits

   splog, 'Number of arc lines: ', nlamp

   if (keyword_set(color)) then begin
      if (color EQ 'blue') then begin
         junk = where(lambda LT 4000., ct)
         splog, 'Number of arc lines below 4000 Ang: ', ct
      endif else if (color EQ 'red') then begin
         junk = where(lambda GT 8000., ct)
         splog, 'Number of arc lines above 8000 Ang: ', ct
      endif
   endif

   for k=0, nlamp-1 do $
    splog, 'Arcline ',k,': lambda=',lambda[k], $
     ' Ang, median dev=', mnarr[k], $
     ' pix, sigma dev=', sgarr[k], ' pix ', rejline[k], $
     format='(A,I3,A,F8.2,A,F7.3,A,F7.3,A,A)'

   ;---------------------------------------------------------------------------
   ; Set multi-plot format
   !p.multi = [0,1,2]

   igood = where(rejline EQ '', ngood)
   ibad = where(rejline NE '', nbad)

   ; Determine the plot limits in wavelength
   ; Pad wavelength range by an additional 5% on either end

   xrange = 10^[min([minwave,maxwave]), max([minwave,maxwave])]
   xrange = [1.05*xrange[0]-0.05*xrange[1], 1.05*xrange[1]-0.05*xrange[0]]

   ;----------
   ; Plot panel of where fits to line centers fall in wavelength

   djs_plot, [0], [0], /nodata, xrange=xrange, yrange=[-10,nfiber+10], $
    xstyle=1, ystyle=1, $
    xtitle='\lambda [A] + 300 * Deviation', ytitle='Fiber Number', $
    title=title
   fibernum = findgen(nfiber)+1
   if (nbad GT 0) then $
    for k=0, nbad-1 do $
     djs_oplot, 300*xdif[igfib,ibad[k]]+lambda[ibad[k]], fibernum[igfib], color='red'
   if (ngood GT 0) then $
    for k=0, ngood-1 do $
     djs_oplot, 300*xdif[igfib,igood[k]]+lambda[igood[k]], fibernum[igfib]

   djs_oplot, 10^minwave, fibernum, color='green', thick=2
   djs_oplot, 10^maxwave, fibernum, color='green', thick=2

   ;----------
   ; Make plot of deviations

   yrange = [-0.15,0.15]
   djs_plot, lambda[igood], mnarr[igood], xrange=xrange, yrange=yrange, $
    psym=6, xstyle=1, ystyle=1, $
    xtitle='\lambda [A]', ytitle='Deviation [pix]'
   if (nbad GT 0) then $
    djs_oplot, lambda[ibad], mnarr[ibad], psym=6, color='red'
   djs_oplot, xrange, [0,0]
   ylo = (mnarr - sgarr) < 0.9*yrange[1]
   yhi = (mnarr + sgarr) > 0.9*yrange[0]
   if (nbad GT 0) then $
    errplot, lambda[ibad], ylo[ibad], yhi[ibad], $
     color=djs_icolor('red')
   if (ngood GT 0) then $
    errplot, lambda[igood], ylo[igood], yhi[igood]
   djs_xyouts, total([0.95,0.05]*!x.crange), total([0.10,0.90]*!y.crange), $
    'Ncoeff = ' + strtrim(n_elements(wset.coeff[*,0]),2)

   xyouts, 0.95, 0., systime(), /normal, align=1 ; , charsize=0.5

   ; End multi-plot format
   !p.multi = 0

   return
end
;------------------------------------------------------------------------------
