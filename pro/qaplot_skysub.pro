;+
; NAME:
;   qaplot_skysub
;
; PURPOSE:
;   Generate QA plots for sky-subtraction
;
; CALLING SEQUENCE:
;   qaplot_skysub, obj, objivar, objsub, objsubivar, wset, $
;    iskies, [filename= ]
;
; INPUTS:
;   obj        - Image
;   objivar    - Inverse variance for OBJ
;   objsub     - Image after sky-subtraction
;   objsubivar - Inverse variance for image after sky-subtraction
;   wset       - Wavelength solution
;   iskies     - List of good sky fibers
;
; OPTIONAL KEYWORDS:
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
;   djs_median
;   djs_oplot
;   djs_plot
;
; INTERNAL SUPPORT ROUTINES:
;   skyplot
;
; REVISION HISTORY:
;   30-Dec-1999  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro skyplot, skywave, skyflux, skyfluxsub, skyivar, xrange=xrange, $
 filename=filename

   ; Only plot if the wavelength range spans at least all of XRANGE
   if (10^min(skywave) LT xrange[0] AND 10^max(skywave) GT xrange[1]) then begin

      ; Set multi-plot format
      !p.multi = [0,1,2]

      ii = where(skywave GE alog10(xrange[0]) AND skywave LE alog10(xrange[1]) $
       AND skyivar GT 0)
      ii = ii[ sort(skywave[ii]) ] ; Sort by wavelength
      xaxis = 10^skywave[ii]
      yerr = 1. / sqrt(skyivar[ii])

      djs_plot, xaxis, skyflux[ii], psym=3, xrange=xrange, $
       xtitle='\lambda [A]', ytitle='Flux [electrons]', $
       title = 'Sky Fibers for '+filename
      djs_oploterr, xaxis, skyflux[ii], yerr=yerr
      djs_oplot, xaxis, skyflux[ii]-skyfluxsub[ii], color='red'

      djs_plot, xaxis, skyfluxsub[ii], psym=3, xrange=xrange, $
       xtitle='\lambda [A]', ytitle='Sky-Sub Flux [electrons]', $
       title = 'Sky-subtracted Sky Fibers for '+filename
      djs_oploterr, xaxis, skyfluxsub[ii], yerr=yerr
      djs_oplot, xrange, [0,0], color='red'

      !p.multi = 0
   endif

   return
end

;------------------------------------------------------------------------------

pro qaplot_skysub, obj, objivar, objsub, objsubivar, wset, iskies, $
      filename=filename

   if (NOT keyword_set(filename)) then filename = ''

   dims = size(objsub, /dimens)
   ncol = dims[0]
   nrow = dims[1]

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
   ; Plot median chi^2 and median flux level for sky fibers.

   ; Set multi-plot format
   !p.multi = [0,1,2]

   djs_plot, iskies, $
    djs_median( objsub[*,iskies]^2 * objivar[*,iskies], 1), $
    xrange=[1,nrow], xstyle=1, psym=2, charsize=1.5, $
    xtitle = 'Fiber number', ytitle='Median \chi^2', $
    title='Sky Fibers for '+filename
   djs_oplot, iskies, $
    djs_median( objsub[*,iskies]^2 * objsubivar[*,iskies], 1), $
    psym=2, color='blue'
   xyouts, 10, 0.1, 'BLUE=After rescaling variance', charsize=1.5

   ;----------
   ; Report all sky-subtracted sky fibers with |median flux| > 10 electrons
   ; (Add 1 to fiber number to get it 1-indexed.)

   medskyflux = djs_median(objsub[*,iskies], 1)
   ibad = where(abs(medskyflux) GT 10, nbad)
   if (nbad GT 0) then begin
      for i=0, nbad-1 do $
       splog, 'WARNING: Median sky-sub flux in sky fiber ', iskies[ibad[i]]+1, $
        ' = ', medskyflux[ibad[i]], ' electrons', format='(a,i4,a,f6.1,a)'
   endif

   djs_plot, iskies, medskyflux, $
    xrange=[1,nrow], xstyle=1, psym=2, charsize=1.5, $
    xtitle = 'Fiber number', ytitle='Median Sky-Sub Flux'
   djs_oplot, [1,nrow], [0,0], color='red'

   !p.multi= 0

   ;---------------------------------------------------------------------------
   ; Plot near specific sky features

   skywave = wave[*,iskies]
   skyflux = obj[*,iskies]
   skyfluxsub = objsub[*,iskies]
   skyivar = objsubivar[*,iskies]

   skyplot, skywave, skyflux, skyfluxsub, skyivar, filename=filename, $
    xrange=[5570,5590]

   skyplot, skywave, skyflux, skyfluxsub, skyivar, filename=filename, $
    xrange=[8820,8900]

   ;---------------------------------------------------------------------------
   ; Plot relative chi^2's
   ; Evaluate this for all the sky fibers, but then linearly interpolate
   ; to every Angstrom for plotting.

   ii = where(objsubivar[*,iskies] GT 0) ; Select unmasked pixels
   relwave = 10^(wave[*,iskies])[ii]
   relchi = sqrt( (objivar[*,iskies])[ii] / (objsubivar[*,iskies])[ii] )

   ; Sort by wavelength
   ii = sort(relwave)
   relwave = relwave[ii]
   relchi = relchi[ii]

   ; Report median and worst relative errors
   maxchi = max(relchi, imax)
   splog, 'Median sky-residual chi = ', median(relchi)
   splog, 'Max sky-residual chi = ', maxchi, ' at ', $
    relwave[imax], 'Ang', format='(a,f7.2,a,f7.1,a)'
   if (maxchi GT 4.0) then $
    splog, 'WARNING: Max sky-residual chi = ', maxchi, ' at ', $
     relwave[imax], ' Ang', format='(a,f7.2,a,f7.1,a)'

   xmin = floor( min(relwave) )
   xmax = ceil( max(relwave) )
   dx = 1.0  ; space by 1.0 Angstroms

   ; Set multi-plot format
   npanel = 3
   !p.multi = [0,1,npanel]

   for ipanel=0, npanel-1 do begin

      if (ipanel EQ 0) then $
       title='Rescaled Errors for '+filename $
       else title=''
      xrange = xmin + [ipanel, ipanel+1] * (xmax-xmin) / float(npanel)

      xaxis = xrange[0] + dx * findgen( fix((xrange[1]-xrange[0]) / dx) )
      yaxis = interpol(relchi, relwave, xaxis)

      djs_plot, xaxis, yaxis, $
       xrange=xrange, yrange=[0.0,5.0], xstyle=1, ystyle=1, $
       xtitle='\lambda [A]', ytitle='Relative \chi', $
       title=title, charsize=1.5

      if (ipanel EQ 0) then $
       xyouts, 0.95*xrange[0] + 0.05*xrange[1], 6.0, $
        'Number of sky fibers = ' + strtrim(string(nskies),2), $
        charsize=1.5

   endfor

   !p.multi= 0

   ;---------------------------------------------------------------------------

   return
end
;------------------------------------------------------------------------------
