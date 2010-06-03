;+
; NAME:
;   qaplot_skysub
;
; PURPOSE:
;   Generate QA plots for sky-subtraction
;
; CALLING SEQUENCE:
;   qaplot_skysub, obj, objivar, objsub, objsubivar, wset, $
;    iskies, [title= ]
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
;   title      - TITLE of plot
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
;   djs_oploterr
;   djs_plot
;   splog
;   traceset2xy
;
; INTERNAL SUPPORT ROUTINES:
;   skyplot
;
; REVISION HISTORY:
;   30-Dec-1999  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro skyplot, skywave, skyflux, skyfluxsub, skyivar, xrange=xrange, $
 title=title

   ; Only plot if the wavelength range spans at least all of XRANGE
   jj = where(skyivar GT 0)
   if (10^min(skywave[jj]) LT xrange[0] AND 10^max(skywave[jj]) GT xrange[1]) then begin

      ; Set multi-plot format
      !p.multi = [0,1,2]

      ii = where(skywave GE alog10(xrange[0]) AND skywave LE alog10(xrange[1]) $
       AND skyivar GT 0)
      ii = ii[ sort(skywave[ii]) ] ; Sort by wavelength
      xaxis = 10^skywave[ii]
      yerr = 1. / sqrt(skyivar[ii])

      djs_plot, xaxis, skyflux[ii], psym=3, xrange=xrange, $
       xtitle='\lambda [A]', ytitle='Flux [electrons]', $
       title=title+' Sky Fibers'
      djs_oploterr, xaxis, skyflux[ii], yerr=yerr
      djs_oplot, xaxis, skyflux[ii]-skyfluxsub[ii], color='red'

      djs_plot, xaxis, skyfluxsub[ii], psym=3, xrange=xrange, $
       xtitle='\lambda [A]', ytitle='Sky-Sub Flux [electrons]', $
       title=title+' Sky-subtracted Sky Fibers'
      djs_oploterr, xaxis, skyfluxsub[ii], yerr=yerr
      djs_oplot, xrange, [0,0], color='red'

      !p.multi = 0
   endif

   return
end

;------------------------------------------------------------------------------
pro qaplot_skysub, obj, objivar, objsub, objsubivar, wset, iskies, $
 title=title

   if (NOT keyword_set(title)) then title = ''

   dims = size(objsub, /dimens)
   ncol = dims[0]
   nrow = dims[1]

   ;----------
   ; Solve for wavelength of each pixel 

   traceset2xy, wset, pixnorm, wave

   ;----------
   ; Find sky fibers

   nsky = n_elements(iskies)

   if (nsky EQ 0) then begin
      splog, 'No sky fibers!'
      return
   endif

   ;----------
   ; Compute median sky flux and abs(chi)
   medskyflux = fltarr(nsky)
   medskychi = fltarr(nsky)
   medskysubchi = fltarr(nsky)
   for j=0, nsky-1 do begin
      qgood = where(objsubivar[*, iskies[j]] GT 0, ct)
      if (ct GT 0) then begin
         medskyflux[j] = djs_median(objsub[qgood,iskies[j]])
         medskychi[j] = djs_median(abs(objsub[qgood,iskies[j]]) $
          * sqrt(objivar[qgood,iskies[j]]))
         medskysubchi[j] = djs_median(abs(objsub[qgood,iskies[j]]) $
          * sqrt(objsubivar[qgood,iskies[j]]))
      endif
   endfor

   ;---------------------------------------------------------------------------
   ; Plot median chi^2 and median flux level for sky fibers.

   ; Set multi-plot format
   !p.multi = [0,1,2]

   csize = 0.8
   xrange = [0,nrow+1]
   yrange = [0, max(medskychi) > 1]
   djs_plot, [iskies], [medskychi], $
    xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
    psym=2, charsize=csize, $
    xtitle = 'Fiber number', ytitle='Median |\chi|', $
    title=title+' Sky Fibers'
   djs_oplot, [iskies], [medskysubchi], $
    psym=2, color='blue'
   djs_oplot, xrange, [0.67,0.67] ; expected level if Gaussian
   xyouts, 0.94*xrange[1], 0.72, 'Gaussian', align=1.0
   xyouts, 0.06*xrange[1], 0.90*yrange[1], $
    'BLACK=Before rescaling variance', charsize=csize
   xyouts, 0.06*xrange[1], 0.82*yrange[1], $
    'BLUE=After rescaling variance', charsize=csize

   ;----------
   ; Report all sky-subtracted sky fibers with |median flux| > 10 electrons
   ; (Add 1 to fiber number to get it 1-indexed.)

   ibad = where(abs(medskyflux) GT 10, nbad)
   if (nbad GT 0) then begin
      for i=0, nbad-1 do $
       splog, 'Warning: Median sky-sub flux in sky fiber ', iskies[ibad[i]]+1, $
        ' = ', medskyflux[ibad[i]], ' electrons'
   endif

   djs_plot, [iskies], [medskyflux], $
    xrange=[1,nrow], xstyle=1, psym=2, charsize=csize, $
    xtitle = 'Fiber number', ytitle='Median Sky-Sub Flux'
   djs_oplot, [1,nrow], [0,0], color='red'

   !p.multi= 0

   ;---------------------------------------------------------------------------
   ; Plot near specific sky features

   skywave = wave[*,iskies]
   skyflux = obj[*,iskies]
   skyfluxsub = objsub[*,iskies]
   skyivar = objsubivar[*,iskies]
   skyplot, skywave, skyflux, skyfluxsub, skyivar, title=title, $
    xrange=[5570,5590]

   skyplot, skywave, skyflux, skyfluxsub, skyivar, title=title, $
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
    relwave[imax], ' Ang'
   if (maxchi GT 5.0) then $
    splog, 'WARNING: Max sky-residual chi = ', maxchi, ' at ', $
     relwave[imax], ' Ang'

   xmin = floor( min(relwave) )
   xmax = ceil( max(relwave) )
   dx = 1.0  ; space by 1.0 Angstroms

   ; Set multi-plot format
   npanel = 3
   !p.multi = [0,1,npanel]

   for ipanel=0, npanel-1 do begin

      if (ipanel EQ 0) then $
       title=title+' Rescaled Errors' $
      else title=''
      xrange = xmin + [ipanel, ipanel+1] * (xmax-xmin) / float(npanel)

      xaxis = xrange[0] + dx * findgen( fix((xrange[1]-xrange[0]) / dx) )
      yaxis = interpol(relchi, relwave, xaxis)

      djs_plot, xaxis, yaxis, $
       xrange=xrange, yrange=[0.0,5.0], xstyle=1, ystyle=1, $
       xtitle='\lambda [A]', ytitle='Relative \chi', $
       title=title, charsize=csize

      if (ipanel EQ 0) then $
       xyouts, 0.95*xrange[0] + 0.05*xrange[1], 6.0, $
        'Number of sky fibers = ' + strtrim(string(nsky),2), $
        charsize=csize

   endfor
   !p.multi= 0

   ;---------------------------------------------------------------------------

   return
end
;------------------------------------------------------------------------------
