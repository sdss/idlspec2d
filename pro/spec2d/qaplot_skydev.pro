;+
; NAME:
;   qaplot_skydev
;
; PURPOSE:
;   Generate QA plot: scatter in sky line positions
;
; CALLING SEQUENCE:
;   qaplot_skydev, flux, fluxivar, vacset, plugsort, color, $
;    [fibermask=, title=]
;
; INPUTS:
;   flux       - Pre-skysubtracted flux array [Nrow, Ntrace]
;   fluxivar   - Associated inverse variance
;   vacset     - Vacuum Wavelength solution
;   plugsort   - Plugmap struct [Ntrace]
;   color      - string specifying 'red' or 'blue' spectra
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   title      - TITLE of plot
;
; OUTPUTS:
;   Output plots only
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This code has been deprecated in favor of QAPLOT_SKYSHIFT.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_oplot
;   djs_plot
;   trace_gweight
;   traceset2xy
;   traceset2pix
;   djs_median
;
; REVISION HISTORY:
;   22-Jan-2000 Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------

pro qaplot_skydev, flux, fluxivar, vacset, plugsort, color, $
 fibermask=fibermask, title=title

  if (N_params() LT 5) then begin
      print, 'Syntax - qaplot_skydev, flux, fluxivar, vacset, plugsort,$'
      print, '         color, [fibermask=, title= ]'
      return
   endif

        nTrace = (size(flux))[2]
;
;	Here are guesses at strong skylines in red and blue
;       (corrected to vacuum wavelengths)
;
        if (color EQ 'red') then $
	  possible = [5891.8, 6866.2, 7278.6, 7343.2, $
                    7753.0, 7796.5, 7824.0, 7916.0, 7967.1, $
                    7995.6, 8028.2, 8401.6, 8432.5, 8829.6, $
                    8888.2, 8922.0, 8960.6] $
        else if (color EQ 'blue') then $
	  possible = [4047.706,  4359.553, 5462.252, 5578.894, 5891.583]

	logpos = alog10(possible)

	pix = traceset2pix(vacset, logpos)
        ycen = long(pix)
	for i=0,nTrace - 1 do ycen[*,i] = i
	
        skycen = trace_gweight(flux, pix, ycen, $
                   invvar=fluxivar, xerr=skyerr, sigma=1.0)
        skycen = trace_gweight(flux, skycen, ycen, $
                   invvar=fluxivar, xerr=skyerr, sigma=1.0)
        skycen = trace_gweight(flux, skycen, ycen, $
                   invvar=fluxivar, xerr=skyerr, sigma=1.0)

	traceset2xy, vacset, skycen, skywave

	skymed = djs_median(skywave,2)

        veldiff = (skywave-skymed # (fltarr(nTrace) +1.0)) * 3.0e5 * 2.302585

	veldiff = veldiff * (skyerr LT 990)

        djs_plot, possible # (fltarr(nTrace) +1.0), veldiff, ps=3, $
             yr=[-40,40], xtitle='Wavelength of sky lines', $
             ytitle='Offsets from median (km/s)', $
             title=title+' Skyline Positions'

        skies = where(plugsort.objtype EQ 'SKY',nskies)
	if (skies[0] NE -1) then $
          djs_oplot, possible # (fltarr(nskies) +1.0), veldiff[*,skies], $
               ps=4, color='red'

	return
end
