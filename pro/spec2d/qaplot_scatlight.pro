;+
; NAME:
;   qaplot_scatlight
;
; PURPOSE:
;   Generate QA plots for scattered light
;
; CALLING SEQUENCE:
;   qaplot_scatlight, scatfit, yrow, wset=, xcen=, [ fibermask=, title= ]
;
; INPUTS:
;   scatfit    - Scattered light image after fitting [NX,NY/8]
;   yrow       - extracted rows in scatfit 
;
; REQUIRED KEYWORDS:
;   wset       - Wavelength solution
;   xcen       - X positions corresponding to the extracted wavelength sol'n
;                [NY,NTRACE]
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
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
;   djs_median()
;   djs_oploterr
;   djs_plot
;   traceset2xy
;
; REVISION HISTORY:
;   13-Dec-1999  Written by D. Schlegel, Princeton
;   22-Jan-2000  Trying to speed up, S. Burles
;-
;------------------------------------------------------------------------------

pro qaplot_scatlight, scatfit, yrow, wset=wset, xcen=xcen, $
 fibermask=fibermask, title=title

   if (N_params() LT 2) then begin
      print, 'Syntax - qaplot_scatlight, scatfit, yrow, [ wset=, xcen=, $'
      print, '         fibermask=, title= ]'
      return
   endif

   if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber) 
   if (NOT keyword_set(title)) then title = ''

   scatmed = median(scatfit)
   scatmax = max(scatfit)
   if (scatmed GT 20) then $
     splog, 'WARNING: Scattered light median = ', scatmed, ' electrons' $
    else $
     splog, 'Scattered light median = ', scatmed, ' electrons'
   if (scatmax GT 40) then $
     splog, 'WARNING: Scattered light max = ', scatmax, ' electrons' $
    else $
     splog, 'Scattered light max = ', scatmax, ' electrons'

   ;---------------------------------------------------------------------------
   ; Plot contour image of scattered light

   contour, (scatfit[*,yrow])[yrow,*], yrow, yrow, $
    /follow, nlevels = 10, /xstyle, /ystyle, $
    xtitle='X [pix]', ytitle='Y [pix]', $
    title=title, c_charsize=1.5

   ;---------------------------------------------------------------------------
   ; Plot spectrum of scattered light

   ; Compute the wavelengths for all vectors from the trace set
   traceset2xy, wset, ycen, loglam

   ;$  "Extract" Scattered Light, and divide by 2.0 * radius

   lamsort = sort(loglam) 
   dims = size(loglam, /dimens)
   ny = dims[0]
   nfiber = dims[1]
  
   xsafe = (xcen > 1) <  ny - 2
   scatextract = (extract_boxcar(scatfit, xsafe, radius=1.0)) / 2.0

   midindex = lindgen(nfiber)*ny + ny/2 
   binwave = 10^loglam[lamsort[midindex]] 

   scattemp = reform(scatextract[lamsort],ny,nfiber)
   meanval = total(scattemp,1)/ny
   medianval = djs_median(scattemp,1)
   
   minval = fltarr(nfiber)
   maxval = fltarr(nfiber)

   for ibin=0,nfiber - 1 do minval[ibin] = min(scattemp[*,ibin])
   for ibin=0,nfiber - 1 do maxval[ibin] = max(scattemp[*,ibin])

   xrange = 10^[min(loglam), max(loglam)] + [-100,100] ; Pad by 100 Angstroms

   ; Plot - set Y limits according to MAXVAL
   djs_plot, binwave, maxval, /nodata, xrange=xrange, xstyle=1, $
    xtitle='\lambda [A]', ytitle='Brightness [electrons/pix]', $
    title=title
   djs_oploterr, binwave, (minval+maxval)/2.0, yerr=(maxval-minval)/2.0
   djs_oplot, binwave, medianval
   djs_oplot, binwave, meanval, color='red'

   xpos = binwave[fix(nfiber/20)]
   xyouts, xpos, 0.15*!y.crange[1] + 0.85*!y.crange[0], $
    'Error bars denote full range'
   xyouts, xpos, 0.10*!y.crange[1] + 0.90*!y.crange[0], $
    'Min value =' + string(min(medianval),format='(f9.3)')
   xyouts, xpos, 0.05*!y.crange[1] + 0.95*!y.crange[0], $
    'Max value =' + string(max(medianval),format='(f9.3)')
   xyouts, 0.95, 0., systime(), /normal, align=1 ; , charsize=0.5

   return
end
;------------------------------------------------------------------------------
