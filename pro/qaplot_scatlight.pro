;+
; NAME:
;   qaplot_scatlight
;
; PURPOSE:
;   Generate QA plots for scattered light
;
; CALLING SEQUENCE:
;   qaplot_scatlight, scatimage, scatfit, wset=, xcen=, [ fibermask=, $
;    dlambda=dlambda, filename= ]
;
; INPUTS:
;   scatimage  - Extracted scattered light image [NX,NY/8]
;   scatfit    - Scattered light image after fitting [NX,NY]
;
; REQUIRED KEYWORDS:
;   wset       - Wavelength solution
;   xcen       - X positions corresponding to the extracted wavelength sol'n
;                [NY,NTRACE]
;
; OPTIONAL KEYWORDS:
;   fibermask  - Mask of 0 for bad fibers and 1 for good fibers [NFIBER]
;   dlambda    - Bin wavelengths by this number; default to 10 Angstroms
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
;   djs_mean()
;   djs_oploterr
;   djs_plot
;   traceset2xy
;
; REVISION HISTORY:
;   13-Dec-1999  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro qaplot_scatlight, scatimage, scatfit, wset=wset, xcen=xcen, $
 fibermask=fibermask, dlambda=dlambda, filename=filename

   if (NOT keyword_set(dlambda)) then dlambda = 10.0
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber) + 1
   if (NOT keyword_set(filename)) then filename = ''

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

   contour, scatfit, /follow, nlevels = 10, /xstyle, /ystyle, $
    xtitle='X [pix]', ytitle='Y [pix]', $
    title='Scattered light image for '+filename, c_charsize=1.5

   ;---------------------------------------------------------------------------
   ; Plot spectrum of scattered light

   ; Compute the wavelengths for all vectors from the trace set
   traceset2xy, wset, ycen, loglam

   dims = size(loglam, /dimens)
   ny = dims[0]
   nfiber = dims[1]

   dims = size(scatfit, /dimens)
   nx = dims[0]

   xnear = fix(xcen+0.5) ; Nearest pixels to XCEN

   xrange = 10^[min(loglam), max(loglam)]

   ;----------
   ; Divide the plotting range up into bins in wavelength spaced by DLAMBDA

   nbin = fix( (xrange[1]-xrange[0]) / dlambda ) + 1
   meanval = fltarr(nbin)
   minval = fltarr(nbin)
   maxval = fltarr(nbin)
   ioutlier = [0]
   for ibin=0, nbin-1 do begin
      j = where(loglam GE alog10(xrange[0]+dlambda*ibin) $
            AND loglam LT alog10(xrange[0]+dlambda*(ibin+1)) $
            AND xnear GE 0 AND xnear LT nx , ct)
      if (ct GT 0) then begin
         meanval[ibin] = djs_mean( scatfit[xcen[j],ycen[j]] )
         minval[ibin] = min( scatfit[xcen[j],ycen[j]] )
         maxval[ibin] = max( scatfit[xcen[j],ycen[j]] )
      endif

      ; Burles counter of bin number...
      print, format='(i4,a4,i4,a1,$)', ibin, ' of ', nbin, string(13b)
   endfor

   ; Plot - set Y limits according to MAXVAL
   xaxis = xrange[0] + (findgen(nbin)+0.5) * dlambda
   djs_plot, xaxis, maxval, /nodata, xrange=xrange, xstyle=1, $
    xtitle='log \lambda [A]', ytitle='Brightness [electrons/pix]', $
    title='Scattered Light Spectrum for '+filename
   djs_oplot, xaxis, meanval
   djs_oploterr, xaxis, (minval+maxval)/2.0, yerr=(maxval-minval)/2.0

   xpos = xaxis[fix(nbin/20)]
   xyouts, xpos, 0.15*!y.crange[1] + 0.85*!y.crange[0], $
    'Error bars denote full range'
   xyouts, xpos, 0.10*!y.crange[1] + 0.90*!y.crange[0], $
    'Min value =' + string(min(meanval),format='(f9.3)')
   xyouts, xpos, 0.05*!y.crange[1] + 0.95*!y.crange[0], $
    'Max value =' + string(max(meanval),format='(f9.3)')
   xyouts, 0.95, 0., systime(), /normal, align=1 ; , charsize=0.5

   return
end
;------------------------------------------------------------------------------
