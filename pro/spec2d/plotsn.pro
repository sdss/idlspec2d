;+
; NAME:
;   plotsn
;
; PURPOSE:
;   Plot S/N and residuals in up to 3 bands
;
; CALLING SEQUENCE:
;   plotsn, snvec, plugmap, [ bands=, plotmag=, fitmag=, plottitle=, plotfile= ]
;
; INPUTS:
;   snvec      - S/N array [nbands, nfibers]
;   plugmap    - Plugmap structure [nfibers]
;
; OPTIONAL KEYWORDS:
;   bands      - Index of bands to fit; default to [1,2,3] for g,r,i bands.
;   plotmag    - Magnitude range for plotting; default to [15,21].
;   fitmag     - Magnitude range over which to fit (S/N) as function of mag;
;                default to [-1,0.5] mags around the fiducial magntidues
;                at which we evaluate the fit.
;                default to the same as PLOTMAG.
;   plottitle  - Title for top of plot
;   plotfile   - Optional plot file
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The fiducial magnitudes at which to evaluate the fit to (S/N) are
;   [19.0, 19.2, 19.25, 18.9, 18.0] in u,g,r,i,z bands.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_arrow
;   djs_icolor()
;   djs_oplot
;   djs_plot
;
; REVISION HISTORY:
;   15-Apr-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
pro plotsn, snvec, plug, bands=bands, plotmag=plotmag, fitmag=fitmag, $
 plottitle=plottitle, plotfile=plotfile

   if (size(snvec,/n_dim) NE 2) then return
   if (NOT keyword_set(bands)) then bands=[1, 2, 3]
   if (NOT keyword_set(plotmag)) then plotmag = [15.0, 21.0]

   ; Set fiducial magnitudes about which to fit and evaluate the fit
   snmag = [19.0, 19.2, 19.25, 18.9, 18.0]

   nbands = n_elements(bands)
   nfibers = n_elements(plug)

   if ((size(snvec,/dimen))[0] NE nbands) then return
   if ((size(snvec,/dimen))[1] NE nfibers) then return

   bandnames = ["u'","g'","r'","i'","z'"]
   slopelabel = " * "+bandnames
   snlabel = '(S/N)^2 @ '+bandnames+' ='+string(snmag,format='(f6.2)')

   iobj = where(strtrim(plug.objtype,2) NE 'SKY', nobj)
   if (nobj LT 3) then return

   ; The following variables are defined only for non-SKY objects...
   plugc = plug[iobj]
   spectroid = (plugc.fiberid - 1) / 320 ; Set to 0 or 1
   s1 = where(spectroid EQ 0)
   s2 = where(spectroid EQ 1)

   ;----------
   ; Open plot file

   if (keyword_set(plotfile)) then begin
      if (nbands LE 2) then ysize=7.0 $
       else ysize=10.0
      set_plot, 'ps'
      device, filename=plotfile, /color, /portrait, $
       xsize=8.0, ysize=ysize, xoffset=0.25, yoffset=0.5, /inch
   endif

   oldmulti = !p.multi  
   !p.multi = [0,2,nbands]

   fiducialfits = [[6.12, -0.28], $    ; u' fit
                   [6.12, -0.28], $    ; g' fit
                   [6.18, -0.28], $    ; r' fit
                   [6.05, -0.28], $    ; i' fit
                   [6.05, -0.28]]      ; z' fit

   ;---------------------------------------------------------------------------
   ; Loop over each band in the plot
   ;---------------------------------------------------------------------------

   for iband=0, nbands-1 do begin

      ;------------------------------------------------------------------------
      ; PLOT 1: (S/N) vs. magnitude
      ;------------------------------------------------------------------------

      ;----------
      ; Select the data points: S/N and magnitude for all objects that
      ; are not SKY fibers.

      snc = snvec[iband,iobj]
      mag = plugc.mag[bands[iband]]

      ;----------
      ; Set up the plot

      if (iband LT nbands-1) then begin
         xchars = 0.001
         ymargin = [0,2]
         xtitle = 'mag'
      endif else begin
         xchars = 1.0
         ymargin = [3,1]
         xtitle = ''
      endelse
      xmargin = [8,1]
      ychars = 1.0
      symsize = 0.65

      plot, mag, snc, /nodata, /ylog, $
       xchars=xchars, ychars=ychars, xrange=plotmag, $
       xtitle=xtitle, ytitle='S/N in '+bandnames[bands[iband]]+'-band', $
       xmargin=xmargin, ymargin=ymargin, /xstyle, yrange=[0.5,100], /ystyle

      if (iband EQ 0 AND keyword_set(plottitle)) then $
       xyouts, plotmag[1], 110.0, plottitle, align=0.5

      ;----------
      ; Plot the fiducial line (a thick blue line)

      djs_oplot, plotmag, 10^poly(plotmag, fiducialfits[*,bands[iband]]), $
       color='blue', thick=5

      ;----------
      ; Plot a fit to the data done across the entire plotting range.
      ; Also, identify which points fall above and below this fit.

      afit = fitsn(mag, snc, sigma=sigma, fitmag=plotmag)
      logsnc = alog10(snc > 0.01)
      diff = logsnc - poly(mag, afit) ; Residuals from this fit

      ;----------
      ; Plot the data points, (S/N) vs. magnitude.
      ; Color code green for positive residuals, red for negative.
      ; Use different symbols for each spectrograph.
      ; Plot these first so that the lines are visibly plotted on top.

      psymvec = (spectroid EQ 0) * 7 + (spectroid EQ 1) * 6 
      colorvec = replicate('red', nobj)
      ipos = where(diff GE 0)
      if (ipos[0] NE -1) then colorvec[ipos] = 'green'
      djs_oplot, mag, snc > 0.6, psym=psymvec, symsize=symsize, $
       color=colorvec

      ; Now overplot the fit line
      if (keyword_set(afit)) then $
       djs_oplot, plotmag, 10^poly(plotmag, afit)

      ;----------
      ; Plot a fit to the data in the range specified by FITMAG,
      ; independently for each spectrograph.
      ; Also, draw an arrow that terminates at the S/N at the magnitude
      ; where we measure the canonical (S/N)^2.

      if (keyword_set(fitmag)) then $
       myfitmag = fitmag $
      else $
       myfitmag = snmag[bands[iband]] + [-1.0,0.5]

      snoise2 = fltarr(2)
      xloc = snmag[bands[iband]]
      if (s1[0] NE -1) then begin
         afit1 = fitsn(mag[s1], snc[s1], fitmag=myfitmag)
         if (keyword_set(afit1)) then begin
            snoise2[0] = 10^(2.0 * poly(snmag[bands[iband]],afit1)) 
            yloc = 10^poly(xloc,afit1)
            djs_arrow, xloc, 20, xloc, yloc, /data
         endif
      endif
      if (s2[0] NE -1) then begin
         afit2 = fitsn(mag[s2], snc[s2], fitmag=myfitmag)
         if (keyword_set(afit2)) then begin
            snoise2[1] = 10^(2.0 * poly(snmag[bands[iband]],afit2)) 
            yloc = 10^poly(xloc,afit2)
            djs_arrow, xloc, 20, xloc, yloc, /data
         endif
      endif

      ;----------
      ; Label the plot

      if (keyword_set(afit)) then $
       xyouts, plotmag[0]+0.5, 0.9, string(format='(a,f5.2,f6.2,a)', $
        'log S/N = ', afit, slopelabel[bands[iband]])

      if (keyword_set(sigma)) then $
       xyouts, plotmag[0]+0.5, 1.5, string(format='(a,f4.2)','Stdev=', sigma)

      djs_xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.3], $
       10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.93], $
       snlabel[bands[iband]]

      djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.50], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.85], psym=7, $
         symsize=symsize
      djs_xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.53], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.83], $
         string(format='("Spec1: ", f5.1)', snoise2[0])

      djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.50], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.75], psym=6, $
         symsize=symsize
      djs_xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.53], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.73], $
         string(format='("Spec2: ", f5.1)', snoise2[1])

      splog, snlabel[bands[iband]], snoise2

      ;------------------------------------------------------------------------
      ; PLOT 2: Throughput deviations plotted on the focal plane
      ;------------------------------------------------------------------------

      good = where(snc GT 0.0 AND mag GE plotmag[0] $
       AND mag LE plotmag[1], ngood)

      plot, [0], [0], /nodata, xchars=xchars, ychars=ychars, $
       xtitle='X [mm]', ytitle='Y [mm]', $
       xrange=[-320,320], yrange=[-320,320], xstyle=1, ystyle=1, $
       xmargin=xmargin, ymargin=ymargin
      if (ngood GT 0) then begin
         colorvec = (diff[good] GE 0) * djs_icolor('green') $
          + (diff[good] LT 0) * djs_icolor('red')
         symvec = (abs(diff[good]) < 1.0) * 3
         djs_oplot, plugc[good].xfocal, plugc[good].yfocal, $
          symsize=symvec, color=colorvec, psym=2
      endif

   endfor

   !p.multi  = oldmulti

   ;----------
   ; Close plot file

   if (keyword_set(plotfile)) then begin
      device, /close
      set_plot, 'x'
   endif

   return
end
