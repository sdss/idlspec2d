;+
; NAME:
;   plotsn
;
; PURPOSE:
;   Plot S/N and residuals in up to 3 bands
;
; CALLING SEQUENCE:
;   plotsn, snvec, plugmap, [ bands=, plotmag=, fitmag=, plottitle=, $
;    plotfile=, synthmag=, snplate= ]
;
; INPUTS:
;   snvec      - S/N array [nbands, nfibers]
;   plugmap    - Plugmap structure [nfibers]
;
; OPTIONAL KEYWORDS:
;   bands      - Index of bands to fit; default to [1,2,3] for g,r,i bands.
;   fitmag     - Magnitude range over which to fit (S/N) as function of mag;
;                default to those used by FITSN().
;   plotmag    - Magnitude range for plotting; default to [15,21], but
;                extend the range to include FITMAG if necessary.
;   plottitle  - Title for top of plot
;   plotfile   - Optional plot file
;   synthmag   - Vector of synthetic magnitudes dimensionsed [5,640],
;                but only the central three mags (gri) are used;
;                if set, then make more than the first page of plots
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   snplate    - Best fit (S/N)^2 at fiducial magnitude
;
; COMMENTS:
;   The fiducial magnitudes at which to evaluate the fit to (S/N) are
;   [20.0, 20.2, 20.25, 19.9, 19.0] in u,g,r,i,z bands.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_arrow
;   djs_icolor()
;   djs_iterstat
;   djs_oplot
;   djs_plot
;   djs_xyouts
;   fitsn()
;   splog
;
; REVISION HISTORY:
;   28-Jan-2003  Add E(B-V) keyword to allow for spectra which have 
;                foreground reddening removed
;   20-Oct-2002  Plot range on synthmag vs fiber mag plot changed by C. Tremonti
;   15-Apr-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
pro plotsn1, plugc, synthmag, i1, i2, plottitle=plottitle, objtype=objtype

   pmulti = !p.multi
   ymargin = !y.margin
   yomargin = !y.omargin

   !p.multi = [0,2,5]
   !y.margin = [1,0]
   !y.omargin = [5,3]
   symsize = 0.5
   psym = 1
   csize = 1.5
   textsize = 1.0
   xrange = [14.,24.]
   yrange = [-0.6,0.6]

   magdiff = synthmag - plugc.mag

   for ipanel=0, 4 do begin
      case ipanel of
      0: begin
         yplot = magdiff[1,*]
         ytext = 'g-mag'
         end
      1: begin
         yplot = magdiff[2,*]
         ytext = 'r-mag'
         end
      2: begin
         yplot = magdiff[3,*]
         ytext = 'i-mag'
         end
      3: begin
         yplot = magdiff[1,*] - magdiff[2,*]
         ytext = '(g-r) color'
         end
      4: begin
         yplot = magdiff[2,*] - magdiff[3,*]
         ytext = '(r-i) color'
         end
      endcase

      if (ipanel EQ 0) then thistitle = plottitle+'  '+objtype $
       else thistitle = ''
      if (ipanel EQ 4) then xtitle = 'PHOTO Magnitude (r)' $
       else xtitle = ''
      if (ipanel EQ 4) then xtickname = '' $
       else xtickname = strarr(20)+' '

      plot, xrange, [0,0], $
       xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
       xtitle=xtitle, ytitle='(Spectro - PHOTO) mag', $
       title=thistitle, charsize=csize, xtickname=xtickname
      xyouts, xrange[0]+1, 0.75*yrange[1], 'Spectro-1: '+ytext, $
       charsize=textsize
      if (i1[0] NE -1) then begin
         djs_oplot, plugc[i1].mag[2], yplot[i1], psym=psym, symsize=symsize
         djs_iterstat, yplot[i1], median=mn, sigma=sig
         mntext = string(mn, format='(" Median= ",f6.3)')
         devtext = string(sig, format='(" Stdev= ",f6.3)')
         xyouts, 0.5*xrange[0]+0.5*xrange[1], 0.40*yrange[0], charsize=textsize, $
          mntext
         xyouts, 0.5*xrange[0]+0.5*xrange[1], 0.70*yrange[0], charsize=textsize, $
          devtext
         splog, 'Spectro-1: ' + objtype + ' ' + ytext + mntext + devtext
      endif

      plot, xrange, [0,0], $
       xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
       xtitle=xtitle, ytitle='(Spectro - PHOTO) mag', $
       title=thistitle, charsize=csize, xtickname=xtickname
      xyouts, xrange[0]+1, 0.75*yrange[1], 'Spectro-2: '+ytext, $
       charsize=textsize
      if (i2[0] NE -1) then begin
         djs_oplot, plugc[i2].mag[2], yplot[i2], psym=psym, symsize=symsize
         djs_iterstat, yplot[i2], median=mn, sigma=sig
         mntext = string(mn, format='(" Median= ",f6.3)')
         devtext = string(sig, format='(" Stdev= ",f6.3)')
         xyouts, 0.5*xrange[0]+0.5*xrange[1], 0.40*yrange[0], charsize=textsize, $
          mntext
         xyouts, 0.5*xrange[0]+0.5*xrange[1], 0.70*yrange[0], charsize=textsize, $
          devtext
         splog, 'Spectro-2: ' + objtype + ' ' + ytext + mntext + devtext
      endif
   endfor

   !p.multi = pmulti
   !y.margin = ymargin
   !y.omargin = yomargin

   return
end
;------------------------------------------------------------------------------
pro plotsn, snvec, plug, bands=bands, plotmag=plotmag, fitmag=fitmag, $
 plottitle=plottitle, plotfile=plotfile, synthmag=synthmag, snplate=snplate

   if (size(snvec,/n_dim) NE 2) then return
   if (NOT keyword_set(bands)) then bands=[1, 2, 3]
   if (NOT keyword_set(plotmag)) then plotmag = [15.0, 21.0]

   ; Set fiducial magnitudes about which to fit and evaluate the fit
   snmag = [20.0, 20.2, 20.25, 19.9, 19.0]

   nbands = n_elements(bands)
   nfibers = n_elements(plug)

   if ((size(snvec,/dimen))[0] NE nbands) then return
   if ((size(snvec,/dimen))[1] NE nfibers) then return

   bandnames = ['u','g','r','i','z']
   slopelabel = " * "+bandnames
   snlabel = '(S/N)^2 @ '+bandnames+' ='+string(snmag,format='(f7.2)')

   ; Use the CALIBFLUX magnitudes instead of MAG, if they exist
   ; from the call to READPLUGMAP() that generated this structure.
   ; Any very small fluxes have MAG=0 so that they are ignored in the plots.
   plugc = plug
   if (tag_exist(plugc,'CALIBFLUX')) then begin
      minflux = 0.1
      plugc.mag = (22.5 - 2.5*alog10(plugc.calibflux > minflux)) $
       * (plugc.calibflux GT minflux)
   endif

   nobj = n_elements(plugc)
   qgood = strtrim(plugc.objtype,2) NE 'SKY' AND plugc.mag[2] GT 0
   igood = where(qgood, ngood)
   if (ngood LT 3) then begin
      splog, 'Warning: Too few non-sky objects to plot'
      return
   endif

   ; The following variables are defined only for non-SKY objects...
   s1 = where(plugc.spectrographid EQ 1 AND qgood)
   s2 = where(plugc.spectrographid EQ 2 AND qgood)

   ;----------
   ; Open plot file

   if (keyword_set(plotfile)) then begin
      if (nbands LE 2) then ysize=7.0 $
       else ysize=9.5
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

   snplate = fltarr(2,nbands)

   pmulti = !p.multi
   ymargin = !y.margin
   yomargin = !y.omargin

   !p.multi = [0,2,nbands]
   !y.margin = [1,0]
   !y.omargin = [5,3]
   csize = 1.5
   textsize = 1.0

   ;---------------------------------------------------------------------------
   ; Loop over each band in the plot
   ;---------------------------------------------------------------------------

   for iband=0, nbands-1 do begin

      ;------------------------------------------------------------------------
      ; 1st PAGE PLOT 1: (S/N) vs. magnitude
      ;------------------------------------------------------------------------

      thismag = plugc.mag[bands[iband]]

      ;----------
      ; Fit the data as S/N vs. mag

      afit = fitsn(thismag[igood], snvec[iband,igood], sigma=sigma, $
       colorband=bandnames[bands[iband]], fitmag=fitmag)
      logsnc = alog10(snvec[iband,*] > 0.01)
      sndiff = logsnc - poly(thismag, afit) ; Residuals from this fit

      ;----------
      ; Extend the plotting range if necessary to include FITMAG

      plotmag[0] = plotmag[0] < fitmag[0]
      plotmag[1] = plotmag[1] > fitmag[1]

      ;----------
      ; Set up the plot

      if (iband LT nbands-1) then begin
         xtickname = strarr(20)+' '
         xtitle1 = ''
         xtitle2 = ''
      endif else begin
         xtickname = ''
         xtitle1 = 'mag'
         xtitle2 = 'X [mm]'
      endelse
      symsize = 0.65

      plot, thismag[igood], snvec[iband,igood], /nodata, /ylog, $
       xtickname=xtickname, xrange=plotmag, $
       xtitle=xtitle1, ytitle='S/N in '+bandnames[bands[iband]]+'-band', $
       /xstyle, yrange=[0.5,100], /ystyle, charsize=csize

      if (iband EQ 0 AND keyword_set(plottitle)) then $
       xyouts, plotmag[1], 106.0, 'S/N for '+plottitle, align=0.5, $
        charsize=csize

      ;----------
      ; Plot the fiducial line (a thick blue line)

      djs_oplot, plotmag, 10^poly(plotmag, fiducialfits[*,bands[iband]]), $
       color='blue', thick=5

      ;----------
      ; Plot the data points, (S/N) vs. magnitude.
      ; Identify which points fall above and below this fit,
      ;  color code green for positive residuals, red for negative.
      ; Use different symbols for each spectrograph.
      ; Plot these first so that the lines are visibly plotted on top.

      psymvec = (plugc.spectrographid EQ 1) * 7 $
       + (plugc.spectrographid EQ 2) * 6 
      colorvec = replicate('red', nobj)
      ipos = where(sndiff[igood] GE 0)
      if (ipos[0] NE -1) then colorvec[ipos] = 'green'
      djs_oplot, thismag[igood], snvec[iband,igood] > 0.6, $
       psym=psymvec, symsize=symsize, color=colorvec

      ; Now overplot the fit line
      if (keyword_set(afit)) then $
       djs_oplot, plotmag, 10^poly(plotmag, afit)

      ;----------
      ; Plot a fit to the data in the range specified by FITMAG,
      ; independently for each spectrograph.
      ; Also, draw an arrow that terminates at the S/N at the magnitude
      ; where we measure the canonical (S/N)^2.

      if (keyword_set(fitmag)) then myfitmag = fitmag

      snoise2 = fltarr(2)
      xloc = snmag[bands[iband]]
      if (s1[0] NE -1) then begin
         afit1 = fitsn(thismag[s1], snvec[iband,s1], fitmag=myfitmag, $
          colorband=bandnames[bands[iband]])
         if (keyword_set(afit1)) then begin
            snoise2[0] = 10^(2.0 * poly(snmag[bands[iband]],afit1)) 
            yloc = 10^poly(xloc,afit1)
            djs_arrow, xloc, 20, xloc, yloc, /data
         endif
      endif
      if (s2[0] NE -1) then begin
         afit2 = fitsn(thismag[s2], snvec[iband,s2], fitmag=myfitmag, $
          colorband=bandnames[bands[iband]])
         if (keyword_set(afit2)) then begin
            snoise2[1] = 10^(2.0 * poly(snmag[bands[iband]],afit2)) 
            yloc = 10^poly(xloc,afit2)
            djs_arrow, xloc, 20, xloc, yloc, /data
         endif
      endif

      ylimits = 10^[0.80 * !y.crange[0] + 0.20 * !y.crange[1], $
                 0.35 * !y.crange[0] + 0.65 * !y.crange[1] ]
      oplot, [myfitmag[0],myfitmag[0]], ylimits, linestyle=1
      oplot, [myfitmag[1],myfitmag[1]], ylimits, linestyle=1

      ;----------
      ; Label the plot

      if (keyword_set(afit)) then $
       xyouts, plotmag[0]+0.5, 0.9, string(format='(a,f6.3,f7.3,a)', $
        'log S/N = ', afit, slopelabel[bands[iband]]), charsize=textsize

      if (keyword_set(sigma)) then $
       xyouts, plotmag[0]+0.5, 1.5, string(format='(a,f4.2)','Stdev=', sigma), $
        charsize=textsize

      djs_xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.3], $
       10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.93], $
       snlabel[bands[iband]], charsize=textsize

      djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.57], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.85], psym=7, $
         symsize=symsize
      djs_xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.60], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.83], $
         string(format='("Spec1: ", f5.1)', snoise2[0]), charsize=textsize

      djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.57], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.75], psym=6, $
         symsize=symsize
      djs_xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.60], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.73], $
         string(format='("Spec2: ", f5.1)', snoise2[1]), charsize=textsize

      splog, snlabel[bands[iband]], snoise2, format='(a20, 2(f10.3))'

      snplate[*,iband] = snoise2

      ;------------------------------------------------------------------------
      ; 1st PAGE PLOT 2: Throughput deviations plotted on the focal plane
      ;------------------------------------------------------------------------

      plot, [0], [0], /nodata, xtickname=xtickname, $
       xtitle=xtitle2, ytitle='Y [mm]', $
       xrange=[-320,320], yrange=[-320,320], xstyle=1, ystyle=1, charsize=csize
      if (ngood GT 0) then begin
         colorvec = (sndiff[igood] GE 0) * djs_icolor('green') $
          + (sndiff[igood] LT 0) * djs_icolor('red')
         symvec = (abs(sndiff[igood]) * 5 < 2) > 0.1
         djs_oplot, plugc[igood].xfocal, plugc[igood].yfocal, $
          symsize=symvec, color=colorvec, psym=2
      endif
   endfor

   !p.multi = pmulti
   !y.margin = ymargin
   !y.omargin = yomargin

   if (NOT keyword_set(synthmag)) then return

   ;---------------------------------------------------------------------------
   ; PAGE 2: Plot spectro-mag vs. PHOTO-mag
   ; Loop over each band in the plot
   ;---------------------------------------------------------------------------

   pmulti = !p.multi
   ymargin = !y.margin
   yomargin = !y.omargin

   !p.multi = [0,2,nbands]
   !y.margin = [1,0]
   !y.omargin = [5,3]

   magdiff = synthmag - plugc.mag
   for iband=0, nbands-1 do begin

      if (iband LT nbands-1) then begin
         xtickname = strarr(20)+' '
         xtitle1 = ''
         xtitle2 = ''
      endif else begin
         xtickname = ''
         xtitle1 = 'Fiber ID'
         xtitle2 = 'X [mm]'
      endelse
      psym = 1

      plot, [0], [0], /nodata, $
       xtickname=xtickname, xrange=[0,641], $
       xtitle=xtitle1, ytitle='(Spectro-PHOTO) '+bandnames[bands[iband]]+'-mag', $
       /xstyle, yrange=[-0.6,0.6], /ystyle, charsize=csize
      oplot, !x.crange, [0,0]
      if (ngood GT 0) then begin
         thisdiff = magdiff[bands[iband],igood]
         colorvec = (thisdiff LT 0) * djs_icolor('green') $
          + (thisdiff GE 0) * djs_icolor('red')
         symvec = (abs(thisdiff) * 5 < 2) > 0.1
         djs_oplot, igood+1, thisdiff, $
          symsize=symvec, color=colorvec, psym=psym
      endif

      if (iband EQ 0 AND keyword_set(plottitle)) then $
       xyouts, !x.crange[1], 1.03*!y.crange[1]-0.03*!y.crange[0], $
        'Throughput for '+plottitle, align=0.5, charsize=csize

      plot, [0], [0], /nodata, xtickname=xtickname, $
       xtitle=xtitle2, ytitle='Y [mm]', $
       xrange=[-320,320], yrange=[-320,320], /xstyle, /ystyle, charsize=csize
      if (ngood GT 0) then begin
         djs_oplot, plugc[igood].xfocal, plugc[igood].yfocal, $
          symsize=symvec, color=colorvec, psym=psym
      endif
   endfor

   !p.multi = pmulti
   !y.margin = ymargin
   !y.omargin = yomargin

   ;------------------------------------------------------------------------
   ; Plots of spectro mags vs. PHOTO mags
   ;------------------------------------------------------------------------

   qstd = strmatch(plugc.objtype, '*STD*') 
   i1 = where(plugc.spectrographid EQ 1 AND qstd AND plugc.mag[2] NE 0)
   i2 = where(plugc.spectrographid EQ 2 AND qstd AND plugc.mag[2] NE 0)
   plotsn1, plugc, synthmag, i1, i2, plottitle=plottitle, objtype='Std-stars'

   qgal = strmatch(plugc.objtype, 'GALAXY*')
   i1 = where(plugc.spectrographid EQ 1 AND qgal AND plugc.mag[2] NE 0)
   i2 = where(plugc.spectrographid EQ 2 AND qgal AND plugc.mag[2] NE 0)
   plotsn1, plugc, synthmag, i1, i2, plottitle=plottitle, objtype='Galaxies'

   qstellar = qgal EQ 0
   i1 = where(plugc.spectrographid EQ 1 AND qstellar AND plugc.mag[2] NE 0)
   i2 = where(plugc.spectrographid EQ 2 AND qstellar AND plugc.mag[2] NE 0)
   plotsn1, plugc, synthmag, i1, i2, plottitle=plottitle, objtype='Stars+QSOs'

   ;----------
   ; Close plot file

   if (keyword_set(plotfile)) then begin
      device, /close
      set_plot, 'x'
   endif

   return
end
;------------------------------------------------------------------------------
