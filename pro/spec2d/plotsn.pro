;+
; NAME:
;   plotsn
;
; PURPOSE:
;   Plot S/N and residuals in up to 3 bands
;
; CALLING SEQUENCE:
;   plotsn, snvec, plugmap, [ filter=, plotmag=, plottitle=, snmin=, $
;    plotfile=, synthmag=, snplate=, dered_snplate=, specsnlimit=, _EXTRA= ]
;
; INPUTS:
;   snvec      - S/N array [NBAND,NFIBER]
;   plugmap    - Plugmap structure [NFIBER]
;
; OPTIONAL KEYWORDS:
;   filter     - Filter names; default to 'g','r','i'
;   plotmag    - Magnitude range for plotting; default to [16,23], but
;                extend the range to include FITMAG if necessary.
;   snmin      - Minimum S/N value to use in fits; default to 0.5
;   plottitle  - Title for top of plot
;   plotfile   - Optional plot file
;   synthmag   - Vector of synthetic magnitudes dimensionsed [5,NFIBER];
;                if set, then make more than the first page of plots
;   _EXTRA     - Keywords for FITSN, such as SNCODE
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   snplate    - Best fit (S/N)^2 at fiducial magnitude(s); array of [2,NBAND]
;                to give the number in each spectrograph and each filter
;   dered_snplate - Best fit (S/N)^2 extinction corrected like SOS pipeline
;   specsnlimit- Returned from FITSN
;
; COMMENTS:
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
; INTERNAL SUPPORT ROUTINES:
;   plotsn_good()
;   plotsn1
;
; REVISION HISTORY:
;   28-Jan-2003  Add E(B-V) keyword to allow for spectra which have 
;                foreground reddening removed
;   20-Oct-2002  Plot range on synthmag vs fiber mag plot changed by C. Tremonti
;   15-Apr-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
; Select good objects, and within FITMAG is specified
function plotsn_good, plugmap, jband, snvec, iband, igood, s1, s2, $
 fitmag=fitmag, snmin=snmin

   qgood = strtrim(plugmap.objtype,2) NE 'SKY' AND plugmap.mag[jband] GT 0 $
    AND snvec[iband,*] GT snmin
   if (keyword_set(fitmag)) then $
    qgood *= (plugmap.mag[jband] GT fitmag[0] $
     AND plugmap.mag[jband] LT fitmag[1])
   igood = where(qgood, ngood)
   if (ngood LT 3) then $
    splog, 'Warning: Too few non-sky objects to plot in band #', iband

   ; The following variables are defined only for non-SKY objects...
   s1 = where(plugmap.spectrographid EQ 1 AND qgood)
   s2 = where(plugmap.spectrographid EQ 2 AND qgood)

   return, ngood
end
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
         djs_oplot, [plugc[i1].mag[2]], [yplot[i1]], $
          psym=psym, symsize=symsize
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
         djs_oplot, [plugc[i2].mag[2]], [yplot[i2]], $
          psym=psym, symsize=symsize
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
pro plotsn, snvec1, plugmap1, filter=filter1, plotmag=plotmag1, snmin=snmin1, $
 plottitle=plottitle, plotfile=plotfile, synthmag=synthmag, snplate=snplate, $
 dered_snplate=dered_snplate, specsnlimit=specsnlimit, _EXTRA=KeywordsForFitSN

   if (keyword_set(filter1)) then filter = filter1 $
    else filter = ['g','r','i']
   if (keyword_set(plotmag1)) then plotmag = plotmag1 $
    else plotmag = [16.0, 23.0]
   if (keyword_set(snmin1)) then snmin = snmin1 $
    else snmin = 0.5

   nband = n_elements(filter)
   nfibers = n_elements(plugmap1)

   if (size(snvec1,/n_dimen) EQ 1) then snvec = reform(snvec1,1,nfibers) $
    else snvec = snvec1
   if (n_elements(snvec) NE nband*nfibers) then $
    message, 'Dimension of SNVEC not consistent with PLUGMAP'

   ; Use the CALIBFLUX magnitudes instead of MAG, if they exist
   ; from the call to READPLUGMAP() that generated this structure.
   ; For small fluxes, set MAG=0 so that they are ignored in the plots.
   ; Also, if CALIBFLUX_IVAR exist for any of the objects, then set MAG=0
   ; for the objects w/out it (e.g., for objects without calibObj photometry).
   plugmap = plugmap1
   if (tag_exist(plugmap,'CALIBFLUX')) then begin
      minflux = 0.1
      plugmap.mag =(22.5 - 2.5*alog10(plugmap.calibflux > minflux)) $
       * (plugmap.calibflux GT minflux)
;      if (total(plugmap.calibflux_ivar GT 0) GT 0) then $
;       plugmap.mag *= (plugmap.calibflux_ivar GT 0)
   endif

   ; Define jband as the index into the plugmap magnitudes,
   ; for example, jband=[1,2,3] if g,r,i filters
   jband = lonarr(nband)
   for iband=0, nband-1 do begin
      jband[iband] = where(filter[iband] EQ ['u','g','r','i','z'])
   endfor

   ;----------
   ; Open plot file

   if (keyword_set(plotfile)) then begin
      if (nband LE 2) then ysize=7.0 $
       else ysize=9.5
      dfpsplot, plotfile, /color, square=(nband LE 2)
   endif

   oldmulti = !p.multi  
   !p.multi = [0,2,nband]

   snplate = fltarr(2,nband)
   dered_snplate = fltarr(2,nband)

   pmulti = !p.multi
   ymargin = !y.margin
   yomargin = !y.omargin

   !p.multi = [0,2,nband]
   !y.margin = [1,0]
   !y.omargin = [5,3]
   csize = 1.0
   textsize = 1.0

   ;---------------------------------------------------------------------------
   ; Loop over each band in the plot
   ;---------------------------------------------------------------------------

   for iband=0, nband-1 do begin

      thismag = plugmap.mag[jband[iband]]

      ngood = plotsn_good(plugmap, jband[iband], snvec, iband, igood, s1, s2, $
       snmin=snmin)

      ;----------
      ; Fit the data as S/N vs. mag

      if (ngood GE 3) then $
       afit = fitsn(thismag[igood], snvec[iband,igood], sigma=sigma, $
        filter=filter[iband], specsnlimit=specsnlimit1, $
        _EXTRA=KeywordsForFitSN) $
      else $
       afit = fitsn([0], [0], sigma=sigma, $
        filter=filter[iband], specsnlimit=specsnlimit1, $
        _EXTRA=KeywordsForFitSN)
      if (iband EQ 0) then specsnlimit = specsnlimit1 $
       else specsnlimit = [specsnlimit, specsnlimit1]
      fitmag = specsnlimit1.fitmag
      snmag = specsnlimit1.snmag
      ; Residuals from this fit
      sndiff = alog10(snvec[iband,*]>0.01) - poly(thismag, afit)

      if (ngood GT 3) then begin
         ;---------------------------------------------------------------------
         ; 1st PAGE PLOT 1: (S/N) vs. magnitude
         ;---------------------------------------------------------------------

         ;----------
         ; Extend the plotting range if necessary to include FITMAG and SNMAG

         plotmag[0] = plotmag[0] < (fitmag[0]-0.5)
         plotmag[1] = plotmag[1] > (fitmag[1]+0.5)
         plotmag[0] = plotmag[0] < (snmag-0.5)
         plotmag[1] = plotmag[1] > (snmag+0.5)
         snlabel = '(S/N)^2 @ '+filter[iband]+' ='+string(snmag,format='(f7.2)')

         ;----------
         ; Set up the plot

         if (iband LT nband-1) then begin
            xtickname = strarr(20)+' ' ; Disable X axis on all but bottom plot
            xtickname = ''
            xtitle1 = ''
            xtitle2 = ''
         endif else begin
            xtickname = ''
            xtitle1 = 'mag'
            xtitle2 = 'X [mm]'
         endelse
         symsize = 0.4

         plot, thismag[igood], snvec[iband,igood], /nodata, /ylog, $
;          xtickname=xtickname, $
          xrange=plotmag, $
          xtitle=xtitle1, ytitle='S/N in '+filter[iband]+'-band', $
          /xstyle, yrange=[0.5,100], /ystyle, charsize=csize

         if (iband EQ 0 AND keyword_set(plottitle)) then $
          xyouts, plotmag[1], 106.0, 'S/N for '+plottitle, align=0.5, $
           charsize=csize

         ;----------
         ; Plot the fiducial line (black line)

;         djs_oplot, plotmag, 10^poly(plotmag, specsnlimit1.fiducial_coeff)

         ;----------
         ; Plot the data points, (S/N) vs. magnitude.

         ; Identify which points fall above and below this fit.
         ; Color code green for positive residuals, red for negative.
         ; Use different symbols for each spectrograph.
;         psymvec = (plugmap.spectrographid EQ 1) * 7 $
;          + (plugmap.spectrographid EQ 2) * 6 
;         colorvec = replicate('red', nfibers)
;         ipos = where(sndiff GE 0)
;         if (ipos[0] NE -1) then colorvec[ipos] = 'green'

         ; Color code spec1=cyan X, spec2=magenta square
         ; Use different symbols for each spectrograph.
         psymvec = (plugmap.spectrographid EQ 1) * 7 $
          + (plugmap.spectrographid EQ 2) * 6 
         colorvec = replicate('cyan', nfibers)
         ipos = where(plugmap.spectrographid EQ 2)
         if (ipos[0] NE -1) then colorvec[ipos] = 'magenta'

         ; Sort the points randomly, such that not all the spec2 points
         ; are plotted on top of the spec1 points
         isort = sort(randomu(1234,n_elements(thismag)))

         ; Plot the points that would fall off the bottom of the plot
         ; as points near the bottom
         djs_oplot, [thismag[isort]], $
          [snvec[iband,isort] > 1.1*(snmin>10^!y.crange[0])], $
          psym=psymvec[isort], symsize=symsize, color=colorvec[isort]

         ; Now overplot the fit line
;         if (keyword_set(afit)) then $
;          djs_oplot, plotmag, 10^poly(plotmag, afit)

         ;----------
         ; Plot a fit to the data in the range specified by FITMAG,
         ; independently for each spectrograph.
         ; Also, draw an arrow that terminates at the S/N at the magnitude
         ; where we measure the canonical (S/N)^2.

         afit1 = 0
         sig1 = 0
         if (s1[0] NE -1) then begin
            afit1 = fitsn(thismag[s1], snvec[iband,s1], $
             filter=filter[iband], _EXTRA=KeywordsForFitSN, $
             sigma=sig1, sn2=sn2, dered_sn2=dered_sn2)
          snplate[0,iband] = sn2
          dered_snplate[0,iband] = dered_sn2
         endif
         sig2 = 0
         afit2 = 0
         if (s2[0] NE -1) then begin
            afit2 = fitsn(thismag[s2], snvec[iband,s2], $
             filter=filter[iband], _EXTRA=KeywordsForFitSN, $
              sigma=sig2, sn2=sn2, dered_sn2=dered_sn2)
           snplate[1,iband] = sn2
           dered_snplate[1,iband] = dered_sn2
         endif

         ylimits = 10^[0.80 * !y.crange[0] + 0.20 * !y.crange[1], $
                    0.35 * !y.crange[0] + 0.65 * !y.crange[1] ]
         oplot, fitmag[0]+[0,0], ylimits, linestyle=1, thick=3
         oplot, fitmag[1]+[0,0], ylimits, linestyle=1, thick=3

         ;----------
         ; Label the plot

         ; Overplot the fit lines
         if (keyword_set(afit1)) then begin
            djs_oplot, plotmag, 10^poly(plotmag, afit1)
            xyouts, plotmag[0]+0.5, 1.35, string(format='(a,f6.3,f7.3,a)', $
             'log S/N = ', afit1, ' * '+filter[iband]), charsize=textsize
         endif
         if (keyword_set(afit2)) then begin
            djs_oplot, plotmag, 10^poly(plotmag, afit2), color='magenta'
            xyouts, plotmag[0]+0.5, 0.90, string(format='(a,f6.3,f7.3,a)', $
             'log S/N = ', afit2, ' * '+filter[iband]), charsize=textsize
         endif
         ; Overplot arrows at fiducial mag
         if (snplate[0,iband] GT 0) then $
          djs_arrow, snmag, sqrt(snplate[0,iband])*2, snmag, $
           sqrt(snplate[0,iband]), /data
         if (snplate[1,iband] GT 0) then $
          djs_arrow, snmag, sqrt(snplate[1,iband])*2, snmag, $
           sqrt(snplate[1,iband]), /data

         if (keyword_set(sig1)) then $
          xyouts, plotmag[0]+0.5, 0.60, $
           string(format='(a,f4.2,f6.2)','Stdev=', sig1, sig2), $
           charsize=textsize

         djs_xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.3], $
          10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.93], $
          snlabel, charsize=textsize

         djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.57], $
          10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.85], psym=7, $
          symsize=symsize
         djs_xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.60], $
          10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.83], $
          string(format='("Spec1: ", f5.1)', snplate[0,iband]), charsize=textsize

         djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.57], $
          10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.75], psym=6, $
          symsize=symsize, color='magenta'
         djs_xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.60], $
          10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.73], $
          string(format='("Spec2: ", f5.1)', snplate[1,iband]), charsize=textsize

         splog, snlabel, snplate[*,iband], format='(a20, 2(f10.3))'

         ;---------------------------------------------------------------------
         ; 1st PAGE PLOT 2: Throughput deviations plotted on the focal plane
         ;---------------------------------------------------------------------

         plot, [0], [0], /nodata, xtickname=xtickname, $
          xtitle=xtitle2, ytitle='Y [mm]', $
          xrange=[-350,350], yrange=[-350,350], xstyle=1, ystyle=1, charsize=csize
         ; Limit this plot to only those objects which are good
         ng2 = plotsn_good(plugmap, jband[iband], snvec, iband, ig2, s1, s2, $
          snmin=snmin)
         if (ng2 GT 0) then begin
            colorvec = (sndiff[ig2] GE 0) * djs_icolor('green') $
             + (sndiff[ig2] LT 0) * djs_icolor('red')
            symvec = (abs(sndiff[ig2]) * 5 < 2) > 0.2
            djs_oplot, [plugmap[ig2].xfocal], [plugmap[ig2].yfocal], $
             symsize=symvec, color=colorvec, psym=2
         endif

         djs_oplot, [-290], [290], symsize=2, psym=2, color='green'
         djs_xyouts, -260, 290, '+0.4 mag', color='green'
         djs_oplot, [-290], [250], symsize=2, psym=2, color='red'
         djs_xyouts, -260, 250, '-0.4 mag', color='red'

      endif
   endfor

   !p.multi = pmulti
   !y.margin = ymargin
   !y.omargin = yomargin

   if (NOT keyword_set(synthmag)) then begin
      if (keyword_set(plotfile)) then dfpsclose
      return
   endif

   ;---------------------------------------------------------------------------
   ; PAGE 2: Plot spectro-mag vs. PHOTO-mag
   ; Loop over each band in the plot
   ;---------------------------------------------------------------------------

   pmulti = !p.multi
   ymargin = !y.margin
   yomargin = !y.omargin

   !p.multi = [0,2,nband]
   !y.margin = [1,0]
   !y.omargin = [5,3]

   magdiff = synthmag - plugmap.mag
   radius = sqrt(plugmap.xfocal^2 + plugmap.yfocal^2)
   for iband=0, nband-1 do begin

      ngood = plotsn_good(plugmap, jband[iband], snvec, iband, igood, s1, s2, $
       snmin=snmin) ; fitmag=specsnlimit[iband].fitmag)
      if (ngood GE 3) then begin

      if (iband LT nband-1) then begin
         xtickname = strarr(20)+' '
         xtitle1 = ''
         xtitle2 = ''
      endif else begin
         xtickname = ''
         xtitle1 = 'Radius [mm]'
         xtitle2 = 'X [mm]'
      endelse
      psym = strmatch(plugmap.objtype, 'QSO*') * 1 $ ; Plus
       + strmatch(plugmap.objtype, 'GALAXY*') * 6 ; Square
      psym = psym + (psym EQ 0) * 2 ; Asterisk for anything else

      plot, [0], [0], /nodata, $
       xtickname=xtickname, xrange=[0,350], xtitle=xtitle1, $
       ytitle='(Spectro-PHOTO) '+filter[iband]+'-mag', $
       /xstyle, yrange=[-0.6,0.6], /ystyle, charsize=csize
      oplot, !x.crange, [0,0]
      if (ngood GT 0) then begin
         thisdiff = transpose(magdiff[jband[iband],igood])
         colorvec = (thisdiff LT 0) * djs_icolor('green') $
          + (thisdiff GE 0) * djs_icolor('red')
         symvec = (abs(thisdiff) * 5 < 2) > 0.1
         djs_oplot, [radius[igood]], [thisdiff], $
          symsize=0.5, color=colorvec, psym=psym[igood]
      endif
      oplot, [15], [0.5], psym=6
      oplot, [15], [0.4], psym=1
      oplot, [15], [0.3], psym=2
      xyouts, 30, 0.5, 'Galaxy'
      xyouts, 30, 0.4, 'QSO'
      xyouts, 30, 0.3, 'other'

      if (iband EQ 0 AND keyword_set(plottitle)) then $
       xyouts, !x.crange[1], 1.03*!y.crange[1]-0.03*!y.crange[0], $
        'Throughput for '+plottitle, align=0.5, charsize=csize

      plot, [0], [0], /nodata, xtickname=xtickname, $
       xtitle=xtitle2, ytitle='Y [mm]', $
       xrange=[-350,350], yrange=[-350,350], /xstyle, /ystyle, charsize=csize
      if (ngood GT 0) then begin
         djs_oplot, [plugmap[igood].xfocal], [plugmap[igood].yfocal], $
          symsize=symvec, color=colorvec, psym=psym[igood]
      endif

      djs_oplot, [-290], [290], symsize=2, psym=2, color='green'
      djs_xyouts, -260, 290, '+0.4 mag', color='green'
      djs_oplot, [-290], [250], symsize=2, psym=2, color='red'
      djs_xyouts, -260, 250, '-0.4 mag', color='red'

   endif
   endfor

   !p.multi = pmulti
   !y.margin = ymargin
   !y.omargin = yomargin

   ;------------------------------------------------------------------------
   ; Plots of spectro mags vs. PHOTO mags
   ;------------------------------------------------------------------------

   qstd = strmatch(plugmap.objtype, '*STD*') 
   i1 = where(plugmap.spectrographid EQ 1 AND qstd AND plugmap.mag[2] NE 0)
   i2 = where(plugmap.spectrographid EQ 2 AND qstd AND plugmap.mag[2] NE 0)
   plotsn1, plugmap, synthmag, i1, i2, plottitle=plottitle, objtype='Std-stars'

   qgal = strmatch(plugmap.objtype, 'GALAXY*')
   i1 = where(plugmap.spectrographid EQ 1 AND qgal AND plugmap.mag[2] NE 0)
   i2 = where(plugmap.spectrographid EQ 2 AND qgal AND plugmap.mag[2] NE 0)
   plotsn1, plugmap, synthmag, i1, i2, plottitle=plottitle, objtype='Galaxies'

   qstellar = qgal EQ 0
   i1 = where(plugmap.spectrographid EQ 1 AND qstellar AND plugmap.mag[2] NE 0)
   i2 = where(plugmap.spectrographid EQ 2 AND qstellar AND plugmap.mag[2] NE 0)
   plotsn1, plugmap, synthmag, i1, i2, plottitle=plottitle, objtype='Stars+QSOs'

   if (keyword_set(plotfile)) then dfpsclose

   return
end
;------------------------------------------------------------------------------
