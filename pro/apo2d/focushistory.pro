;+
; NAME:
;   focushistory
;
; PURPOSE:
;   Plot the history of wavelength focus based upon the reduced spPlate files.
;
; CALLING SEQUENCE:
;   focushistory, [ mjdrange= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjdrange   - 2-element vector of plotting range in MJD; default to
;                using all MJDs.
;
; OUTPUT:
;
; COMMENTS:
;   A single PostScript file is created "Focushistory-$MJDSTART-$MJDEND.ps".
;   The plate list must exist for the PLATELIST procedure, and the spPlate
;   files should be in $BOSS_SPECTRO_REDUX.
;
; EXAMPLES:
;   Make plots of history of wavelength focus for all time:
;     IDL> focushistory
;
;   Make plots of history of wavelength focus between MJD 52000 and MJD 52100:
;     IDL> focusistory, mjdrange=[5200,52100]
;
; BUGS:
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_icolor()
;   fileandpath()
;   mjd2datelist
;   platelist
;   readspec
;   splog
;   yanny_free
;   yanny_par()
;   yanny_read
;
; REVISION HISTORY:
;   21-Mar-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro focushistory, mjdrange=mjdrange1

   fthresh = 1.00 ; Threshhold for labeling out-of-focus plates
   rescale = [0.88,1.00,0.98,1.00]

   ;----------
   ; Get the list of plates and trim to reduced plates

   platelist, plist=plist
   if (NOT keyword_set(plist)) then begin
      print, 'No reduced plates found'
      return
   endif

   idone = where(strtrim(plist.statuscombine,2) EQ 'Done')
   if (idone[0] EQ -1) then begin
      print, 'No reduced plates found'
      return
   endif
   plist = plist[idone]

   ;----------
   ; Trim the list of plates to the specified MJD range

   if (keyword_set(mjdrange1)) then begin
      if (n_elements(mjdrange1) NE 2) then begin
         print, 'MJDRANGE must be 2-element vector, e.g. MJDRANGE=[52000,52100]'
         return
      endif
      indx = where(plist.mjd GE mjdrange1[0] AND plist.mjd LE mjdrange1[1], nfile)
      if (nfile EQ 0) then begin
         print, 'No reduced plates within specified MJD range'
         return
      endif
      print, 'Trimming to ', nfile, ' plates within specified MJD range'
      plist = plist[indx]
      mjdrange = mjdrange1
   endif else begin
      mjdrange = [min(plist.mjd),max(plist.mjd)]
   endelse

   ;----------
   ; Initialize the log file

   logfile = string(min(mjdrange), max(mjdrange), $
    format='("Focushistory-",i5.5,"-",i5.5,".log")')
   splog, file=logfile
   splog, 'PLATE  MJD   Foc-b1 Foc-r1 Foc-b2 Foc-r2 FocMax Temp '
   splog, '-----  ----- ------ ------ ------ ------ ------ -----'

   ;----------
   ; Read through data for each plate

   plist = plist[sort(plist.mjd)] ; Sort by MJD
   nfile = n_elements(plist)
   wsigarr = fltarr(4,nfile)
   wsigdev = fltarr(4,nfile)
   wsigmin = fltarr(4,nfile)
   wsigmax = fltarr(4,nfile)
   wsigworst = fltarr(nfile)

   for ifile=0, nfile-1 do begin
;      print, 'Reading file ', ifile+1, ' of ', nfile, ': ', $
;       plist[ifile].plate, '/', plist[ifile].mjd
      readspec, plist[ifile].plate, mjd=plist[ifile].mjd, 1, wave=wave
      readspec, plist[ifile].plate, mjd=plist[ifile].mjd, disp=dispimg
      dims = size(dispimg,/dimens)

      indx1 = where(wave GT 4000 AND wave LT 5500)
      indx2 = where(wave GT 6500 AND wave LT 9000)

      ; Mean dispersions for b1,r1,b2,r2
      i_b1 = where(dispimg[indx1,0:319] GT 0)
      djs_iterstat, (dispimg[indx1,0:319])[i_b1], mean=mn_b1, sigma=sg_b1
      min_b1 = min((dispimg[indx1,0:319])[i_b1], max=max_b1)

      i_r1 = where(dispimg[indx2,0:319] GT 0)
      djs_iterstat, (dispimg[indx2,0:319])[i_r1], mean=mn_r1, sigma=sg_r1
      min_r1 = min((dispimg[indx2,0:319])[i_r1], max=max_r1)

      i_b2 = where(dispimg[indx1,320:639] GT 0)
      djs_iterstat, (dispimg[indx1,320:639])[i_b2], mean=mn_b2, sigma=sg_b2
      min_b2 = min((dispimg[indx1,320:639])[i_b2], max=max_b2)

      i_r2 = where(dispimg[indx2,320:639] GT 0)
      djs_iterstat, (dispimg[indx2,320:639])[i_r2], mean=mn_r2, sigma=sg_r2
      min_r2 = min((dispimg[indx2,320:639])[i_r2], max=max_r2)

      wsigarr[*,ifile] = [mn_b1, mn_r1, mn_b2, mn_r2] * rescale
      wsigdev[*,ifile] = [sg_b1, sg_r1, sg_b2, sg_r2] * rescale
      wsigmin[*,ifile] = [min_b1, min_r1, min_b2, min_r2] * rescale
      wsigmax[*,ifile] = [max_b1, max_r1, max_b2, max_r2] * rescale

      wsigworst[ifile] = max(wsigarr[*,ifile])

      splog, plist[ifile].plate, plist[ifile].mjd, $
       wsigarr[*,ifile], wsigworst[ifile], plist[ifile].airtemp, $
       format='(i5,i7,5f7.2,f6.1)'
   endfor

   splog, /close

   ;----------
   ; Initialize the plot file

   mjdplot = plist.tai / (24.D*3600.D)

   plotfile = string(min(mjdrange), max(mjdrange), $
    format='("Focushistory-",i5.5,"-",i5.5,".ps")')
   dfpsplot, plotfile, /color

   xrange = minmax(plist.mjd) + [-30,30]
   yrange = [0.60,1.40]
   mjd2datelist, min(plist.mjd), max(plist.mjd), step='year', $
    mjdlist=mjdlist, datelist=datelist
   nplot = n_elements(mjdlist) - 1

   ; Set multi-plot format
   !p.multi = [0,1,2]

   ;----------
   ; Make one plotting panel per calander year

   for iplot=0, nplot-1 do begin

      plot, [0], [0], /nodata, $
       xrange=mjdlist[[iplot,iplot+1]], /xstyle, yrange=yrange, /ystyle, $
       xtickformat='(i5)', $
       xtitle='MJD', ytitle='Worst Focus [pix]', $
       title='Focus History (Year=' + strmid(datelist[iplot],7)+')'
      oplot, !x.crange, [0.85,0.85]

      mjd2datelist, mjdlist[iplot], mjdlist[iplot+1], step='month', $
       mjdlist=mjd1, datelist=date1
      for i=1, n_elements(mjd1)-2 do $
       oplot, [mjd1[i],mjd1[i]], !y.crange, linestyle=1
      for i=1, n_elements(mjd1)-2 do $
       xyouts, mjd1[i], 0.95*!y.crange[0] + 0.05*!y.crange[1], $
        strmid(date1[i],0,6), orient=90, align=0

      indx = where(mjdplot GE mjdlist[iplot] AND mjdplot LE mjdlist[iplot+1])
      if (indx[0] NE -1) then begin
         oplot, mjdplot[indx], wsigworst[indx], psym=4
         ibad = where(wsigworst[indx] GT fthresh)
         if (ibad[0] NE -1) then $
          xyouts, mjdplot[indx[ibad]], wsigworst[indx[ibad]], $
           ' '+strtrim(string(plist[indx[ibad]].plate),2)
      endif
   endfor

   ;----------
   ; Make a final panel with the histogram of focus values per CCD

   binsz = 0.01
   camname = ['b1','r1','b2','r2']
   colorvec = ['blue','red','green','magenta']
   plothist, [wsigworst,wsigworst], /nodata, $ ; Double count to get big Y limit
    bin=binsz, xrange=[0.7,1.5], /xstyle, $
    xtitle='Focus [pix]', ytitle='Number of Plates', $
    title='Distribution of Focus Values'
   dy = 0.05 * (!y.crange[1] - !y.crange[0])
   yplot = !y.crange[1] - 2*dy
   for iccd=0, 3 do begin
      plothist, wsigarr[iccd,*], bin=binsz, /overplot, $
       xrange=!x.crange, yrange=!y.crange, xstyle=5, ystyle=4, $
       color=djs_icolor(colorvec[iccd])
      text = 'Camera='+camname[iccd] + ' (scaled by ' $
       + string(rescale[iccd], format='(f4.2)') + ')'
      xyouts, 1.10, yplot, text, color=djs_icolor(colorvec[iccd])
      yplot = yplot - dy
   endfor

   ;----------
   ; Plot worst focus vs. air temp

   indx = where(plist.airtemp GT -20 AND plist.airtemp NE 0)
   if (indx[0] NE -1) then $
    plot, [plist[indx].airtemp], [wsigworst[indx]], /ynozero, psym=4, $
     xtitle='Air Temperature [deg C]', ytitle='Worst Focus', $
     title='Focus vs. Temperature'

   ;----------
   ; Close the plot file

   dfpsclose
   !p.multi = 0

   ;----------
   ; Create a save-set with these data

   ssfile = string(min(mjdrange), max(mjdrange), $
    format='("Focushistory-",i5.5,"-",i5.5,".ss")')
   save, file=ssfile

   return
end
;------------------------------------------------------------------------------
