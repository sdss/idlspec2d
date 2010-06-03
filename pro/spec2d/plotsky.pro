pro plotsky, plate, expnum, camname=camname, plotfile=plotfile1, $
 _EXTRA=KeywordsForPlot

   if (NOT keyword_set(camname)) then camname = 'r1'
   if (NOT keyword_set(plate) OR NOT keyword_set(expnum)) then return

   filename = string(camname, expnum, $
    format='("spFrame-",a2,"-",i8.8,".fits*")')
   fullname = (findfile(filepath(filename, root_dir=getenv('BOSS_SPECTRO_REDUX'), $
    subdir=string(plate,format='(i4.4)')), count=ct))[0]
   if (ct EQ 0) then begin
      splog, 'File not found ' + filename
      return
   endif

   if (keyword_set(plotfile1)) then begin
      if (size(plotfile1,/tname) EQ 'STRING') then plotfile = plotfile1 $
       else plotfile = string(plate, camname, expnum, $
        format='("sky-", i4.4, "-", a2, "-", i8.8, ".ps")')
   endif

   ifiber = 160

   hdr = headfits(fullname)
   exptime = sxpar(hdr,'EXPTIME') > 1
   mjd = sxpar(hdr,'MJD') > 1
   wset = mrdfits(fullname, 3)
   skyimg = mrdfits(fullname, 6)
   traceset2xy, wset, xx, yy
   wave = 10.d0^yy[*,ifiber]
   skyvec = skyimg[*,ifiber] / exptime

   plottitle = 'Sky for PLATE=' + strtrim(string(plate),2) $
    + '/' + string(mjd,format='(i5)') $
    + ' Cam=' + camname $
    + ' Exp=' + strtrim(string(expnum),2)
   xtitle = 'Wavelength [Ang]'
   ytitle = 'Sky Flux [e-/sec]'
   xrange = minmax(wave)
   yrange = [0, 1.1*max(skyvec)]
   xplot = total(xrange * [0.95,0.05])
   yplot = total(yrange * [0.10,0.90])
   expstring = 'EXPTIME = ' + string(exptime, format='(i4)') + ' sec'

   csize = 1.4
   if (keyword_set(plotfile)) then begin
      dfpsplot, plotfile, /square
      plot, wave, skyvec, charsize=csize, $
       xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
       xtitle=xtitle, ytitle=ytitle, title=plottitle, $
       _EXTRA=KeywordsForPlot
      xyouts, xplot, yplot, expstring, charsize=csize
      dfpsclose
   endif else begin
      splot, wave, skyvec, charsize=csize, $
       xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
       xtitle=xtitle, ytitle=ytitle, title=plottitle, $
       _EXTRA=KeywordsForPlot
      sxyouts, xplot, yplot, expstring, charsize=csize
   endelse

   return
end

