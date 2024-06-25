;+
; NAME:
;   apo_plotbias
;
; PURPOSE:
;   Plot the histogram of bias values for all 4 cameras of a single exposure
;
; CALLING SEQUENCE:
;   apo_plotbias, expnum, [ plotfile= ]
;
; INPUTS:
;   expnum     - Exposure number
;
; OPTIONAL INPUTS:
;   plotfile   - Plot file; if set, then send plot to this PostScript file
;                rather than to the default (X) display.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The histogram of bias values is plotted for all (4) camera files that
;   match the given exposure number.
;
;   A fiducial line is drawn as a thick blue line.  This line approximates
;   what we expect to see for each camera.
;
;   If $BOSS_SPECTRO_DATA is not set, then it is assumed to be
;     /data/spectro
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_filepath()
;   djs_icolor()
;   djs_xyouts
;   fileandpath()
;   headfits()
;   plothist
;   quickbias
;   sdssproc
;   struct_append()
;   sxpar()
;
; REVISION HISTORY:
;   06-Dec-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro apo_plotbias, expnum, plotfile=plotfile

   if (n_params() LT 1) then begin
      print, 'Syntax - apoplotbias, expnum, [plotfile= ]'
      return
   endif

   rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   if (NOT keyword_set(rawdata_dir)) then $
    rawdata_dir = '/data/spectro'

   expstr = string(expnum, format='(i8.8)')
   filename = 'sdR-*-' + expstr +  '.fit*'
   fullname = findfile( djs_filepath(filename, root_dir=rawdata_dir, $
    subdirectory='*'), count=nfile)

   ;----------
   ; Read the MJD from the header of the first file

   hdr = headfits(fullname[0])
   thismjd = sxpar(hdr, 'MJD')
   mjdstr = strtrim(string(thismjd),2)

   ;----------
   ; Start the plot

   xrange = [-1,5]

   yrange = 10.^[0,7]
   bin = 0.2
   csize = 2.0
   plotcolor = ['default', 'blue', 'green', 'red']
   nplotcolor = n_elements(plotcolor)

   if (keyword_set(plotfile)) then $
    dfpsplot, plotfile, /color

   plot, [0], [0], /nodata, charsize=csize, $
    xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, /ylog, $
    xtitle='Log(ADU)', ytitle='Number of pixels', $
    title='HISTOGRAM OF BIAS VALUES (MJD='+mjdstr+')'

   ;----------
   ; Overplot a fiducial line as a thick blue curve

   nplot = 1000
   xplot = (xrange[1] - xrange[0]) * findgen(nplot) / nplot + xrange[0]
   xmean = [0.6, 2.0]
   xsig = [0.5, 2.0]
   ymax = [6.0, 2.0]
   yplot = fltarr(nplot)
   for j=0, n_elements(xmean)-1 do $
    yplot = yplot + 10.^(ymax[j] - ((xplot - xmean[j]) / xsig[j])^2)
   oplot, xplot, yplot, color=djs_icolor('blue'), thick=3

   ;----------
   ; Loop through each of the 4 files (one for each camera), and plot
   ; the histogram of values.

   for ifile=0, nfile-1 do begin
      sdssproc, fullname[ifile], imgflux, imgivar
      igood = where(imgivar NE 0, ngood)

      icolor = djs_icolor(plotcolor[ifile MOD nplotcolor])
      if (ngood GT 10) then $
       plothist, alog10(imgflux[igood]>1), charsize=csize, $
        bin=bin, /noerase, color=icolor, $
        xrange=xrange, yrange=yrange, xstyle=5, ystyle=5, /ylog

      djs_xyouts, 2.0, 10.^(6.5 - 0.5*ifile), charsize=csize, $
       fileandpath(fullname[ifile]), color=icolor
   endfor

   ;----------
   ; Close the plot file

   if (keyword_set(plotfile)) then $
    dfpsclose

   ;----------
   ; Compute and print statistics

   print
   print, 'Now re-computing percentiles...'
   for ifile=0, nfile-1 do $
    rstruct = struct_append(rstruct, quickbias(fullname[ifile]))

   print
   print, 'Filename           ', $
    '  02%     05%     10%     50%     90%     95%     98%   '
   print, '-------------------', $
    '  ------  ------  ------  ------  ------  ------  ------'
   for ifile=0, nfile-1 do $
    print, fileandpath(fullname[ifile]), $
     rstruct[ifile].percentile[2], $
     rstruct[ifile].percentile[5], $
     rstruct[ifile].percentile[10], $
     rstruct[ifile].percentile[50], $
     rstruct[ifile].percentile[90], $
     rstruct[ifile].percentile[95], $
     rstruct[ifile].percentile[98], $
     format='(a19,7f8.1)'
   print

   return
end
;------------------------------------------------------------------------------
