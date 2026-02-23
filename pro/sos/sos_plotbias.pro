;+
; NAME:
;   sos_plotbias
;
; PURPOSE:
;   Plot the histogram of bias values for all 4 cameras of a single exposure
;
; CALLING SEQUENCE:
;   sos_plotbias, expnum, [ plotfile= ]
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
pro sos_plotbias, expnum, plotfile=plotfile, tolog=tolog

   if (n_params() LT 1) then begin
      print, 'Syntax - apoplotbias, expnum, [plotfile= ]'
      return
   endif

   if strmatch(getenv('OBSERVATORY'), 'LCO',/fold_case) then begin
    BOSS_SPECTRO_DATA = 'BOSS_SPECTRO_DATA_S'
   endif else BOSS_SPECTRO_DATA = 'BOSS_SPECTRO_DATA_N'
   rawdata_dir = getenv(BOSS_SPECTRO_DATA)

   if (NOT keyword_set(rawdata_dir)) then $
    rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   
   if (NOT keyword_set(rawdata_dir)) then $
    rawdata_dir = '/data/spectro'

   expstr = string(expnum, format='(i8.8)')
   filename = 'sdR-*-' + expstr +  '.fit*'
   fullname = findfile( djs_filepath(filename, root_dir=rawdata_dir, $
    subdirectory='*'), count=nfile)

   ;----------
   ; Read the MJD from the header of the first file
   
   if nfile eq 0 then return
   hdr = headfits(fullname[0])
   thismjd = sxpar(hdr, 'MJD')
   mjdstr = strtrim(string(thismjd),2)

   ;----------
   ; Start the plot

   if strmatch(sxpar(hdr,'CARTID'), '*FPS-S*', /fold_case) then obs='LCO' else obs='APO'
   if strmatch(obs, 'apo',/fold_case) eq 1 then begin
        xrange = [-1,5]
        xmean = [-.05, 1.5]
        xsig = [0.5, 2.0]
        ymax = [6.0, 1.2]
   endif else begin
        xrange = [-1,5]
        xmean = [0.1, 0.8]
        xsig = [0.4, 1.5]
        ymax = [6.0, 1.5]
   endelse


   yrange = 10.^[0,7]
   bin = 0.2
   csize = 1.8
   plotcolor = ['default', 'blue', 'green', 'red']
   nplotcolor = n_elements(plotcolor)

   if (keyword_set(plotfile)) then $
    dfpsplot, plotfile, /color

   plot, [0], [0], /nodata, charsize=csize, $
    xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, /ylog, $
    xtitle='', ytitle='Number of pixels', $
    title='HISTOGRAM OF BIAS VALUES (MJD='+mjdstr+')'
   XYOUTS, .95, .04, 'Log(ADU)', ALIGN=1.0, /normal, charsize=csize

   ;----------
   ; Overplot a fiducial line as a thick blue curve

   nplot = 1000
   xplot = (xrange[1] - xrange[0]) * findgen(nplot) / nplot + xrange[0]
   yplot = fltarr(nplot)
   for j=0, n_elements(xmean)-1 do $
    yplot = yplot + 10.^(ymax[j] - ((xplot - xmean[j]) / xsig[j])^2)
   oplot, xplot, yplot, color=djs_icolor('blue'), thick=3

   ;----------
   ; Loop through each of the 4 files (one for each camera), and plot   ; the histogram of values.

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

   if keyword_set(tolog) then begin
       splog, 'Now re-computing percentiles...'
   endif else print, 'Now re-computing percentiles...'
   for ifile=0, nfile-1 do $
        rstruct = struct_append(rstruct, quickbias(fullname[ifile]))

   if keyword_set(plotfile) then begin
    text_xpos =0.01 ; -0.01 ;xrange[0]-1  ; Start from the left side of the plot area
    text_ypos = 0.05; yrange[0] / 10.0  ; Place below the y-axis range (log scale)
    offset = 0.009 ;0.0012  ; Vertical spacing between text lines
    csize = .8
    ; Write headers to the PostScript file
    djs_xyouts, text_xpos, text_ypos, string('Filename',FORMAT='(a21)')+ ' '+string('02%',FORMAT='(a8)')+' '+$
                                      string('05%',FORMAT='(a8)') +' '+ string('10%',FORMAT='(a8)')+' '+$
                                      string('50%',FORMAT='(a8)')+ string('90%',FORMAT='(a8)')+$
                                      string('95%',FORMAT='(a8)')+ string('98%',FORMAT='(a8)'), charsize=csize, /normal, ALIGNMENT=0
    text_ypos -= offset
    djs_xyouts, -.01, text_ypos, '---------------  -----  -----  ----  -----  -----  -----  ----', charsize=csize, /normal, ALIGNMENT=0
    text_ypos -= offset

    ; Loop through each file and add its statistics to the PostScript file
    for ifile=0, nfile-1 do begin
        djs_xyouts, text_xpos, text_ypos, $
            string(fileandpath(fullname[ifile]), FORMAT='(a19)') + ' '+$
            string(rstruct[ifile].percentile[1], FORMAT='(F8.1)') +' '+$
            string(rstruct[ifile].percentile[4], FORMAT='(F8.1)') +' '+$
            string(rstruct[ifile].percentile[9], FORMAT='(F8.1)')+' '+$
            string(rstruct[ifile].percentile[49], FORMAT='(F8.1)')+' '+$
            string(rstruct[ifile].percentile[89], FORMAT='(F8.1)')+' '+$
            string(rstruct[ifile].percentile[94], FORMAT='(F8.1)')+' '+$
            string(rstruct[ifile].percentile[97], FORMAT='(F8.1)'), $
           ; format='(a19,7f8.1)', $
            charsize=csize, /normal, ALIGNMENT=0
        text_ypos -= offset
    endfor
   endif
   ;----------
   ; Close the plot file

   if (keyword_set(plotfile)) then $
    dfpsclose

   ;----------
   ; Compute and print statistics

   if keyword_set(tolog) then begin
       splog, 'Filename           ', $
        '  02%     05%     10%     50%     90%     95%     98%   '
       splog, '-------------------', $
        '  ------  ------  ------  ------  ------  ------  ------'
       for ifile=0, nfile-1 do $
        splog, fileandpath(fullname[ifile]), $
         rstruct[ifile].percentile[1], $
         rstruct[ifile].percentile[4], $
         rstruct[ifile].percentile[9], $
         rstruct[ifile].percentile[49], $
         rstruct[ifile].percentile[89], $
         rstruct[ifile].percentile[94], $
         rstruct[ifile].percentile[97], $
         format='(a19,7f8.1)'
   endif else begin
       print
       print, 'Filename           ', $
        '  02%     05%     10%     50%     90%     95%     98%   '
       print, '-------------------', $
        '  ------  ------  ------  ------  ------  ------  ------'
       for ifile=0, nfile-1 do $
        print, fileandpath(fullname[ifile]), $
         rstruct[ifile].percentile[1], $
         rstruct[ifile].percentile[4], $
         rstruct[ifile].percentile[9], $
         rstruct[ifile].percentile[49], $
         rstruct[ifile].percentile[89], $
         rstruct[ifile].percentile[94], $
         rstruct[ifile].percentile[97], $
         format='(a19,7f8.1)'
       print
   endelse
   return
end
;------------------------------------------------------------------------------
