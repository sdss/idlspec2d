;+
; NAME:
;   plotsn
;
; PURPOSE:
;   Plot S/N and residuals in up to 3 bands
;
; CALLING SEQUENCE:
;   plotsn, sn, plugmap, [ bands=, title=, plotfile= ]
;
; INPUTS:
;   sn         - 2d S/N array [nbands, nfibers]
;   plugmap    - plugmap structure [nfibers]
;
; OPTIONAL KEYWORDS:
;   bands      - Label of each band 
;   title      - Add title to top of plot
;   plotfile   - Optional plot file
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
;   s_plotdev
;
; REVISION HISTORY:
;   15-Apr-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
pro plotsn, sn, plug, bands=bands, title=title, plotfile=plotfile

   if (NOT keyword_set(bands)) then bands=[1, 2, 3]

   ;----------
   ; Close plot file

   if (keyword_set(plotfile)) then begin
      set_plot, 'ps'
      device, filename=plotfile, /color
   endif

   bandnames = ["u'","g'","r'","i'","z'"]
   slopelabel = " * "+bandnames
   snmag = [19.0, 19.2, 19.25, 18.9, 18.0]
   snlabel = '(S/N)^2 @ '+bandnames+' ='+string(snmag,format='(f6.2)')

   nbands = n_elements(bands)
   nfibers = n_elements(plug)

   if (size(sn,/n_dim) NE 2) then return

   if ((size(sn,/dimen))[0] NE nbands) then return
   if ((size(sn,/dimen))[1] NE nfibers) then return

   nonsky = where(strtrim(plug.objtype,2) NE 'SKY', nnonsky)
   if (nnonsky LT 3) then return

   plugc = plug[nonsky]

   bright = bytarr(nnonsky) + 1

   spectroid = (plugc.fiberid - 1) / 320

   s1 = where(spectroid EQ 0)
   s2 = where(spectroid EQ 1)

   oldmulti = !p.multi  
   !p.multi = [0,2,nbands]
   !p.charsize = 1.6

   minmag = 15.0
   maxmag = 21.0
   xfit = [minmag, maxmag]

   fiducialfits = [[6.12, -0.28], $    ; u' fit
                   [6.12, -0.28], $    ; g' fit
                   [6.18, -0.28], $    ; r' fit
                   [6.05, -0.28], $    ; i' fit
                   [6.05, -0.28]]      ; z' fit

   for i=0, nbands-1 do begin

      ; Let's first plot S/N vs magnitude

      snc = sn[i,nonsky]
      mag = plugc.mag[bands[i]]
      a = fitsn(mag, snc, sigma=sigma, minmag=minmag, maxmag=maxmag)

      logsnc = alog10(snc > 0.01)
      diff = logsnc - poly(mag, a)

      limit = 10^poly(mag, fiducialfits[*,bands[i]])
      low = where(snc LT limit AND mag LT maxmag, nlow)

      bright = bright * (diff GT 0.2)

      xchars=0.001
      ymargin=[1,0]
      if (i EQ nbands-1) then begin 
        xchars=2.0
        ymargin=[4,0]
      endif

      plot, mag, snc, /nodata, /ylog, $
           ychars=2.0, xchars =xchars, xr=[minmag,maxmag], $
           ymargin = ymargin, /xs, yr=[0.8,100], /ys

      djs_oplot, xfit, 10^poly(xfit, a)
      djs_oplot, xfit, 10^poly(xfit, fiducialfits[*,bands[i]]), $
           color='blue', thick=5
      if (i EQ 0 AND keyword_set(title)) then $
         xyouts, minmag+0.5, 110.0, title, charsize=1.2

      if (s1[0] NE -1) then begin
          upper = where(diff[s1] GT 0)
          if (upper[0] NE -1) then djs_oplot, mag[s1[upper]], $
             snc[s1[upper]] > 0.9, ps=1, syms=0.7, color='green'
          lower = where(diff[s1] LE 0)
          if (lower[0] NE -1) then djs_oplot, mag[s1[lower]], $
             snc[s1[lower]] > 0.9, ps=1, syms=0.7, color='red'
      endif
      if (s2[0] NE -1) then begin
          upper = where(diff[s2] GT 0)
          if (upper[0] NE -1) then djs_oplot, mag[s2[upper]], $
             snc[s2[upper]] > 0.9, ps=4, syms=0.7, color='green'
          lower = where(diff[s2] LE 0)
          if (lower[0] NE -1) then djs_oplot, mag[s2[lower]], $
             snc[s2[lower]] > 0.9, ps=4, syms=0.7, color='red'
      endif

      xyouts, minmag+0.5, 3.0, string(format='(a,f5.2,f6.2,a)','log S/N = ', $
           a, slopelabel[bands[i]])

      if (size(sigma, /tname) NE 'UNDEFINED') then $
      xyouts, minmag+0.5, 2.0, string(format='(a,f4.2)','Stdev=', sigma)

      djs_xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.6], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.94], $
         snlabel[bands[i]]

      djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.6], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.89], ps=1, $
         syms=0.7, color='red'
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.62], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.88], $
         'Spec 1', color=djs_icolor('red')

      djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.6], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.83], ps=4, syms=0.7
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.62], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.82], 'Spec 2'

      sn2 = fltarr(2)

      if (s1[0] NE -1) then begin
        a1 = fitsn(mag[s1], snc[s1], minmag=snmag[bands[i]]-0.5, $
                  maxmag=snmag[bands[i]]+0.5)
        sn2[0] = 10^(2.0 * poly(snmag[bands[i]],a1)) 
      endif
      if (s2[0] NE -1) then begin
        a2 = fitsn(mag[s2], snc[s2], minmag=snmag[bands[i]]-0.5, $
                  maxmag=snmag[bands[i]]+0.5)
        sn2[1] = 10^(2.0 * poly(snmag[bands[i]],a2)) 
      endif

      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.8], $
           10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.88], $
           string(format='(f5.1)', sn2[0]), color=djs_icolor('red')
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.8], $
           10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.82], $
           string(format='(f5.1)', sn2[1])

      print, snlabel[bands[i]], sn2

      good = where(snc GT 0.0 AND mag GT minmag, ngood)

      if (ngood GT 0) then $
        s_plotdev, plugc[good].xfocal, plugc[good].yfocal, $
            diff[good], 3.0, ychars=2.0, xcharsize = xchars, $
            ymargin=ymargin, cap=1.0 $
        else plot, plug.xfocal, plug.yfocal, ychars=2.0, xchars=xchars, $
            ymargin=ymargin, /nodata
      if (nlow GT 0) then djs_oplot, plugc[low].xfocal, plugc[low].yfocal, ps=1


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
