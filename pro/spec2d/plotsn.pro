;+
; NAME:
;   plotsn
;
; PURPOSE:
;   Plot S/N and residuals in up to 3 bands
;
; CALLING SEQUENCE:
;   plotsn, sn, plugmap, bands=bands, title=title
;
; INPUTS:
;   sn         - 2d S/N array [nbands, nfibers]
;   plugmap    - plugmap structure [nfibers]
;
; OPTIONAL KEYWORDS:
;   bands      - Label of each band 
;   title      - Add title to top of plot
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
pro plotsn, sn, plug, bands=bands, title=title

   if (NOT keyword_set(bands)) then bands=[1, 2, 3]

   bandnames = ["* u'","* g'","* r'","* i'","* z'"]
   snmag = [19.0, 19.0, 19.0, 19.0, 19.0]
   snof1 = [60.0, 50.0, 50.0, 100.0, 60.0]/100.0

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
         xyouts, minmag+0.5, 110.0, title, charsize=2.0

      if (s1[0] NE -1) then djs_oplot, mag[s1], snc[s1] > 0.9, $
           ps=1, syms=0.7, color='red'
      if (s2[0] NE -1) then djs_oplot, mag[s2], snc[s2] > 0.9, ps=4, syms=0.7

      xyouts, minmag+0.5, 3.0, string(format='(a,f5.2,f6.2,a)','log S/N = ', $
           a, bandnames[bands[i]])

      if (size(sigma, /tname) NE 'UNDEFINED') then $
      xyouts, minmag+0.5, 2.0, string(format='(a,f5.3)','Deviation: ', sigma)

      djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.6], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.91], ps=1, $
         syms=0.7, color='red'
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.62], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.9], $
         'Spec 1', color=djs_icolor('red')

      djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.6], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.85], ps=4, syms=0.7
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.62], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.84], 'Spec 2'

      minutes = fltarr(2)

      if (s1[0] NE -1) then begin
        a1 = fitsn(mag[s1], snc[s1], minmag=snmag[bands[i]]-0.5, $
                  maxmag=snmag[bands[i]]+0.5)
        minutes[0] = 10^(2.0 * poly(snmag[bands[i]],a1)) * snof1[bands[i]]
      endif
      if (s2[0] NE -1) then begin
        a2 = fitsn(mag[s2], snc[s2], minmag=snmag[bands[i]]-0.5, $
                  maxmag=snmag[bands[i]]+0.5)
        minutes[1] = 10^(2.0 * poly(snmag[bands[i]],a2)) * snof1[bands[i]]
      endif

      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.8], $
           10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.9], $
           string(format='(f5.1)', minutes[0]), color=djs_icolor('red')
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.8], $
           10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.84], $
           string(format='(f5.1)', minutes[1])


      good = where(snc GT 0.0 AND mag GT minmag, ngood)

      if (ngood GT 0) then $
        s_plotdev, plugc[good].xfocal, plugc[good].yfocal, $
            diff[good], 3.0, ychars=2.0, xcharsize = xchars, $
            ymargin=ymargin, cap=1.0 $
        else plot, plug.xfocal, plug.yfocal, ychars=2.0, xchars=xchars, $
            ymargin=ymargin, /nodata
      if (nlow GT 0) then djs_oplot, plugc[low].xfocal, plug[low].yfocal, ps=1


  endfor

  !p.multi  = oldmulti

return
end
