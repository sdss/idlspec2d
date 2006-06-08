; Plot one of the best throughput curves, for example plate 406/51817
; exposure #6791
;------------------------------------------------------------------------------
pro plot_bestthru

   plate = 406
   mjd = 51817
   expnum = 6791
   plotfile = 'plot_bestthru.ps'

   ;----------

   npix = 4000
   loglam = 3.5700 + dindgen(npix)*1d-4
   eff = fltarr(npix,4)
   camname = ['b1','b2','r1','r2']
   for icam=0, 3 do begin
      foo = spthroughput(plate, camname=camname[icam], expnum=expnum, $
       loglam=loglam, efficiency=eff1, /median)
      eff[*,icam] = eff1
   endfor

   dfpsplot, plotfile, /color, /encap, /landscape
   djs_plot, [0], [0], /nodata, xrange=[3700,9300], yrange=[0,0.22], $
    /xstyle, /ystyle, xtitle=textoidl('Wavelength [\AA]'), $
    ytitle='Throughput', charsize=1.4
   for icam=0, 3 do begin
      ii = where(eff[*,icam] GT 0)
      thiscolor = strmatch(camname[icam],'b*') ? 'blue' : 'red'
      lstyle = strmatch(camname[icam],'*1') ? 0 : 1
      djs_oplot, 10^loglam[ii], eff[ii,icam], color=thiscolor, $
       linestyle=lstyle
   endfor
   dfpsclose

   return
end
;------------------------------------------------------------------------------
