pro rline_plotsamp

   restore, 'rline_sample.ss'

   dfpsplot, 'rline_samp.ps', /color

   psize = 2.0
   pthick = 2

   xrange = [0.0, 0.6]
   yrange = [15.5,21.5]
   plot, [0], [0], /nodata, $
    position=[0.15,0.50,0.95,0.92], $
    ytitle='r-band mag', $
    xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
    charsize=psize, charthick=pthick, $
    xtickname=replicate(' ',30), $
    title='SDSS Luminous Red Galaxy Sample'

   oplot, sampzans.z, sampplug.mag[2], ps=3

   binsz = 0.01
   nbin = long( (xrange[1]-xrange[0]) / binsz ) + 1
   xhist = findgen(nbin) * binsz
   yhist = lonarr(nbin)
   for i=0, n_elements(sampzans.z)-1 do begin
      ibin = long(sampzans[i].z / binsz + 0.5)
      if (ibin GE 0 AND ibin LT nbin) then yhist[ibin] = yhist[ibin] + 1
   endfor

   yrange = [0,1.05*max(yhist)]
   plot, [0], [0], /nodata, /noerase, $
    position=[0.15,0.15,0.95,0.50], $
    xtitle='Redshift', ytitle='Number', $
    xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
    charsize=psize, charthick=pthick

;   oplot, xhist, yhist, psym=10, thick=2
   plothist, sampzans.z, bin=binsz, thick=2, /fill, /overplot

   djs_xyouts, 0.4, 0.8*!y.crange[1], $
    'N_{gal} = '+strtrim(string(n_elements(sampzans)),2), charsize=psize

   dfpsclose

   return
end
