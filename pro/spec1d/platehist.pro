; Make a plot of the number of spectra as a function of time (MJD)
pro platehist

   platelist, plist=plist

   ; The following trims to all survey quality data, including repeats
   plist = plist[ where(strmatch(plist.platequality,'good*') $
    OR strmatch(plist.platequality,'marginal*')) ]

   nbin = max(plist.mjd)-min(plist.mjd)+1
   mjdvec = min(plist.mjd) + lindgen(nbin)
   totvec = lonarr(nbin)
   galvec = lonarr(nbin)
   qsovec = lonarr(nbin)
   starvec = lonarr(nbin)
   skyvec = lonarr(nbin)
   unkvec = lonarr(nbin)
;   for i=0, nbin-1 do $
;    totvec[i] = total(640 - plist[where(plist.mjd LE mjdvec[i])].n_sky)
   for i=0, nbin-1 do $
    totvec[i] = total(640 * n_elements(where(plist.mjd LE mjdvec[i])))
   for i=0, nbin-1 do $
    galvec[i] = total(plist[where(plist.mjd LE mjdvec[i])].n_galaxy)
   for i=0, nbin-1 do $
    qsovec[i] = total(plist[where(plist.mjd LE mjdvec[i])].n_qso)
   for i=0, nbin-1 do $
    starvec[i] = total(plist[where(plist.mjd LE mjdvec[i])].n_star)
   for i=0, nbin-1 do $
    skyvec[i] = total(plist[where(plist.mjd LE mjdvec[i])].n_sky)
   for i=0, nbin-1 do $
    unkvec[i] = total(plist[where(plist.mjd LE mjdvec[i])].n_unknown)

   mjd2datelist, min(mjdvec), max(mjdvec), step='year', $
    mjdlist=mjdlist, datelist=datelist

   csize = 1.6

   dfpsplot, 'platehist.ps', /color, /square
   djs_plot, minmax(mjdlist), minmax(totvec), /nodata, charsize=csize, $
    xtickformat='(i10)', /xstyle, $
    xtitle='Modified Julian Date', ytitle='Cumulative Number', $
    title='SDSS Survey Quality Spectra'
   djs_oplot, mjdvec, totvec, psym=10
   djs_oplot, mjdvec, galvec, psym=10, color='red'
   djs_oplot, mjdvec, qsovec, psym=10, color='green'
   djs_oplot, mjdvec, starvec, psym=10, color='blue'
   djs_oplot, mjdvec, skyvec, psym=10, color='magenta'
   djs_oplot, mjdvec, unkvec, psym=10, color='yellow'

   xyouts, mjdvec[nbin-1], totvec[nbin-1], $
    string(totvec[nbin-1], format='("Total (",i6,")")'), $
    charsize=csize, align=0.5
   xyouts, mjdvec[nbin-1], galvec[nbin-1], $
    string(galvec[nbin-1], format='("Galaxies (",i6,")")'), $
    charsize=csize, align=0.5
   xyouts, mjdvec[nbin-1], qsovec[nbin-1], $
    string(qsovec[nbin-1], format='("QSOs (",i6,")")'), $
    charsize=csize, align=0.5
   xyouts, mjdvec[nbin-1], starvec[nbin-1], $
    string(starvec[nbin-1], format='("Stars (",i6,")")'), $
    charsize=csize, align=0.5
   xyouts, mjdvec[nbin-1], skyvec[nbin-1], $
    string(skyvec[nbin-1], format='("Sky (",i6,")")'), $
    charsize=csize, align=0.5
   xyouts, mjdvec[nbin-1], unkvec[nbin-1], $
    string(unkvec[nbin-1], format='("Unclassified (",i6,")")'), $
    charsize=csize, align=0.5

   for i=0, n_elements(mjdlist)-1 do begin
      if (i NE 0 AND i NE n_elements(mjdlist)-1) then $
       djs_oplot, [mjdlist[i],mjdlist[i]], !y.crange, linestyle=1
      xoff = 0.02 * (!x.crange[1] - !x.crange[0])
      xplot = mjdlist[i] + 2*xoff
      xyouts, xplot, total(!y.crange * [0.60,0.40]), $
       datelist[i], orient=90, charsize=csize, align=0.5
   endfor

   dfpsclose

   return
end
