; Make a plot of the number of spectra as a function of time (MJD)
pro platehist

   platelist, plist=plist

   ;  The following would trim to unique survey quality plates
;   plist = plist[where(plist.qsurvey)]

   ; The following trims to all survey quality data, including repeats
   minsn2 = 13.0
   plist = plist[where(plist.sn2_g1 GT minsn2 AND plist.sn2_g2 GT minsn2 $
    AND plist.sn2_i1 GT minsn2 AND plist.sn2_i2 GT minsn2)]

   nbin = max(plist.mjd)-min(plist.mjd)+1
   mjdvec = min(plist.mjd) + lindgen(nbin)
   totvec = lonarr(nbin)
   galvec = lonarr(nbin)
   qsovec = lonarr(nbin)
   starvec = lonarr(nbin)
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
   xyouts, mjdvec[nbin-1], unkvec[nbin-1], $
    string(unkvec[nbin-1], format='("Unclassified (",i6,")")'), $
    charsize=csize, align=0.5

   for i=1, n_elements(mjdlist)-1 do begin
      djs_oplot, [mjdlist[i],mjdlist[i]], !y.crange, linestyle=1
      xyouts, mjdlist[i]-18, total(!y.crange * [0.60,0.40]), $
       datelist[i], orient=90, charsize=csize, align=0.5
   endfor

   ; Jan 1, 2000 = MJD 51544
   ; Jul 1, 2000 = MJD 51726
   ; Jan 1, 2001 = MJD 51910
   ; Jul 1, 2001 = MJD 52092 ?
;   mjdlist = [51544, 51726, 51910, 52092]
;   datelist = ['1 Jan 2000', '1 July 2000', '1 Jan 2001', '1 July 2001']
;   for i=0, n_elements(mjdlist)-1 do begin
;      djs_oplot, [mjdlist[i],mjdlist[i]], !y.crange, linestyle=1
;      xyouts, mjdlist[i]-18, total(!y.crange * [0.60,0.40]), $
;       datelist[i], orient=90, charsize=csize, align=0.5
;   endfor

   dfpsclose

   return
end
