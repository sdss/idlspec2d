; Plot the histograms of spectrophotometry vs. photo residuals
pro plot_sphoto1, mdiff, istd, igal, iqso, xtitle=xtitle, plotfile=plotfile

   dfpsplot, plotfile, /color, /square

   binsz = 0.002
   minval = -1.0
   maxval = 1.0
   csize = 1.7

   ; Plot standard stars
   djs_iterstat, mdiff[istd], mean=mn1, sigma=sig1
   binsz = 200. / n_elements(istd)
   xaxis = minval + binsz * findgen((maxval-minval)/binsz + 1)
   mhist = histogram(mdiff[istd], bin=binsz, min=minval, max=maxval)
   djs_plot, xaxis, mhist, psym=10, xtitle=xtitle, ytitle='Number', $
    /xstyle, charsize=csize, title='SPECTROPHOTOMETRY '+xtitle
   oplot, [0,0], !y.crange
   djs_xyouts, 0.9*minval, 0.9*!y.crange[1], 'Standard stars', charsize=csize
   djs_xyouts, 0.1*maxval, 0.9*!y.crange[1], charsize=csize, $
    string(mn1, sig1, format='(f6.3, " offset, ", f5.3," RMS")')

   ; Plot galaxies
   djs_iterstat, mdiff[igal], mean=mn1, sigma=sig1
   binsz = 200. / n_elements(igal)
   xaxis = minval + binsz * findgen((maxval-minval)/binsz + 1)
   mhist = histogram(mdiff[igal], bin=binsz, min=minval, max=maxval)
   djs_oplot, xaxis, mhist, psym=10, color='red'
   djs_xyouts, 0.9*minval, 0.8*!y.crange[1], 'Galaxies', color='red', charsize=csize
   djs_xyouts, 0.1*maxval, 0.8*!y.crange[1], color='red', charsize=csize, $
    string(mn1, sig1, format='(f6.3, " offset, ", f5.3," RMS")')

   ; Plot QSOs
   djs_iterstat, mdiff[iqso], mean=mn1, sigma=sig1
   binsz = 200. / n_elements(iqso)
   xaxis = minval + binsz * findgen((maxval-minval)/binsz + 1)
   mhist = histogram(mdiff[iqso], bin=binsz, min=minval, max=maxval)
   djs_oplot, xaxis, mhist, psym=10, color='blue'
   djs_xyouts, 0.9*minval, 0.7*!y.crange[1], 'QSOs', color='blue', charsize=csize
   djs_xyouts, 0.1*maxval, 0.7*!y.crange[1], color='blue', charsize=csize, $
    string(mn1, sig1, format='(f6.3, " offset, ", f5.3," RMS")')

   dfpsclose

   return
end

pro plot_sphoto

   spfile = filepath('spAll.fits', $
    root_dir='/u/dss/spectro_v5')
   minflux = 1.0

   columns=[ 'PLATE', 'MJD', 'FIBERID', $
    'CLASS', 'SUBCLASS', $
    'CALIBFLUX', 'SPECTROFLUX', $
    'PRIMTARGET', 'SECTARGET', 'Z', 'ZWARNING']

   objs = hogg_mrdfits(spfile, 1, $
    columns=columns, nrowchunk=10000L)

   qgood = total(objs.calibflux[1:3] GT minflux,1) EQ 3 $
    AND total(objs.spectroflux[1:3] GT minflux,1) EQ 3 $
    AND objs.zwarning EQ 0
   qstd = (objs.sectarget AND sdss_flagval('TTARGET','REDDEN_STD')) NE 0 $
    OR (objs.sectarget AND sdss_flagval('TTARGET','SPECTROPHOTO_STD')) NE 0
   qstd = qstd AND strmatch(objs.class,'STAR*')
   istd = where(qstd, nstd)

   qgal = (objs.primtarget AND sdss_flagval('TARGET','GALAXY')) NE 0 $
    OR (objs.primtarget AND sdss_flagval('TARGET','GALAXY_RED')) NE 0 $
    OR (objs.primtarget AND sdss_flagval('TARGET','GALAXY_RED_II')) NE 0
   qgal = qgal AND strmatch(objs.class,'GALAXY*')
   igal = where(qgal, ngal)

   qqso = (objs.primtarget AND sdss_flagval('TARGET','QSO_HIZ')) NE 0 $
    OR (objs.primtarget AND sdss_flagval('TARGET','QSO_CAP')) NE 0 $
    OR (objs.primtarget AND sdss_flagval('TARGET','QSO_SKIRT')) NE 0
   qqso = qqso AND strmatch(objs.class,'QSO*')
   iqso = where(qqso, nqso)

   smag = 22.5 - 2.5 * alog10(objs.spectroflux>minflux)
   pmag = 22.5 - 2.5 * alog10(objs.calibflux>minflux)
   magdiff = smag - pmag
   pcolor = pmag[0:3,*] - pmag[1:4,*]
   scolor = smag[0:3,*] - smag[1:4,*]
   colordiff = scolor - pcolor

   plot_sphoto1, transpose(magdiff[1,*]), istd, igal, iqso, $
    xtitle='g_{spectro} - g_{photo}', plotfile='sphoto_g.ps'
   plot_sphoto1, transpose(magdiff[2,*]), istd, igal, iqso, $
    xtitle='r_{spectro} - r_{photo}', plotfile='sphoto_r.ps'
   plot_sphoto1, transpose(magdiff[3,*]), istd, igal, iqso, $
    xtitle='i_{spectro} - i_{photo}', plotfile='sphoto_i.ps'
   plot_sphoto1, transpose(colordiff[1,*]), istd, igal, iqso, $
    xtitle='(g-r)_{spectro} - (g-r)_{photo}', plotfile='sphoto_gr.ps'
   plot_sphoto1, transpose(colordiff[2,*]), istd, igal, iqso, $
    xtitle='(r-i)_{spectro} - (r-i)_{photo}', plotfile='sphoto_ri.ps'

   return
end
