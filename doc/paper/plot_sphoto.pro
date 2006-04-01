; Plot the histograms of spectrophotometry vs. photo residuals
;------------------------------------------------------------------------------
pro plot_sphoto1, mdiff, istd, igal, iqso, xtitle=xtitle, plotfile=plotfile

   dfpsplot, plotfile, /color, /square

   binsz = 0.002
   minval = -1.0
   maxval = 1.0
   csize = 1.5

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
   splog, 'Stars ' + xtitle, ' mean=', mn1, ' sigma=', sig1

   ; Plot galaxies
   djs_iterstat, mdiff[igal], mean=mn1, sigma=sig1
   binsz = 200. / n_elements(igal)
   xaxis = minval + binsz * findgen((maxval-minval)/binsz + 1)
   mhist = histogram(mdiff[igal], bin=binsz, min=minval, max=maxval)
   djs_oplot, xaxis, mhist, psym=10, color='red'
   djs_xyouts, 0.9*minval, 0.8*!y.crange[1], 'Galaxies', color='red', charsize=csize
   djs_xyouts, 0.1*maxval, 0.8*!y.crange[1], color='red', charsize=csize, $
    string(mn1, sig1, format='(f6.3, " offset, ", f5.3," RMS")')
   splog, 'Galaxies ' + xtitle, ' mean=', mn1, ' sigma=', sig1

   ; Plot QSOs
   djs_iterstat, mdiff[iqso], mean=mn1, sigma=sig1
   binsz = 200. / n_elements(iqso)
   xaxis = minval + binsz * findgen((maxval-minval)/binsz + 1)
   mhist = histogram(mdiff[iqso], bin=binsz, min=minval, max=maxval)
   djs_oplot, xaxis, mhist, psym=10, color='blue'
   djs_xyouts, 0.9*minval, 0.7*!y.crange[1], 'QSOs', color='blue', charsize=csize
   djs_xyouts, 0.1*maxval, 0.7*!y.crange[1], color='blue', charsize=csize, $
    string(mn1, sig1, format='(f6.3, " offset, ", f5.3," RMS")')
   splog, 'QSOs ' + xtitle, ' mean=', mn1, ' sigma=', sig1

   dfpsclose

   return
end

;------------------------------------------------------------------------------
pro plot_sphoto, synflux=synflux

   minflux = 1.0

   ;----------
   ; Get the list of EDR-DR5 plates

   public = ['EDR','DR1','DR2','DR3','DR4','DR5']
   platelist, plist=plist
   qkeep = bytarr(n_elements(plist))
   for i=0, n_elements(public)-1 do $ 
    qkeep = qkeep OR strmatch(plist.public,public[i])
   qkeep = qkeep AND strmatch(plist.statuscombine,'Done*')
   plist = plist[where(qkeep, nplate)]

   columns=[ 'PLATE', 'MJD', 'FIBERID', $
    'CLASS', 'SUBCLASS', $
    'CALIBFLUX', 'SPECTROFLUX', 'SPECTROSYNFLUX', $
    'PRIMTARGET', 'SECTARGET', 'Z', 'ZWARNING']

   for iplate=0L, nplate-1L do begin
      print, format='("Plate ",i4," of ",i4,a1,$)', iplate, nplate, string(13b)
      readspec, plist[iplate].plate, mjd=plist[iplate].mjd, $
       zans=zans1, plug=plug1, /silent
      zans1 = struct_selecttags(zans1, select_tags=columns)
      plug1 = struct_selecttags(plug1, select_tags=columns)
      plug1 = struct_selecttags(plug1, except_tags=tag_names(zans1))
      obj1 = struct_addtags(zans1, plug1)
      if (iplate EQ 0) then objs = replicate(obj1[0], nplate*640L)
      copy_struct_inx, obj1, objs, index_to=640L*iplate+lindgen(640)
   endfor

   if (keyword_set(synflux)) then objs.spectroflux = objs.spectrosynflux

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
   plot_sphoto1, transpose(colordiff[1,*]-colordiff[2,*]), istd, igal, iqso, $
    xtitle='(g-i)_{spectro} - (g-i)_{photo}', plotfile='sphoto_gi.ps'

   return
end
;------------------------------------------------------------------------------
