; Plot velocity errors from M67 plate (plate 321)
pro plotm67

   dfpsplot, 'm67vels.ps'

   cspeed = 3.e5
   csize = 2.0

   readspec, 321, zans=zans, plug=plug
   indx = where(plug.expl GT -900 and zans.zwarning EQ 0 $
    AND strmatch(zans.class,'STAR*'))
   zans = zans[indx]
   plug = plug[indx]

   cz_dave = zans.z * cspeed
   cz_err = zans.z_err * cspeed
   cz_cat = plug.expl
   vdiff = cz_dave - cz_cat

   plot, cz_cat, cz_dave, psym=4, charsize=csize, $
    xtitle='Catalog cz [km/s]', ytitle='SDSS cz [km/s]', $
    title='Velocities in M67 (Plate 321)'
   djs_oploterr, cz_cat, cz_dave, yerr=cz_err
   djs_oplot, !x.crange, !x.crange

   ibad = where(abs(vdiff) GT 20)
   djs_xyouts, cz_cat[ibad], cz_dave[ibad], $
    ' '+strtrim(zans.subclass,2), charsize=csize

   dfpsclose

   return
end

