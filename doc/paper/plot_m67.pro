; Plot velocities from M67 compared to catalog velocities.
;------------------------------------------------------------------------------
pro plot_m67

   setenv, 'SPECTRO_DATA=/nfs/baryon8/data/spectro/2d_test'

   xrange = [-40,120]
   yrange = [-40,120]
   thick = 3
   csize = 2

   readspec, 321, mjd=51612, zans=zans, plug=plug

   cspeed = 2.99792458d5
   vel1 = plug.expl
   vel2 = zans.elodie_z * cspeed
   rmag = 22.5 - 2.5 * alog10(zans.spectroflux[2] > 1)
;   rmag = plug.mag[2]

   indx = where(strmatch(zans.class,'STAR*') AND zans.zwarning EQ 0 $
    AND plug.expl GT -900 AND plug.expl NE 0 AND rmag LT 20)

   dfpsplot, 'm67.ps', /square
   djs_plot, xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
    xtitle='Mathieu et al. velocity [km/s]', ytitle='SDSS velocity [km/s]', $
    vel1[indx], vel2[indx], psym=4, $
    charsize=csize, charthick=thick, thick=thick
   djs_oplot, !x.crange, !y.crange, $
    charsize=csize, charthick=thick, thick=thick
   dfpsclose

   return
end
;------------------------------------------------------------------------------
