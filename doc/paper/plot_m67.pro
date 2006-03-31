; Plot velocities from M67 compared to catalog velocities.
;------------------------------------------------------------------------------
pro plot_m67

   xrange = [-40,120]
   yrange = [-40,120]

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
    vel1[indx], vel2[indx], psym=4, charsize=2
   djs_oplot, !x.crange, !y.crange
   dfpsclose

   return
end
;------------------------------------------------------------------------------
