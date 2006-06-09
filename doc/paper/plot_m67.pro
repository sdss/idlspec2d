; Plot velocities from M67 compared to catalog velocities.
;------------------------------------------------------------------------------
pro plot_m67

;   setenv, 'SPECTRO_DATA=/nfs/baryon8/data/spectro/2d_test'
   setenv, 'SPECTRO_DATA=/nfs/baryon8/data/spectro/2d_120sec'

   xrange = [-30,60]
   yrange = [-25,25]
   thick = 3
   csize = 2

   readspec, 321, mjd=51612, zans=zans, plug=plug

   cspeed = 2.99792458d5
   vel1 = plug.expl
   vel2 = zans.elodie_z * cspeed
   rmag = 22.5 - 2.5 * alog10(zans.spectroflux[2] > 1)
;   rmag = plug.mag[2]

   indx = where(strmatch(zans.class,'STAR*') $
;    AND zans.zwarning EQ 0 $
    AND plug.expl GT -900 AND plug.expl NE 0 AND rmag LT 20)

   dfpsplot, 'plot_m67.ps', /square
   djs_plot, xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
    xtitle='Mathieu et al. velocity [km/s]', $
    ytitle='SDSS - Mathieu velocity [km/s]', $
    vel1[indx], vel2[indx]-vel1[indx], psym=4, $
    charsize=csize, charthick=thick, thick=thick, xthick=thick, ythick=thick
   djs_oplot, !x.crange, [0,0], charsize=csize, charthick=thick, thick=thick
   dfpsclose

   splog, 'Number of stars plotted = ', n_elements(indx)
   splog, 'Min/max velocity difference = ', minmax(vel1[indx] - vel2[indx])
   splog, 'Median velocity difference = ', median(vel1[indx] - vel2[indx])
   splog, 'Stdev velocity difference = ', djsig(vel1[indx] - vel2[indx])

   ibad = where(abs(vel1[indx] - vel2[indx]) GT 30, nbad)
   for i=0, nbad-1 do $
    splog, 'Fiber ', plug[indx[ibad[i]]].fiberid, $
     ' at RA=', plug[indx[ibad[i]]].ra, $
     ' Dec=', plug[indx[ibad[i]]].dec, ' v(CAT)=', vel1[indx[ibad[i]]], $
     ' v(Elodie)=', vel2[indx[ibad[i]]]

   return
end
;------------------------------------------------------------------------------
