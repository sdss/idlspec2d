pro plotspec, plate, fiberid, mjd=mjd

   readspec, plate, fiberid, flux=objflux, flerr=objerr, loglam=loglam, $
    zans=zans
   wave = 10^loglam
   synflux = synthspec(zans, loglam=loglam)

   csize = 2
   yrange = minmax(synflux)
   ymin = 1.2 * yrange[0] - 0.2 * yrange[1] < 0
   ymax = -0.2 * yrange[0] + 1.2 * yrange[1]

   title = 'Plate ' + strtrim(string(plate),2) $
    + '  Fiber ' + strtrim(string(fiberid),2) $
    + '  MJD=' + strtrim(string(zans.mjd),2)
   splot, wave, objflux, yrange=[ymin,ymax], title=title, charsize=csize
   soplot, wave, objerr, color='red'
   soplot, wave, synflux, color='blue', lw=2

   xpos = 0.9 * !x.crange[0] + 0.1 * !x.crange[1]
   ypos = 0.1 * !y.crange[0] + 0.9 * !y.crange[1]
   sxyouts, xpos, ypos, zans.class + ' ' + zans.subclass $
    + '  z=' + strtrim(string(zans.z),2), $
    charsize=csize

   ypos = 0.2 * !y.crange[0] + 0.8 * !y.crange[1]
   sxyouts, xpos, ypos, $
    ' X^2_r =' + strtrim(string(zans.chi2 / (zans.dof>1), format='(f5.2)'),2), $
    charsize=csize
end

