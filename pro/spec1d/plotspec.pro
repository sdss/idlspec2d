pro plotspec, plate, fiberid, mjd=mjd

   readspec, plate, fiberid, flux=objflux, flerr=objerr, loglam=loglam, $
    plug=plug, zans=zans
   wave = 10^loglam
   synflux = synthspec(zans, loglam=loglam)

   primtarget = sdss_flagname('TARGET', plug.primtarget, /concat)
   sectarget = sdss_flagname('TTARGET', plug.sectarget, /concat)

   csize = 2
   yrange = minmax(synflux)
   ymin = 1.2 * yrange[0] - 0.2 * yrange[1] < 0
   ymax = -0.2 * yrange[0] + 1.2 * yrange[1]
   if (ymax EQ ymin) then ymax = ymin + 1

   title = 'Plate ' + strtrim(string(plate),2) $
    + '  Fiber ' + strtrim(string(fiberid),2) $
    + '  MJD=' + strtrim(string(zans.mjd),2)
   splot, wave, objflux, yrange=[ymin,ymax], $
    xtitle='Wavelength [Ang]', $
    ytitle=TeXtoIDL('Flux [10^{-17} erg/cm/s/Ang]'), $
    title=title, charsize=csize
   soplot, wave, objerr, color='red'
   soplot, wave, synflux, color='blue', lw=2

   xpos = 0.9 * !x.crange[0] + 0.1 * !x.crange[1]
   dypos = 0.05 * (!y.crange[0] - !y.crange[1])
   ypos = !y.crange[1] + 1.5 * dypos
   sxyouts, xpos, ypos, zans.class + ' ' + zans.subclass $
    + '  z=' + strtrim(string(zans.z),2), $
    charsize=csize

   ypos = ypos + dypos
   sxyouts, xpos, ypos, $
    TeXtoIDL('X^2_r =' + strtrim(string(zans.chi2 / (zans.dof>1), format='(f5.2)'),2)), $
    charsize=csize

   if (keyword_set(primtarget)) then begin
      ypos = ypos + dypos
      sxyouts, xpos, ypos, 'PRIMTARGET = ' + primtarget, charsize=csize
   endif

   if (keyword_set(sectarget)) then begin
      ypos = ypos + dypos
      sxyouts, xpos, ypos, 'SECTARGET = ' + sectarget, charsize=csize
   endif

end
