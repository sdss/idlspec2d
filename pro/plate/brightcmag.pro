; Procedure for making color-magnitude plots for the four
; star cluster plates.  Do this for Stephane Charlot.

pro brightcmag1, plugfile, plottitle, plotfile

   plugfile = findfile('/peyton/scr/spectro0/data/platedb/plate-0797-0811/' $
    + plugfile)
   plugfile = plugfile[ sort(plugfile) ]
   ; Don't read the last pointing, which is only for the fundamental star
   plugfile = plugfile[0:n_elements(plugfile)-2]
   for ifile=0, n_elements(plugfile)-1 do begin
print, 'Reading ', plugfile[ifile]
      yanny_read, plugfile[ifile], pp
      thisplug = *pp[0]
      thisplug = thisplug[ where(thisplug.mag[1] NE 25 $
       AND strtrim(thisplug.holetype,2) EQ 'OBJECT') ]
print, 'Mag range = ', minmax(thisplug.mag[1])
      if (NOT keyword_set(plug)) then plug = thisplug $
       else plug = [plug, thisplug]
      yanny_free, pp
   endfor

   ; Re-derive the "Tycho magnitudes", which we had converted to SDSS
   ; magnitudes according to the prescription in sdss-spectro/766.

   b_tycho = 1.445990 * plug.mag[1] - 0.447249 * plug.mag[2] + 0.17953
   v_tycho = 0.467941 * plug.mag[1] + 0.533559 * plug.mag[2] - 0.03654

   ; Read the Tycho stars
   tycdat = tyc_read(/small, epoch=1998)
   racen = median(plug.ra)
   deccen = median(plug.dec)
   adiff = djs_diff_angle(tycdat.radeg, tycdat.dedeg, racen, deccen)
   tycdat = tycdat[ where(adiff LT 3.0) ]

   ; For every star in the plug-map, find the match in Tycho
   junk = djs_angle_match(plug.ra, plug.dec, $
    tycdat.radeg, tycdat.dedeg, dtheta=2.0/3600, mindx=mindx)
   igood = where(mindx NE -1, ngood)
   if (ngood GT 0) then begin
      v_tycho[igood] = tycdat[mindx[igood]].vmag
      b_tycho[igood] = tycdat[mindx[igood]].bmv + v_tycho[igood]
   endif

   dfpsplot, plotfile, /square
   bv_tycho = b_tycho - v_tycho
   plot, bv_tycho>(-1.8), b_tycho, psym=4, yrange=reverse(minmax(b_tycho)), $
    xtitle='(B-V)', ytitle='B', title=plottitle, charsize=2
   dfpsclose
end

pro brightcmag

   plugfile = 'plPlugMapP-0798*.par'
   plottitle = 'Praesepe'
   plotfile = 'cmag-praesepe.ps'
   brightcmag1, plugfile, plottitle, plotfile

   plugfile = 'plPlugMapP-0799*.par'
   plottitle = 'Pleiades'
   plotfile = 'cmag-pleiades.ps'
   brightcmag1, plugfile, plottitle, plotfile

   plugfile = 'plPlugMapP-0800*.par'
   plottitle = 'NGC 752'
   plotfile = 'cmag-n752.ps'
   brightcmag1, plugfile, plottitle, plotfile

   plugfile = 'plPlugMapP-0801*.par'
   plottitle = 'NGC 1817'
   plotfile = 'cmag-n1817.ps'
   brightcmag1, plugfile, plottitle, plotfile

end
