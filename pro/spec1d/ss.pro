pro ss, fiber

   readspec, 306, fiber, flux=objflux

;   hdr = headfits('/scr1/data/2d_test/0306new2/spPlate-0306-51690.fits')
;   zans = mrdfits('/scr1/data/2d_test/0306new2/spZ-0306-51690.fits',1)

   hdr = headfits('/home/data/2d_test/0306/spPlate-0306-51690.fits')
   zans = mrdfits('/home/data/2d_test/0306/spZ-0306-51690.fits',1)
   nper = n_elements(zans) / 640L
   zans = zans[nper*(fiber-1)]

   synflux = synthspec(zans, hdr=hdr)

   yrange = minmax(synflux)
   ymin = 1.2 * yrange[0] - 0.2 * yrange[1]
   ymax = -0.2 * yrange[0] + 1.2 * yrange[1]

   splot, objflux, yrange=[ymin,ymax]
   soplot, synflux, color='red', lw=2

   xpos = 0.9 * !x.crange[0] + 0.1 * !x.crange[1]
   ypos = 0.1 * !y.crange[0] + 0.9 * !y.crange[1]
   sxyouts, xpos, ypos, zans.class + ' ' + zans.subclass $
    + '  z=' + strtrim(string(zans.z),2), $
    charsize=3
end

