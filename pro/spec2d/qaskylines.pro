pro qaskylines, flux, fluxivar, wset, plugsort

        nTrace = (size(flux))[2]
;
;	Here are guesses at strong skylines in red 
;       (vacuum wavelengths)
;
	possible = [5891.8, 5897.8, 6866.2, 7278.6, 7343.2, $
                    7753.0, 7796.5, 7824.0, 7916.0, 7967.1, $
                    7995.6, 8028.2, 8401.6, 8432.5, 8829.6, $
                    8888.2, 8922.0, 8960.6]

	logpos = alog10(possible)

	pix = traceset2pix(wset, logpos)
        ycen = long(pix)
	for i=0,nTrace - 1 do ycen[*,i] = i
	
        skycen = trace_fweight(flux, pix, ycen, $
                   radius=2.0, invvar=fluxivar, xerr=skyerr)

	traceset2xy, wset, skycen, skywave

	skymed = djs_median(skywave,2)

        veldiff = (skywave-skymed # (fltarr(nTrace) +1.0)) * 3.0e5 * 2.302585

        djs_plot, possible # (fltarr(nTrace) +1.0), veldiff, ps=3, $
             yr=[-100,100], xtitle='Wavelength of sky lines', $
             ytitle='Offsets from median (km/s)'

        skies = where(plugsort.objtype EQ 'SKY',nskies)
	if (skies[0] NE -1) then $
          djs_oplot, possible # (fltarr(nskies) +1.0), veldiff[*,skies], $
               ps=4, color='red'

	return
end
