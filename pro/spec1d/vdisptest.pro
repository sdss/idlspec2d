; Tests for looking at the velocity dispersion fits for plate 406.
; Loop from the highest-S/N galaxies to the lowest.
pro vdisptest

eigenfile='spEigenElodie.fits'
eigendir='.'
columns=lindgen(4)
   readspec, 406, mjd=52238, flux=objflux, invvar=objivar, zans=zans
   hdr = headfits('spPlate-0406-52238.fits')
   igal = where(strtrim(zans.class,2) EQ 'GALAXY' AND zans.zwarning EQ 0, $
    ngal)
   ; Sort by S/N
   igal = igal[ reverse(sort(zans[igal].sn_median)) ]
   loglam = sxpar(hdr,'COEFF0') + lindgen(sxpar(hdr,'NAXIS1')) $
    * sxpar(hdr,'COEFF1')

   for jj=0, ngal-1 do begin
      wave = 10^loglam / (1 + zans[igal[jj]].z)
      vdispfit, objflux[*,igal[jj]], objivar[*,igal[jj]], $
       hdr=hdr, zobj=zans[igal[jj]].z, $
       eigenfile=eigenfile, eigendir=eigendir, columns=columns, $
       sigma=sigma, sigerr=sigerr, yfit=yfit
      ipix = where(objivar[*,igal[jj]] GT 0 AND yfit NE 0) > 0
      ymax = 1.25 * max(median(objflux[ipix,igal[jj]],101))
      ymin = -0.2 * ymax
      splot, wave[ipix], objflux[ipix,igal[jj]], yrange=[ymin,ymax], $
       color='default', xtitle='Rest-frame Wavelength'
      soplot, wave[ipix], yfit[ipix], color='red'
      soplot, !x.crange, [0,0]
      soplot, wave[ipix], objflux[ipix,igal[jj]]-yfit[ipix], color='red'
print,'Fiber = ', zans[igal[jj]].fiberid, ' sigma=', sigma, ' +/- ', sigerr
stop
plothist, (objflux[ipix,igal[jj]]-yfit[ipix])*sqrt(objivar[ipix,igal[jj]]), $
 bin=0.1,xr=[-10,10]
wait,0.4
   endfor


return
end
