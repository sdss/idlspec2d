; Tests for looking at the velocity dispersion fits for plate 406.
; Loop from the highest-S/N galaxies to the lowest.
pro vdisptest

eigenfile='spEigenElodie.fits'
columns=lindgen(24)
debug = 1
brightsort = 0

   if (keyword_set(debug)) then doplot = 1 $
    else doplot = 0

   readspec, 406, mjd=52238, flux=objflux, invvar=objivar, zans=zans, $
    loglam=loglam
   loglam = loglam[*,0] ; same for every object
   if (keyword_set(brightsort)) then begin
      igal = where(strtrim(zans.class,2) EQ 'GALAXY' AND zans.zwarning EQ 0, $
       ngal)
      igal = igal[ reverse(sort(zans[igal].sn_median)) ] ; Sort by S/N
   endif else begin
      ngal = (n_elements(zans))
      igal = lindgen(ngal)
   endelse

   for jj=0, ngal-1 do begin
      wave = 10^loglam / (1 + zans[igal[jj]].z)
      vdans1 = vdispfit(objflux[*,igal[jj]], objivar[*,igal[jj]], $
       loglam, zobj=zans[igal[jj]].z, $
       eigenfile=eigenfile, eigendir=eigendir, columns=columns, $
       yfit=yfit1, /doplot, /debug)
       if (jj EQ 0) then begin
          vdans = vdans1
          yfit = yfit1
       endif else begin
          vdans = [[vdans],[vdans1]]
          yfit = [[yfit],[yfit1]]
       endelse

      if (keyword_set(doplot)) then begin
         ipix = where(objivar[*,igal[jj]] GT 0 AND yfit1 NE 0) > 0
         ymax = 1.25 * max(median(objflux[ipix,igal[jj]],101))
         ymin = -0.2 * ymax
         djs_plot, wave[ipix], objflux[ipix,igal[jj]], yrange=[ymin,ymax], $
          color='default', xtitle='Rest-frame Wavelength'
         djs_oplot, wave[ipix], yfit1[ipix], color='red'
         djs_oplot, !x.crange, [0,0]
         djs_oplot, wave[ipix], objflux[ipix,igal[jj]]-yfit1[ipix], color='red'
print,'Fiber = ', zans[igal[jj]].fiberid, ' sigma=', vdans1.vdisp, $
 ' +/- ', vdans1.vdisp_err

         if (keyword_set(debug)) then begin
            print, 'Press any key...'
            cc = strupcase(get_kbrd(1))
         endif

         chivec = (objflux[ipix,igal[jj]]-yfit1[ipix]) $
          * sqrt(objivar[ipix,igal[jj]])
         plothist, chivec, bin=0.1,xr=[-10,10]
print, 'Median chi=',median(abs(chivec))

         if (keyword_set(debug)) then begin
            print, 'Press any key...'
            cc = strupcase(get_kbrd(1))
         endif
      endif
   endfor

return
end
