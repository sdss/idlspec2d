
; Compute the refraction terms for QSOs.

pro qsorefract, ztab=ztab, urefract=urefract, grefract=grefract, $
 ugcolor=ug, grcolor=gr, ricolor=ri, izcolor=iz

   eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
   plotfile = 'qso-refract.ps'

   ;----------
   ; Find the template QSO spectrum file

;   eigenfile = 'spEigenQSO-*.fits'
   eigenfile = 'spEigenQSO.fits'
   allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
   thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]
   splog, 'Selecting EIGENFILE=' + thisfile

   ;----------
   ; Read the QSO template spectrum

   qsoflux = readfits(thisfile, qhdr)
   qsoflux = smooth(qsoflux[*,0],25) ; Just take the mean spectrum (smoothed)
   naxis1 = sxpar(qhdr, 'NAXIS1')
   loglam0 = sxpar(qhdr, 'COEFF0')
   dloglam0 = sxpar(qhdr, 'COEFF1')
   qsowave = 10.d0^(loglam0 + lindgen(naxis1) * dloglam0)

   ;----------
   ; Find the template stellar spectrum file

   eigenfile = 'spEigenStar-*.fits'
   allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
   thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]
   splog, 'Selecting EIGENFILE=' + thisfile

   ;----------
   ; Read the stellar template spectrum

;   starflux = readfits(thisfile, shdr)
;   naxis1 = sxpar(shdr, 'NAXIS1')
;   loglam0 = sxpar(shdr, 'COEFF0')
;   dloglam0 = sxpar(shdr, 'COEFF1')
;   starwave = 10.d0^(loglam0 + lindgen(naxis1) * dloglam0)
;   starname = sxpar(shdr,'NAME0')
;   sxdelpar, shdr, 'NAME0'
;   starname = [starname, sxpar(shdr,'NAME*')]
;   nstar = n_elements(starname)

   ;----------
   ; Tabulate the index of refraction of air
   ; from Meggers & Peters 1919, ApJ 50, 56.

   rwave = 2000. + findgen(33) * 250.
   rval = [3255.82,3107.87,3014.05,2951.08,2906.85,2874.63,2850.43, $
    2831.79,2817.12,2805.36,2795.78,2787.88,2781.27,2775.69,2770.94, $
    2766.85,2763.31,2760.22,2757.51,2755.11,2752.99,2751.10,2749.40, $
    2747.88,2746.50,2745.25,2744.12,2743.09,2742.14,2741.28,2740.48, $
    2739.75,2739.07] * 1.0d-7

   ;----------
   ; Measure refraction terms for QSOs

   ztab = 0.0 + 0.05 * findgen(80)
   nz = n_elements(ztab)
   rtabqso = fltarr(nz,5)
   qsomag = fltarr(nz,5)

   for iz=0, nz-1 do begin
      ; Redshift the QSO spectrum
      tmpwave = qsowave * (1.0 + ztab[iz])
      tmpflux = qsoflux

      ; Apply Burles' Lyman forest opacity estimation
      bluewardlya = (qsowave LT 1215.67)
      bluewardlymanlimit = (920.0 - qsowave) > 0
      z_temp = 2.5
      zratio = (1 + ztab[iz]) / (1 + z_temp)
      lyacorrection = 0.10 * (zratio^3.0 - 1.0) * bluewardlya $
                   + 0.02 * (zratio^2.0 - 1.0) * bluewardlymanlimit
      tmpflux = tmpflux * exp(-lyacorrection)

      ; Interpolate the index of refraction curve
      tmprval = interpol(rval, rwave, tmpwave)

      ; Integrate over the SDSS filters -- counts photons, not energy
      flambda2fnu = tmpwave^2 / 2.99792e18
      num1 = filter_thru(tmpflux * tmprval * flambda2fnu, $
       waveimg=tmpwave, /toair)
      num2 = filter_thru(tmpflux * flambda2fnu, $
       waveimg=tmpwave, /toair)
      rtabqso[iz,*] = num1 / num2
      qsomag[iz,*] = num2
   endfor

   ;----------
   ; Measure refraction terms for stars

;   rtabstar = fltarr(nstar,5)
;   tmprval = interpol(rval, rwave, starwave)
;   for istar=0, nstar-1 do begin
;      num1 = filter_thru(starflux[*,istar] * tmprval, waveimg=starwave, /toair)
;      num2 = filter_thru(starflux[*,istar], waveimg=starwave, /toair)
;      rtabstar[istar,*] = num1 / num2
;   endfor

   ;----------
   ; Now make the dtheta-dtheta plot

   csize = 1.6
   scale = !radeg * 3600. ; Re-scale plots to arcsec at 45 deg from zenith

   if (keyword_set(plotfile)) then dfpsplot, plotfile, /square

   ; Compare against ave. r+i position
   rnorm = 0.5 * (rtabqso[*,2] + rtabqso[*,3])
   xqso = (rtabqso[*,0] - rnorm) * scale
   yqso = (rtabqso[*,1] - rnorm) * scale
   djs_plot, xqso, yqso, psym=-4, charsize=csize, /ynozero, $
    title='QSO Refraction-Refraction Plot', $
    xtitle='u-band refraction [arcsec]', ytitle='g-band refraction [arcsec]'
   for j=0, nz-1, 10 do $
    xyouts, xqso[j], yqso[j], 'z='+string(ztab[j],format='(f3.1)'), $
    charsize=csize

;   ii = [37,35,0,2,8,9,11,13,14,16,27,31,34]
;   rnorm = 0.5 * (rtabstar[*,2] + rtabstar[*,3])
;   xstar = (rtabstar[*,0] - rnorm) * scale
;   ystar = (rtabstar[*,1] - rnorm) * scale
;   for i=0, n_elements(ii)-1 do $
;    djs_xyouts, xstar[ii[i]], ystar[ii[i]], strtrim(starname[ii[i]],2)

   ug = -2.5 * alog10(qsomag[*,0] / qsomag[*,1])
   gr = -2.5 * alog10(qsomag[*,1] / qsomag[*,2])
   ri = -2.5 * alog10(qsomag[*,2] / qsomag[*,3])
   iz = -2.5 * alog10(qsomag[*,3] / qsomag[*,4])

   djs_plot, ug, gr, psym=-4, charsize=csize, /ynozero, $
    title='QSO Color-Color Plot', $
    xtitle='(u-g) [mag]', ytitle='(g-r) [mag]'
   for j=0, nz-1, 10 do $
    xyouts, ug[j], gr[j], 'z='+string(ztab[j],format='(f3.1)'), $
    charsize=csize

   djs_plot, ug, xqso, psym=-4, charsize=csize, /ynozero, $
    title='QSO Color-Refraction Plot', $
    xtitle='(u-g) [mag]', ytitle='u-band refraction [arcsec]'
   for j=0, nz-1, 10 do $
    xyouts, ug[j], xqso[j], 'z='+string(ztab[j],format='(f3.1)'), $
    charsize=csize

   !p.multi = [0,1,4]
   djs_plot, ztab, ug, /ynozero, psym=-4, $
    xtitle='Redshift', ytitle='(u-g)', charsize=csize
   djs_plot, ztab, gr, /ynozero, psym=-4, $
    xtitle='Redshift', ytitle='(g-r)', charsize=csize
   djs_plot, ztab, xqso, /ynozero, psym=-4, $
    xtitle='Redshift', ytitle='Refract(u) [asec]', charsize=csize
   djs_plot, ztab, yqso, /ynozero, psym=-4, $
    xtitle='Redshift', ytitle='Refract(g) [asec]', charsize=csize
   !p.multi = 0

   if (keyword_set(plotfile)) then dfpsclose

   urefract = xqso
   grefract = yqso
end

