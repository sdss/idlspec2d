
; Compute the refraction terms for QSOs.

pro qsorefract, eigenfile

;   if (NOT keyword_set(eigenfile)) then eigenfile = 'spEigenQSO-*.fits'
   if (NOT keyword_set(eigenfile)) then eigenfile = 'spEigenQSO.fits'

   ;----------
   ; Find the template spectrum file

   eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
   allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
   if (ct EQ 0) then $
    message, 'Unable to find EIGENFILE matching '+eigenfile
   thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]
   splog, 'Selecting EIGENFILE=' + thisfile

   ;----------
   ; Read the template spectrum

   flux = readfits(thisfile, shdr)
   flux = smooth(flux[*,0],25) ; Just take the mean spectrum (smoothed)
   naxis1 = sxpar(shdr, 'NAXIS1')
   loglam0 = sxpar(shdr, 'COEFF0')
   dloglam0 = sxpar(shdr, 'COEFF1')
   wave = 10.d0^(loglam0 + lindgen(naxis1) * dloglam0)
;flux = 0*flux + 1
;flux[where(wave GT 1200 AND wave LT 1300)] = 10

   ;----------
   ; Tabulate the index of refraction of air
   ; from Meggers & Peters 1919, ApJ 50, 56.

   rwave = 2000. + findgen(33) * 250.
   rval = [3255.82,3107.87,3014.05,2951.08,2906.85,2874.63,2850.43, $
    2831.79,2817.12,2805.36,2795.78,2787.88,2781.27,2775.69,2770.94, $
    2766.85,2763.31,2760.22,2757.51,2755.11,2752.99,2751.10,2749.40, $
    2747.88,2746.50,2745.25,2744.12,2743.09,2742.14,2741.28,2740.48, $
    2739.75,2739.07] * 1.0d-7

   ztab = 0.0 + 0.05 * findgen(100)
   nz = n_elements(ztab)
   rtab = fltarr(nz,5)

   for iz=0, nz-1 do begin
      ; Redshift the QSO spectrum
      tmpwave = wave * (1.0 + ztab[iz])
      tmpflux = flux

      ; Interpolate the index of refraction curve
      tmprval = interpol(rval, rwave, tmpwave)

      ; Integrate over the SDSS filters
      num1 = filter_thru(tmpflux * tmprval, waveimg=tmpwave, /toair)
      num2 = filter_thru(tmpflux, waveimg=tmpwave, /toair)
      rtab[iz,*] = num1 / num2

   endfor

   ; Normalize to arc-seconds at a zenith distance of 45 deg
   rtab = rtab * !radeg * 3600.

   ; Now make the dtheta-dtheta plot
   csize = 2.0
   rnorm = 0.5 * (rtab[*,2] + rtab[*,3]) ; Compare against ave. r+i position
   scale = 1. ; Re-scale plots to arcsec
   xplot = (rtab[*,0]-rnorm) * scale
   yplot = (rtab[*,1]-rnorm) * scale
   dfpsplot, 'qso-refract.ps'
   djs_plot, xplot, yplot, psym=-4, charsize=csize, /ynozero, /square, $
    title='QSO Refraction-Refraction Plot', $
    xtitle='u-band refraction [arcsec]', ytitle='g-band refraction [arcsec]'
   for j=0, nz-1, 10 do $
    xyouts, xplot[j], yplot[j], 'z='+string(ztab[j],format='(f3.1)'), $
    charsize=csize
   dfpsclose

end

