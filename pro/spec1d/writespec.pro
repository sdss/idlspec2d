pro writespec, plate, fiberid, mjd=mjd, filename=filename

   readspec, plate, fiberid, flux=objflux, flerr=objerr, loglam=loglam, $
    zans=zans
   wave = 10^loglam
   synflux = synthspec(zans, loglam=loglam)

   platestr = string(plate, format='(i4.4)')
   mjdstr = string(zans.mjd, format='(i5.5)')
   fibstr = string(zans.fiberid, format='(i3.3)')

   if (NOT keyword_set(filename)) then $
    filename = 'spec-' + platestr + '-' + mjdstr + '-' + fibstr + '.dat'

   get_lun, olun
   openw, olun, filename

   printf, olun, '# Plate ' + platestr
   printf, olun, '# MJD ' + mjdstr
   printf, olun, '# Fiber ' + fibstr
   printf, olun, '# Wavelength Flux   Error'
   printf, olun, '# [Ang]      [10^(-17) erg/cm/s/Ang]  [10^(-17) erg/cm/s/Ang]'
   for i=0, n_elements(objflux)-1 do begin
      printf, olun, wave[i], objflux[i], objerr[i], $
       format='(f10.3, e12.4, e12.4)'
   endfor

   close, olun
   free_lun, olun

end

