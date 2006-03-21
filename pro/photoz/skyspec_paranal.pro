;+
; NAME:
;   skyspec_paranal
;
; PURPOSE:
;   Return the sky spectrum at Paranal
;
; CALLING SEQUENCE:
;   skyflux = skyspec_paranal(loglam, [ disp= ] )
;
; INPUTS:
;   loglam     - Log10 of wavelengths [vacuum Ang]; need not be uniformly
;                spaced
;
; OPTIONAL INPUTS:
;   disp       - Sigma of the intrumental response (dispersion) in units
;                of the LOGLAM pixels, either as a scalar, or as a vector
;                with the same number of elements as LOGLAM.
;                If not set, then the returned sky spectrum is at the
;                resolution of the Paranal data files, and will be
;                undersampled for any resolution worse than about 100,000.
;
; OUTPUTS:
;   skyflux    - Sky spectrum in units of e-17 erg/s/cm^2/Ang/asec^2
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The lowest sky level is on plate 1292/52736 taken at airmass=1.04.
;   A representative low sky level is plate 406/51817.
;   A high sky level is plate 298/51955.
;
; EXAMPLES:
;
; BUGS:
;   The Paranal sky spectra are missing data from 8550-8610 Ang.
;   There are a few regions where their spectra are not in agreement
;   with the lower-resolution SDSS sky spectra such as at 7600-7650 Ang,
;   probably because the SDSS spectra "correct up" the sky level for
;   the telluric absorption.  Also, the Paranal sky is higher than SDSS
;   by about a factor of 2 blueward of 4800 Ang.
;
; DATA FILES:
;   $IDLSPEC2D_DIR/templates/fluxed_sky_346.fits
;   $IDLSPEC2D_DIR/templates/fluxed_sky_580L.fits
;   $IDLSPEC2D_DIR/templates/fluxed_sky_860L.fits
;   $IDLSPEC2D_DIR/templates/fluxed_sky_437.fits
;   $IDLSPEC2D_DIR/templates/fluxed_sky_580U.fits
;   $IDLSPEC2D_DIR/templates/fluxed_sky_860U.fits
;   $IDLSPEC2D_DIR/templates/fluxed_sky_564U.fits
;   $IDLSPEC2D_DIR/templates/fluxed_sky_800U.fits
;
; PROCEDURES CALLED:
;   djs_maskinterp()
;   mrdfits()
;   populate_image
;   sxpar()
;
; REVISION HISTORY:
;   20-Mar-2006  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
function skyspec_paranal, loglam, disp=disp

   if (n_elements(loglam) LT 1) then $
    message, 'LOGLAM must have at least two elements!'

   npix = n_elements(loglam)
;   loglam0 = loglam[0]
;   dloglam = loglam[1] - loglam[0]

   ; If we need to convolve the spectrum with a Gaussian dispersion,
   ; then pad the spectra that we generate, and trim at the end
   if (keyword_set(disp)) then nhalf = ceil(4 * (max(disp)>1)) $
    else nhalf = 0

   skyflux = fltarr(npix+2*nhalf)
   skymask = fltarr(npix+2*nhalf)

   files = 'fluxed_sky_' + ['346','437','564U','580L','580U','860L','860U'] $
    + '.fits'
   nfile = n_elements(files)
   for ifile=0, nfile-1 do begin
      thisfile = files[ifile]
      ; The factor of 10 below converts to units of e-17 erg/s/cm^2/Ang/asec^2
      thisflux = 10. * mrdfits(filepath(thisfile, $
       root_dir=getenv('IDLSPEC2D_DIR'), subdir='templates'), 0, hdr)
      thiswave = sxpar(hdr,'CRVAL1') $
       + (dindgen(sxpar(hdr,'NAXIS1'))+1-sxpar(hdr,'CRPIX1')) $
       * sxpar(hdr,'CDELT1')
      airtovac, thiswave
      ; Trim out some bad data points in this one file...
      if (thisfile EQ 'fluxed_sky_564U.fits') then begin
         thiswave = thiswave[300:n_elements(thiswave)-1001]
         thisflux = thisflux[300:n_elements(thisflux)-1001]
      endif
      igood = where(thisflux NE 0)
      thiswave = thiswave[igood]
      thisflux = thisflux[igood]
      thisloglam = alog10(thiswave)
;      thispix = (thisloglam - loglam0) / dloglam + nhalf
      thispix = interpol(findgen(npix), loglam, thisloglam) + nhalf
      populate_image, skyflux, thispix, weights=thisflux>0, assign='cic'
      populate_image, skymask, thispix, assign='cic'
   endfor

   skyflux = skyflux * (skymask GT 0) / (skymask + (skymask LE 0))

   ; Interpolate over missing wavelengths
   skyflux = djs_maskinterp(skyflux, skyflux LE 0, /const)

   if (keyword_set(disp)) then begin
      gpix = findgen(2*nhalf+1) - nhalf
      if (n_elements(disp) EQ 1) then begin
         skyflux = convol(skyflux, gauss1(gpix, [0.,disp[0],1.]))
      endif else begin
         if (n_elements(disp) NE npix) then $
          message, 'Number of elements in DISP and LOGLAM disagree!'
         origflux = skyflux
         for i=0, npix-1 do $
          if (disp[i] GT 0) then $
           skyflux[i+nhalf] = total(origflux[i:i+2*nhalf] $
            * gauss1(gpix, [0.,disp[i],1.]))
      endelse
      skyflux = skyflux[nhalf:npix-1+nhalf]
   endif

   return, skyflux
end
;------------------------------------------------------------------------------
