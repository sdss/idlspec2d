;+
; NAME:
;   kurucz_restore
;
; PURPOSE:
;   Read back Kurucz models from FITS and apply an empirical correction
;
; CALLING SEQUENCE:
;   kurucz_restore, kwave, kflux, nkflux=, hdr=, kindx=, smoothpix=
;
; INPUTS:
; 
; OPTIONAL INPUTS:
;   smoothpix - smooth models by this number of pixels -- useful for 
;               comparisons to non-SDSS data
;   
; OUTPUT:
;   kwave - wavelength in Angstroms of Kurucz models (binning is in log-10 A)
;           [npix]
;   kflux - flux of Kuruz models [npix, nmodel]
;
; OPTIONAL OUTPUTS:
;   nkflux - normalized model flux [npix, nmodel]
;   hdr    - header of FITS file 
;   kindx  - structure containing info on each model (Teff, logg, [Fe/H], mags)
;
; COMMENTS:
;   The model file called "kurucz_stds_v5.fit' is assumed to be in 
;   $IDLSPEC2D_DIR/etc.  The models have been produced by the SPECTRUM code 
;   (R.0. Gray & C. J. Corbally, 1994, AJ, 107, 742) using the latest Kurucz
;   model atmospheres.  The 0th HDU contains the flux in ergs/s/cm^2/A -- 
;   the absolute value of the flux is arbitrary.  The spectra have been 
;   convolved to SDSS resolution (approximately) and rebinned to dloglam 
;   = 1e-4.  The 1st HDU contains information about each model such as 
;   effective temperature, surface gravity, and metallicity.  The wavelength 
;   information is in the header.  The second HDU contains information about
;   an empirical correction which should be applied to the models.
;
; BUGS:
;   The empirical correction which is applied is derived from comparing 
;   hot DA white dwarfs in SDSS to models.  This only corrects the low order 
;   shape of the spectra.  A different correction is probably needed for each 
;   T_eff and [Fe/H].
;
; EXAMPLES:
;   kurucz_restore, kwave, kflux
;   plot, kwave, kflux[*, 5]
;
; PROCEDURES CALLED:
;   filter_thru()
;   mrdfits()
;   rectify
;   sxpar()
;   traceset2xy
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   12-Aug-2003  Created by C. Tremonti, Steward Observatory
;
;-
;------------------------------------------------------------------------------

pro kurucz_restore, kwave, kflux, nkflux = nkflux, hdr = hdr, kindx = kindx, $
    smoothpix = smoothpix

kurucz_file = filepath('kurucz_stds_v5.fit', $
              root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')

kflux = mrdfits(kurucz_file, 0, hdr, /silent)  ; flux
kindx = mrdfits(kurucz_file, 1, /silent)
kset = mrdfits(kurucz_file, 2, /silent)

;----------------
; Create Wavelength Array

npix = n_elements(kflux[*,0])
nmod = n_elements(kflux[0,*])
crval = sxpar(hdr, 'CRVAL1')
kwave = 10.0^(lindgen(npix) * 1.0d-4 + crval)

;-----------------
; Apply a correction derived from White Dwarf Spectra

traceset2xy, kset, kwave, kfix
for i = 0, nmod - 1 do begin
  kflux[*,i] = kflux[*,i] / kfix
  fluxfnu = kflux[*,i] * kwave^2 / 2.99792e18
  fthru=filter_thru(fluxfnu, waveimg=kwave, $
                    filter_prefix = 'sdss_jun2001', /toair)
  kindx[i].mag = -2.5*alog10(fthru) - 48.6
endfor

;--------------------
; Smooth to lower resolution if desired

if keyword_set(smoothpix) then begin
  nkpix = long(5*smoothpix)
  kern = exp( -0.5 * (findgen(nkpix*2+1) - nkpix)^2 / smoothpix^2)
  kern = kern / total(kern)
  for imod = 0, nmod - 1 do kflux[*,imod] = convol(kflux[*,imod], kern, /center)
endif

;---------------
; Return rectified fluxes if desired

;nkflux = rectify(kflux, /mask, wave = kwave)
nkflux = rectify(kflux)

end

