
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
; Apply a linear correction to fluxes derived from White Dwarf Spectra

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

