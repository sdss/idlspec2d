;+
; NAME:
;   flux2maggie
;
; PURPOSE:
;   Integrate a spectrum under filter curves to get maggies or magnitudes
;
; CALLING SEQUENCE:
;   spectroflux = flux2maggie( objflux, waveimg=, _EXTRA=EXTRA, mag=mag )
;
; INPUTS:
;   objflux    - Flux image in units of 10^-17 erg/s/cm^2/Ang [NPIX,NOBJ]
;   waveimg    - Wavelength image in Angstroms [NPIX,NOBJ], or this could
;                be a single vector if the wavelength mapping is the same
;                for all traces (note the latter is faster to compute)
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;   EXTRA      - Keywords for FILTER_THRU()
;
; OUTPUTS:
;   spectroflux- Integrated response in all 5 SDSS filters, ordered ugriz;
;                dimensions are [NOBJ,5] or [5] if NOBJ=1
;
; OPTIONAL OUTPUTS:
;   mag        - Magnitudes in all 5 SDSS filters, ordered ugriz;
;                dimensions are [NOBJ,5] or [5] if NOBJ=1
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   filter_thru()
;
; REVISION HISTORY:
;   08-Apr-2009  Written by D. Schlegel, LBNL
;-
;------------------------------------------------------------------------------
function flux2maggie, objflux, waveimg=waveimg, mag=mag, $
 _EXTRA=KeywordsForFilterThru

   dims = size(objflux, /dimens)
   npix = dims[0]
   if (npix EQ 0) then $
    message, 'FLUX not defined'
   if (n_elements(dims) EQ 1) then nobj = 1 $
    else nobj = dims[1]

   dim2 = size(waveimg, /dimens)
   if (dim2[0] NE dims[0]) then $
    message, 'Dimensions of FLUX and WAVE do not agree'
   if (n_elements(dim2) GT 1) then begin
      if (dim2[1] NE dims[1]) then $
       message, 'Dimensions of FLUX and WAVE do not agree'
   endif

   spectroflux = fltarr(5,nobj)

   flambda2fnu = waveimg^2 / 2.99792e18
   fthru = filter_thru(objflux * rebin(flambda2fnu,npix,nobj), $
    waveimg=waveimg, _EXTRA=KeywordsForFilterThru)
   spectroflux = transpose(fthru) * 10^((22.5 + 48.6 - 2.5*17.)/2.5)

   if (arg_present(mag)) then $
    mag = 22.5 - 2.5*alog10(spectroflux)

   return, spectroflux
end
;------------------------------------------------------------------------------
