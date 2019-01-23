;+
; NAME:
;   rm_spec2mag
;
; PURPOSE:
;   Given a single input spectrum, return synthetic AB magnitude in SDSS
;   ugriz bands using filter_thru in idlutils.
;
; CALLING SEQUENCE:
;
;   synmag = rm_spec2mag(objwave, objflux, objivar, synmag_err=synmag_err, ntrial=)
;
; INPUTS:
;
;   objwave  - vaccum wavelength [npix]
;   objflux  - flux array [npix]; 1d-17 erg/s/cm^2/A
;   objivar  - inverse variance [npix]
;
; OUTPUTS:
;
;   synmag   - synthetic SDSS (AB) magnitudes in ugriz [5]
;   synmag_err - magnitude errors estimated using Monte Carlo trials
;
;  OPTIONAL OUTPUTS:
;
;   synflux  - synthetic fluxes in ugriz bands, in units of erg/s/cm2/Hz
;   synfnu_err - error in synflux, in same units

function rm_spec2mag, objwave1, objflux1, objivar1, synmag_err=synmag_err, $
  ntrial=ntrial, synflux=synflux, synfnu_err=synfnu_err

   ; number of MC trials in estimating magnitude errors
   if ~keyword_set(ntrial) then ntrial=50L
   npix = (size(objflux1, /dimens))[0]

   objwave = objwave1 & objflux = objflux1 & objivar = objivar1
   objerr1 = dblarr(npix)
   ind = where(objivar1 gt 1d-5)
   if ind[0] ne -1 then objerr1[ind] = 1./sqrt(objivar1[ind])

   flambda2fnu = objwave*objwave / 2.99792e18
   synflux = transpose(filter_thru(objflux*flambda2fnu, waveimg=objwave, $
       mask=(objivar LE 0),/toair))
   synthmag = fltarr(5)
   igood = where(synflux GT 0, ngood)
   if (ngood GT 0) then $
     synthmag[igood] = -2.5 * alog10(synflux[igood]) - 48.6 + 2.5*17.0
   synflux=synflux*1d-17

   ; now estiamte magnitude errors with MC trials
   objflux = dblarr(npix, ntrial)
   objivar = dblarr(npix, ntrial)
   for i=0L, ntrial-1 do begin
     objivar[*, i] = objivar1
     objflux[*, i] = objflux1 + randomn(seed, npix)*objerr1
   endfor
   flambda2fnu = (objwave*objwave / 2.99792e18) # replicate(1,ntrial)
   synflux1 = transpose(filter_thru(objflux*flambda2fnu, waveimg=objwave, $
    mask=(objivar LE 0),/toair))
   synthmag1 = fltarr(5,ntrial)
   igood = where(synflux1 GT 0, ngood)
   if (ngood GT 0) then begin
      synthmag1[igood] = -2.5 * alog10(synflux1[igood]) - 48.6 + 2.5*17.0
   endif
   synmag_err = (quantile_2d(0.84, synthmag1, dim=2) -  $
       quantile_2d(0.16, synthmag1, dim=2) ) * 0.5
   synfnu_err = (quantile_2d(0.84, synflux1, dim=2) -  $
       quantile_2d(0.16, synflux1, dim=2) ) * 0.5 * 1d-17

   ;message, 'stop' 

   return, synthmag

end
