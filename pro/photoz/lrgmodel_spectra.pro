
;------------------------------------------------------------------------------
pro lrgmodel_spectra, loglam, zarr=zarr1, $
 ageburst=ageburst1, zmetal=zmetal1, OmegaM=OmegaM1, OmegaL=OmegaL1, $
 objflux=objflux, synflux=synflux

   common com_lrgmodel_spectra, allwave, allflux, agevec

   if (keyword_set(zarr1)) then zarr = zarr1 $
    else zarr = 0
   if (keyword_set(ageburst1)) then ageburst = ageburst1 $
    else ageburst = 2.5
   if (keyword_set(zmetal1)) then zmetal = zmetal1 $
    else zmetal = 0.018
   if (n_elements(OmegaM1) GT 0) then OmegaM = OmegaM1[0] $
    else OmegaM = 0.27
   if (n_elements(OmegaL1) GT 0) then OmegaL = OmegaL1[0] $
    else OmegaL = 0.73
   nfpix = n_elements(loglam)

   ;----------
   ; Read the Bruzual-Charlot model spectra FITS files.

   metalstr = ['z008', 'z02', 'z05']
   metalvec = [0.008, 0.02, 0.05]

   ; Read in a model LRG spectra, assuming the same wavelengths for all
   if (NOT keyword_set(allwave)) then begin
      nmetal = n_elements(metalvec)
      eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
      for imetal=0, nmetal-1 do begin
         eigenfile = 'bc03_padova1994_chab_' + metalstr[imetal] + '_ssp.fit.gz'
         bcdat = mrdfits(djs_filepath(eigenfile, root_dir=eigendir), 1)
         if (imetal EQ 0) then begin
            npix = n_elements(bcdat.wave)
            nage = n_elements(bcdat.age)
            allwave = bcdat.wave
            airtovac, allwave ; Convert to vacuum!!
            agevec = bcdat.age / 1e9 ; convert to Gyr
            allflux = fltarr(npix, nage, nmetal)
         endif
         allflux[*,*,imetal] = bcdat.flux
      endfor
      ; Convert to f_nu
; ???
;      flambda2fnu = allwave^2 / 2.99792d18
      for imetal=0, nmetal-1 do $
       for iage=0, nage-1 do $
        allflux[*,iage,imetal] = allflux[*,iage,imetal]
;        allflux[*,iage,imetal] = allflux[*,iage,imetal] * flambda2fnu
   endif

   ;----------
   ; Compute the fluxes as a function of redshift.

   numz = n_elements(zarr)
   if (arg_present(synflux)) then synflux = dblarr(nfilt,numz)
   if (arg_present(objflux)) then objflux = dblarr(nfpix,numz)

   for iz=0L, numz-1 do begin
      print, format='(" Z ",i5," of ",i5,a1,$)', $
        iz, numz, string(13b)

      ; Convert redshift to an age, using the WMAP cosmology,
      ; and say that these galaxies formed when the Universe
      ; was AGEBURST Gyr old
      hubble0 = 71. * 3.1558e7 / 3.0856e19 * 1e9 ; units of Gyr^-1
      thisage = lookback(1000., OmegaM, OmegaL) / hubble0 $
       - lookback((zarr[iz] > 0), OmegaM, OmegaL) / hubble0 - ageburst

      ; Demand that we stay in the bounds (at least the lower bounds)
      ; of the models.  Specifically, at the minimum, we don't want
      ; to be taking the logarithm of zero or negative numbers for
      ; these parameters when we interpolate in log-log space.
      thisage = thisage > agevec[0]
      zmetal = zmetal > metalvec[0]

      ; Linearly interpolate from the templates in log(age),log(zmetal),
      ; vs. log(flux) space.
      i0 = ((reverse(where(agevec LT thisage)))[0] > 0) $
       < (n_elements(agevec)-2)
      j0 = ((reverse(where(metalvec LT zmetal)))[0] > 0) $
       < (n_elements(metalvec)-2)
      i1 = i0 + 1
      j1 = j0 + 1
      agewts = [alog10(agevec[i1]/thisage), -alog10(agevec[i0]/thisage)]
      metwts = [alog10(metalvec[j1]/zmetal), -alog10(metalvec[j0]/zmetal)]
      agewts = agewts / total(agewts)
      metwts = metwts / total(metwts)

      thisflux = 10.d0^( $
       agewts[0] * metwts[0] * alog10(allflux[*,i0,j0]) $
       + agewts[0] * metwts[1] * alog10(allflux[*,i0,j1]) $
       + agewts[1] * metwts[0] * alog10(allflux[*,i1,j0]) $
       + agewts[1] * metwts[1] * alog10(allflux[*,i1,j1]) )

      thiswave = allwave * (1 + zarr[iz])
      thisloglam = alog10(thiswave)

      thispix = interpol(dindgen(nfpix), loglam, thisloglam)
      newmask = fltarr(nfpix)
      newflux = fltarr(nfpix)
      populate_image, newmask, thispix, assign='cic'
      populate_image, newflux, thispix, weights=thisflux, assign='cic'
      qgood = newmask GT 0
      newflux = qgood * newflux / (newmask + (qgood EQ 0))
      newflux = djs_maskinterp(newflux, newmask LE 0, /const)

      if (arg_present(objflux)) then objflux[*,iz] = newflux

; ???
      if (arg_present(synflux)) then begin
         dloglam = 0.5 * abs(shift(filtloglam,-1) - shift(filtloglam,1))
         dloglam[0] = dloglam[1]
         dloglam[nfpix-1] = dloglam[nfpix-2]
         for ifilt=0, nfilt-1 do begin
            sumfilt = total(filtcurve[*,ifilt] * dloglam)
            synflux[ifilt,iz] = total(bigflux1 * filtcurve[*,ifilt] * dloglam) $
             / (sumfilt + (sumfilt LE 0))
         endfor
      endif
   endfor
   print

   return
end
;------------------------------------------------------------------------------
