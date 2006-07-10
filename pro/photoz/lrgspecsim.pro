; old - If set, then assume the SDSS-I telescope fiber sizes,
;       integration times, and efficiencies; otherwise assume 
;       2 arcsec fibers, 7200 sec exposures, and 40% efficiency
; imag - i-band magnitudes for objects at all redshifts
; resol - If set, then resample the spectra at this resolution, assuming
;         3 pixels per resolution element
pro lrgspecsim, old=old, imag=imag, resol=resol

   if (NOT keyword_set(imag)) then imag = 22.2

   tel_area = !pi * (125.^2 - 50.^2) ; cm^2
   fibsize = 2.0 ; fiber size in arcsec
   iseed = 1234
   eigenfile = 'spEigenGal-*.fits'
   eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')

   readspec, 402, 16, mjd=51793, zans=zans, loglam=loglam, wave=wave, $
    flux=flux, sky=sky, synflux=synflux, objhdr=objhdr
   npix = n_elements(loglam)
   zorig = zans.z
   sky = sky * (fibsize/3.)^2 ; Convert sky from 3 arcsec to 2 arcsec fiber
   dwave = [wave[1:npix-1]-wave[0:npix-2], wave[npix-1]-wave[npix-2]] ; Ang

   if (keyword_set(old)) then begin
      exptime = 2700. ; sec
      tel_efficiency = 0.15 / (((wave-6500)/1500)^4+1) ; best at 6500 Ang
      fibsize = 3. ; arcsec
      photons_per_flux = fcalib_default('b1',loglam,exptime) $
       + fcalib_default('r1',loglam,exptime)
   endif else begin
      exptime = 7200. ; sec
      tel_efficiency = 0.25
      fibsize = 2. ; arcsec
      photons_per_flux = tel_area * tel_efficiency * exptime $
       / (6.62e-27 * 3.e10 / (wave * 1e-8)) * 1e-17 / dwave
   endelse

   flambda2fnu = (wave*wave / 2.99792e18)

   zarr = 0.15 + 0.01 * findgen(85)
   numz = n_elements(zarr)
   allflux = fltarr(npix,numz)
   allivar = fltarr(npix,numz)

   for iz=0L, numz-1 do begin
      zthis = zarr[iz]

      ; Simulate the sky and galaxy spectra, and scale to the requested
      ; i-band magnitude
      pixshift = 1e4 * alog10((1.+zthis) / (1.+zorig))
      thisflux = sshift(synflux>0, pixshift)
      res = filter_thru(thisflux*flambda2fnu, waveimg=wave)
      synmag = -2.5 * alog10(res) - 48.6 + 2.5*17.0
      thisflux = thisflux * 10.^((synmag[3] - imag)/2.5)

      ; Compute photon noise, and add it to the data
      qgood = (sky GT 0) AND (photons_per_flux GT 0)
      objerr = sqrt( (thisflux>0) / (photons_per_flux + (qgood EQ 0)))
      skyerr = sqrt( (sky>0) / (photons_per_flux + (qgood EQ 0)))
      toterr = sqrt(objerr^2 + skyerr^2)
      allivar[*,iz] = qgood / (toterr^2 + (qgood EQ 0))
      allflux[*,iz] = thisflux + randomn(iseed, npix) * toterr
   endfor

   ;----------
   ; Read in the galaxy templates

   allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
   if (ct EQ 0) then $
    message, 'Unable to find EIGENFILE matching '+eigenfile
   thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]
   splog, 'Selecting EIGENFILE=' + thisfile
   starflux = readfits(thisfile, shdr,/silent)
   starloglam0 = sxpar(shdr, 'COEFF0')
   stardloglam0 = sxpar(shdr, 'COEFF1')

   ;----------
   ; Bin all spectra

   if (keyword_set(resol)) then begin
      ; c/RESOL should be the full-width of the resolution in km/sec,
      ; with an assumption of 3 pixels to sample that full width.
      rebinfac = 3.e5/(69.1*3) / resol
      gsigma = rebinfac * 3. / 2.355
      orig_sigma = 70.
      if (gsigma GT orig_sigma) then begin
         gkern = gauss_kernel(sqrt(gsigma^2 - orig_sigma^2))
         allflux = convol(allflux, gkern, /center, /edge_truncate, /normalize)
      endif
      rebinfac = round(rebinfac)
      npix = floor(npix/rebinfac)
      allflux = rebinfac * rebin(allflux[0:npix*rebinfac-1,*], npix, numz)
      allvar = (allivar GT 0) / (allivar + (allivar EQ 0))
      allvar = rebinfac * rebin(allvar[0:npix*rebinfac-1,*], npix, numz)
      allivar = (allvar GT 0) / (allvar + (allvar EQ 0))
      sxaddpar, objhdr, 'COEFF0', loglam[0] + 0.5*(rebinfac-1.)*1d-4
      sxaddpar, objhdr, 'COEFF1', rebinfac*1d-4

      dims = size(starflux,/dimens)
      nspix = dims[0]
      if (size(starflux,/n_dimen) EQ 1) then nstar = 1 $
       else nstar = dims[1]
      nsnew = floor(float(nspix)/rebinfac)
      starflux = rebinfac * rebin(starflux[0:nsnew*rebinfac-1,*], nsnew, nstar)
   endif

   ;----------
   ; Now try fitting a redshift

   res_gal = zfind(allflux, allivar, hdr=objhdr, $
    starflux=starflux, starloglam0=starloglam0, stardloglam=stardloglam, $
    npoly=3, zmin=-0.01, zmax=1.00, pspace=2, nfind=5, width=10)

stop
vdiff = (zarr - res_gal[0,*].z) * 3e5
zmax = 2000
splot, zarr, (vdiff < zmax) > (-zmax), psym=4

end

