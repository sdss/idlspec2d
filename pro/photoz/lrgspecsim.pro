; old - If set, then assume the SDSS-I telescope fiber sizes,
;       integration times, and efficiencies; otherwise assume 
;       2 arcsec fibers, 7200 sec exposures, and 40% efficiency
; imag - i-band magnitudes for objects at all redshifts
; resol - If set, then resample the spectra at this resolution, assuming
;         3 pixels per resolution element
pro lrgspecsim, old=old, zarr=zarr, imag=imag, resol=resol, zans=zans

   if (NOT keyword_set(zarr)) then $
    zarr = 0.01 + 0.01 * findgen(100)
   if (NOT keyword_set(imag)) then imag = 0*zarr + 21.2
   if (n_elements(zarr) NE n_elements(imag)) then $
    message, 'Number of elements in ZARR and IMAG must agree!'
old = 1 ; ???

   exptime = 2700.
   nexp = 3 ; number of exposures
   rnoise = 12. ; effective read noise in photons per wavelength bin
   tel_area = !pi * (125.^2 - 65.^2) ; cm^2
   iseed = 1234
   eigenfile = 'spEigenGal-*.fits'
   eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')

   loglam = 3.5780d0 + dindgen(4000)*1d-4
   npix = n_elements(loglam)
   wave = 10^loglam
   dwave = abs([wave[1:npix-1]-wave[0:npix-2], wave[npix-1]-wave[npix-2]]) ; Ang

   ; Return the photons/(10^-17 erg/s) for the given pixels + exposure time
   photons_per_flux = exptime $
    * ( spthroughput(402, camname='b1', expnum=6481, $
     efficiency=eff_b, loglam=loglam, /median) $
    + spthroughput(402, camname='r1', expnum=6481, $
     efficiency=eff_r, loglam=loglam, /median) )
;   tel_efficiency = eff_b + eff_r
;   photons_per_flux= tel_area * tel_efficiency * exptime * (wave * 1e-8) $
;    * 1e-17 * dwave^2 / (6.62e-27 * 3.e10)

   if (keyword_set(old)) then begin
      fibsize = 3. ; arcsec
   endif else begin
      fibsize = 2. ; arcsec
      ; Boost the efficiencies above 7000 Ang to be nearly flat
      photons_per_flux = 1.2 * photons_per_flux*(( (wave-6000)>0)/2000.>1)
   endelse

   ; Read in a representative sky spectrum
   readspec, 402, 16, mjd=51793, loglam=sky_loglam, sky=sky_flux
   sky_flux = sky_flux * (fibsize/3.)^2
   sky = interpol(sky_flux, sky_loglam, loglam)

   ; Simulate LRG spectra, where these evolve with redshift
   flambda2fnu = (wave*wave / 2.99792e18)
   numz = n_elements(zarr)
   lrgmodel_spectra, loglam, zarr=zarr, objflux=allflux
   allivar = 0. * allflux

   ; Broaden by velocity dispersion + instrumental response ???
for iz=0L, numz-1L do allflux[*,iz] = smooth(allflux[*,iz],3)

   for iz=0L, numz-1L do begin
      zthis = zarr[iz]

      ; Scale the galaxy spectrum to the requested i-band magnitude
      thisflux = allflux[*,iz]
      res = filter_thru(thisflux*flambda2fnu, waveimg=wave)
      synmag = -2.5 * alog10(res) - 48.6 + 2.5*17.0
      thisflux = thisflux * 10.^((synmag[3] - imag[iz])/2.5)

      ; Compute photon noise, and add it to the data
      qgood = (sky GT 0) AND (photons_per_flux GT 0)
      objerr = sqrt( (thisflux>0) / (photons_per_flux + (qgood EQ 0)))
      skyerr = sqrt( (sky>0) / (photons_per_flux + (qgood EQ 0)))
      ccderr = rnoise * sqrt(nexp) / (photons_per_flux + (qgood EQ 0))
      ; Following is the floor on the errors for bright objects... 
      adderr = 0.03 * (thisflux > 0) / sqrt(nexp)
      toterr = sqrt(objerr^2 + skyerr^2 + ccderr^2 + adderr^2)
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

   mkhdr, objhdr, 0
   sxaddpar, objhdr, 'COEFF0', loglam[0]
   sxaddpar, objhdr, 'COEFF1', loglam[1]-loglam[0]
   res_gal = zfind(allflux, allivar, hdr=objhdr, $
    starflux=starflux, starloglam0=starloglam0, stardloglam=stardloglam, $
    npoly=3, zmin=-0.01, zmax=1.00, pspace=2, nfind=5, width=10)

stop
vdiff = (zarr - res_gal[0,*].z) * 3e5
verr = res_gal[0,*].z_err * 3e5
vmax = 2000
splot, zarr, (vdiff < vmax) > (-vmax), psym=4

   zans = res_gal[0,*]
   return
end
pro lrgspecsim_test
   readspec, 400+lindgen(1), zans=zans, plug=plug
   indx = where((plug.primtarget AND 2L^5+2L^26) NE 0 $
    AND strmatch(zans.class,'GALAXY*') AND zans.zwarning EQ 0)
   imag = 22.5 - 2.5*alog10(zans.spectrosynflux[3]>1)
   lrgspecsim, zarr=zans[indx].z, imag=imag[indx], /old
end

