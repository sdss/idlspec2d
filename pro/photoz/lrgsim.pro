; area - Survey area; default to 1000 deg^2

;------------------------------------------------------------------------------
function lrgmodel_colors, filtloglam, filtcurve, zarr, $
 ageburst=ageburst1, zmetal=zmetal1

   common com_lrgmodel_photoz, allwave, allflux

   if (keyword_set(ageburst1)) then ageburst = ageburst1 $
    else ageburst = 2.5
   if (keyword_set(zmetal1)) then zmetal = zmetal1 $
    else zmetal = 0.018
   ndimf = size(filtcurve, /n_dimen)
   if (ndimf EQ 1) then nfilt = 1 $
    else nfilt = (size(filtcurve,/dimens))[1]
   if (n_elements(filtloglam) NE (size(filtcurve,/dimens))[0]) then $
    message, 'Number of elements in LOGLAM and FILTCURVE are inconsistent'

   ;----------
   ; Read the Bruzual-Charlot model spectra FITS files.

   metalstr = ['z008', 'z02', 'z05']
   metalvec = [0.008, 0.02, 0.05]
   agestr = [ '5Myr', '25Myr', '100Myr', '290Myr', '640Myr', '900Myr', $
    '1.4Gyr', '2.5Gyr', '5Gyr', '11Gyr' ]
   agevec = [0.005, 0.025, 0.100, 0.290, 0.640, 0.900, 1.4, 2.5, 5.0, 11.0]

   ; Read in an model LRG spectra, assuming the same wavelengths for all
   if (NOT keyword_set(allwave)) then begin
      nage = n_elements(agevec)
      nmetal = n_elements(metalvec)
      eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
      for iage=0, nage-1 do begin
         for imetal=0, nmetal-1 do begin
            eigenfile = 'ssp_' + agestr[iage] + '_' + metalstr[imetal] + '.spec'
            readcol, djs_filepath(eigenfile, root_dir=eigendir), $
             lambda1, flux1, comment='#', format='(F,F)'
            if (iage EQ 0 AND imetal EQ 0) then begin
               npix = n_elements(lambda1)
               allwave = lambda1
               allflux = fltarr(npix, nage, nmetal)
            endif
            ; Convert to f_nu
            flambda2fnu = lambda1^2 / 2.99792d18
            allflux[*,iage,imetal] = flux1 * flambda2fnu
         endfor
      endfor
   endif

   ;----------
   ; Compute the colors as a function of redshift.

   numz = n_elements(zarr)
   synflux = dblarr(nfilt,numz)
   for iz=0L, numz-1 do begin
      print, format='(" Z ",i5," of ",i5,a1,$)', $
        iz, numz, string(13b)

      ; Convert redshift to an age, using the WMAP cosmology,
      ; and say that these galaxies formed when the Universe
      ; was AGEBURST Gyr old
      hubble0 = 71. * 3.1558e7 / 3.0856e19 * 1e9 ; units of Gyr^-1
      thisage = lookback(1000., 0.27, 0.73) / hubble0 $
       - lookback((zarr[iz] > 0), 0.27, 0.73) / hubble0 - ageburst

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
      vactoair, thiswave
      thisloglam = alog10(thiswave)

      bigloglam = [thisloglam, filtloglam]
      bigloglam = bigloglam[uniq(bigloglam,sort(bigloglam))]
      dloglam = bigloglam - shift(bigloglam,1)
      dloglam[0] = dloglam[1]

      linterp, thisloglam, thisflux, bigloglam, bigflux1
      for ifilt=0, nfilt-1 do begin
         linterp, filtloglam, filtcurve[*,ifilt], bigloglam, bigfiltcurve
         sumfilt = total(bigfiltcurve * dloglam)
         synflux[ifilt,iz] = total(bigflux1 * bigfiltcurve * dloglam) $
          / (sumfilt + (sumfilt LE 0))
      endfor

      ; Space equally in log-wavelength
;      bigloglam = 3.2000d0 + dindgen(8800) * 1.d-4
;      bigwave = 10.d0^bigloglam
;      linterp, thiswave, thisflux, bigwave, bigspecflux
;      synflux[*,iz] = filter_thru( bigspecflux, $
;       waveimg=bigwave, /toair)
   endfor
   print

   return, synflux
end
;------------------------------------------------------------------------------
pro lrg_filters, fsystem, wave, filtcurve

; FOLD IN TELESCOPE + ATMOSPHERE THROUGHPUT !!!???

   case strupcase(fsystem) of
   'SDSS' : begin
      wave = 2980.d0 + findgen(8251)
      filtfiles = 'sdss_jun2001_' + ['u','g','r','i','z'] + '_atm.dat'
      nfile = n_elements(filtfiles)
      filtcurve = fltarr(n_elements(wave),nfile)
      for ifile=0, nfile-1 do begin
         filename = filepath(filtfiles[ifile], $
          root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
         readcol, filename, fwave, fthru
         linterp, fwave, fthru, wave, filtcurve1
         filtcurve[*,ifile] = filtcurve1
      end
      end
   'PS1' : begin
         filename = filepath('PS1.dat', $
          root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
         readcol, w_g, f_g, w_r, f_r, w_i, f_i, w_z, f_z, w_y, f_y
         wave = 10. * w_g ; convert from nm to Ang
         filtcurve = [[f_g],[f_r],[f_i],[f_z],[f_y]]
      end
   endcase

   return
end

;------------------------------------------------------------------------------
; Return number per (h^-1 Mpc)^3 per (rest-frame R-petro-mag)
function lrg_phi, Rmag
   smallh = 1.0
   mu = -22.29 + 5.*alog10(smallh)
   sig = 0.219
   nbar = 0.159e-4 * smallh^3
   q = -0.5
   qq = 1./q^2
   const = -0.4 * alog(10.)
   numer = (Rmag - mu) / (sig / const)
   phi = nbar * q * const * qq^qq / (sig * gamma(qq)) $
    * exp(qq * (q*numer - exp(q*numer)))
   return, phi
end

;------------------------------------------------------------------------------
pro lrgsim, area1, fsystem1

   if (keyword_set(area1)) then area = area1 $
    else area = 100.
   if (keyword_set(fsystem1)) then fsystem = fsystem1 $
    else fsystem = 'SDSS'

   OmegaM = 0.3
   OmegaL = 0.7

   mrange = [-24, -20]
   zrange = [0., 2.]
   deltaz = 0.01
   zmax = 2.0
   deltam = 0.02
   cspeed = 3.e5
   iseed = 12345

   ;----------
   ; Construct bins in redshift
   ; Compute dVolume in units of (Mpc/h)^3

   znum = long((zrange[1] - zrange[0]) / deltaz)
   zvec = (findgen(znum) + 0.5) * deltaz
   ; Compute the volume in (Mpc/h)^3 per deg^2 per redshift slice
   dVolume = deltaz * dcomvoldz(zvec, OmegaM, OmegaL) * (cspeed/100)^3 $
    * (!pi/180.)^2 * area
   Dlum = lumdis(zvec, OmegaM, OmegaL)
   dmodulus = 5.*alog10(dlum)

   ;----------
   ; Construct bins in absolute luminosity

   mnum = long((mrange[1] - mrange[0]) / deltam)
   Mabsvec = Mrange[0] + (findgen(mnum) + 0.5) * deltam
   phi = lrg_phi(Mabsvec)

   ;----------
   ; Construct number density in z-M

   zarr = rebin(zvec,znum,mnum)
   Mabsarr = transpose(rebin(Mabsvec,mnum,znum))
   narr = dvolume # phi

;foo = 0*narr ; ???
   ;----------
   ; Construct the random catalog of objects

   nsum = total(narr)
   ntot = long(nsum) > 1
   splog, 'Generating random numbers for ', ntot, ' objects'
   rand = randomu(iseed, ntot) * nsum
   rand = [rand[sort(rand)], nsum+1]
   j = 0L ; index into RAND
   tmpsum = 0.d0
   outdat = replicate(create_struct('z', 0., 'Mabs', 0.), ntot)
   for iz=0L, znum-1L do begin
      print,zvec[iz],string(13b),format='(" z=",f,a1,$)'
      for im=0L, mnum-1L do begin
         tmpsum = tmpsum + narr[iz,im]
         while (rand[j] LT tmpsum) do begin
;foo[iz,im] = foo[iz,im] + 1 ; ???
            outdat[j].z = zarr[iz,im]
            outdat[j].Mabs = Mabsarr[iz,im]
            j = j + 1
         endwhile
      endfor
   endfor
   print
   ; Avoid discretization at the bin boundaries by adding small randomu numbers
   outdat.z = outdat.z + (randomu(iseed,ntot) - 0.5) * deltaz
   outdat.Mabs = outdat.z + (randomu(iseed,ntot) - 0.5) * deltam

   ;----------
   ; Compute the galaxy fluxes in the specified filters

   lrg_filters, fsystem, filtwave, filtcurve
   synflux = lrgmodel_colors(alog10(filtwave), filtcurve, zvec, $
    ageburst=ageburst, zmetal=zmetal)

stop
end
;------------------------------------------------------------------------------
