; area - Survey area; default to 1000 deg^2

;------------------------------------------------------------------------------
; Return number per (h^-1 Mpc)^3 per (rest-frame r-band Petro mags).
; I make things brighter by 0.106 mag as the k+e-correction from
; Blanton's function at z=0.10 to z=0.
function lrg_phi, Rmag
   smallh = 1.0

   ; This lum function is from Blanton et al astro-ph/0210215
   phi_star = 1.49e-2 * smallh^3
   M_star = -20.44 + 5.*alog10(smallh) - 0.106
   alpha = -1.05
   phi = 0.4 * alog(10.) * phi_star * (10.^(0.4*(M_star-Rmag)))^(alpha+1) $
    * exp(-10^(0.4*(M_star-Rmag)))

   return, phi
end

;------------------------------------------------------------------------------
function lrgmodel_colors, filtloglam, filtcurve, zarr, $
 ageburst=ageburst1, zmetal=zmetal1, OmegaM=OmegaM, OmegaL=OmegaL

   common com_lrgmodel_photoz, allwave, allflux, agevec

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
            agevec = bcdat.age / 1e9 ; convert to Gyr
            allflux = fltarr(npix, nage, nmetal)
         endif
         allflux[*,*,imetal] = bcdat.flux
      endfor
      ; Convert to f_nu
      flambda2fnu = allwave^2 / 2.99792d18
      for imetal=0, nmetal-1 do $
       for iage=0, nage-1 do $
        allflux[*,iage,imetal] = allflux[*,iage,imetal] * flambda2fnu
   endif

   ;----------
   ; Compute the fluxes as a function of redshift.

   numz = n_elements(zarr)
   synflux = dblarr(nfilt,numz)
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
pro lrgsim, area1, fsystem1

   if (keyword_set(area1)) then area = area1 $
    else area = 1000.
   if (keyword_set(fsystem1)) then fsystem = fsystem1 $
    else fsystem = 'SDSS'
   nfilt = 5 ; ???

   OmegaM = 0.27
   OmegaL = 0.73

   mrange = [-24, -21.84]
mrange = [-24, -20.5]
   zrange = [0., 0.9]
   deltaz = 0.001
   deltam = 0.01
   cspeed = 3.e5
   iseed = 12345

   ;----------
   ; Construct bins in redshift
   ; Compute dVolume in units of (Mpc/h)^3

   znum = long((zrange[1] - zrange[0]) / deltaz) + 1
   zvec = (findgen(znum) + 0.5) * deltaz + zrange[0]
   ; Compute the volume in (Mpc/h)^3 per redshift slice
   dVolume = deltaz * dcomvoldz(zvec, OmegaM, OmegaL) * (cspeed/100)^3 $
    * (!pi/180.)^2 * area
   Dlum = lumdis(zvec, OmegaM, OmegaL) / lumdis(1e-3/cspeed, OmegaM, OmegaL)
   dmodulus = 5.*alog10(Dlum)

   ;----------
   ; Construct bins in absolute luminosity

   mnum = long((mrange[1] - mrange[0]) / deltam) + 1
   Mabsvec = Mrange[0] + (findgen(mnum) + 0.5) * deltam
   phi = lrg_phi(Mabsvec) * deltam

   ;----------
   ; Construct number density in z-M

   zarr = rebin(zvec,znum,mnum)
   Mabsarr = transpose(rebin(Mabsvec,mnum,znum))
   narr = dVolume # phi

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
   outdat = replicate(create_struct( $
    'z', 0., $
    'Mabs', 0., $
    'dmodulus', 0., $
    'mag', fltarr(nfilt)), ntot)

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
   ; Avoid discretization at the bin boundaries by adding small random numbers
   outdat.z = outdat.z + (randomu(iseed,ntot) - 0.5) * deltaz
   outdat.Mabs = outdat.Mabs + (randomu(iseed,ntot) - 0.5) * deltam

   ; Convert absolute magnitudes to apparent
   linterp, zvec, dmodulus, outdat.z, thismodulus
   outdat.dmodulus = thismodulus

   ;----------
   ; Compute the galaxy fluxes in the specified filters

   lrg_filters, fsystem, filtwave, filtcurve
   synflux = lrgmodel_colors(alog10(filtwave), filtcurve, zvec, $
    ageburst=ageburst, zmetal=zmetal, OmegaM=OmegaM, OmegaL=OmegaL)
   synmag = -2.5 * alog10(synflux)
   ; Normalize these to an absolute mag of zero in r-band at z=0
   synmag = synmag - synmag[2,0]
   for ifilt=0, nfilt-1 do begin
      linterp, zvec, reform(synmag[ifilt,*]), outdat.z, thismag
      outdat.mag[ifilt] = thismag + outdat.Mabs + outdat.dmodulus
   endfor

   ;----------
   ; Select the Padmanabhan photo-z sample, and plot it
   ; Before selection, add typical sky errors to the fluxes,
   ; plus an error that's 1% of the object flux
   ; I make things brighter by 0.07 mag to approximately convert from
   ; Blanton's Petrosian magnitudes to model magnitudes.
   ; Our Cut I should be Petro mags, and Cut II should be model mags.

   pmags = outdat.mag
   nmgy_err = [0.4, 0.2, 0.3, 0.5, 1.7]
   perr = 0 * pmags
   for ifilt=0, 4 do $
    perr[ifilt,*] = randomn(iseed,ntot) $
     * (nmgy_err[ifilt] + 0.01 * 10.^((22.5-reform(pmags[ifilt,*]))/2.5))
   for ifilt=0, 4 do $
    pmags[ifilt,*] = 22.5 $
     - 2.5*alog10( (10.^((22.5-pmags[ifilt,*])/2.5) + perr[ifilt,*]) > 0.01)

   rmag = reform(pmags[2,*])
   imag = reform(pmags[3,*])
   grcolor = reform(pmags[1,*] - pmags[2,*])
   ricolor = reform(pmags[2,*] - pmags[3,*])
   c_perp = ricolor - grcolor/4. - 0.18
   d_perp = ricolor - grcolor/8.
   c_par = 0.7*grcolor + 1.2*(ricolor - 0.18)
   qcut1 = (abs(c_perp) LT 0.2) $
    AND (rmag LT (13.6 + c_par/0.3)) $
    AND (rmag LT 19.7)
   qcut2 = (d_perp GT 0.55) $
    AND (grcolor GT 1.4) $
    AND (imag-0.07 LT (18.3 + 2.*d_perp)) $
    AND (imag-0.07 LT 20)

   splot,outdat.z,outdat.mag[3],ps=3, xtitle='z', ytitle='i-mag'
   i = where(qcut1 EQ 1 AND qcut2 EQ 0)
   soplot,outdat[i].z,outdat[i].mag[3],ps=3,color='green'
   j = where(qcut1 EQ 0 AND qcut2 EQ 1)
   soplot,outdat[j].z,outdat[j].mag[3],ps=3,color='red'
   k = where(qcut1 EQ 1 AND qcut2 EQ 1)
   soplot,outdat[k].z,outdat[k].mag[3],ps=3,color='blue'

   qpadman = (qcut1 OR qcut2) AND outdat.z GT 0.2 ; AND outdat.z LT 0.6
   npadman = total(qpadman)
   splog, 'Number in Padmanabhan sample = ', npadman
   splog, 'Areal density = ', npadman/area, ' / deg^2'
   splog, 'Volume density = ', $
    npadman/total(dvolume[0.20/deltaz:0.60/deltaz]), ' / (h^-1 Mpc)^3'
stop

   ; Plot the redshift distribution, and the same smoothed by 0.03 RMS in z
   zhist = float(histogram(outdat[where(qpadman)].z, binsize=deltaz, $
    min=zrange[0], max=zrange[1]))
   sigerr = 0.03 ; photo-z error
   nghalf = long(0.1 * n_elements(zvec))
   gkern = exp(-0.5*((findgen(nghalf*2+1)-nghalf)*deltaz/sigerr)^2)
   zhist2 = convol(zhist, gkern/total(gkern), /center, /edge_truncate)
   splot, zvec, zhist, xtitle='z', ytitle='Number', psym=10
   soplot, zvec, zhist2, color='red', psym=10

stop
end
;------------------------------------------------------------------------------
