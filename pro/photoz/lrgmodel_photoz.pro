;+
; NAME:
;   lrgmodel_photoz
;
; PURPOSE:
;   Simple photo-z finder for LRGs using Bruzual-Charlot models.
;
; CALLING SEQUENCE:
;   zfit = lrgmodel_photoz(pflux, pflux_ivar, [ /abcorrect, extinction=, $
;    abfudge=, ageburst=, zmetal=, filterlist=, adderr=, z_err=, chi2= ] )
;
; INPUTS:
;   pflux          - Object fluxes in the 5 SDSS filters [5,NOBJ]
;   pflux_ivar     - Inverse variances for FLUX [5,NOBJ]
;
; OPTIONAL INPUTS:
;   abcorrect      - If set, then convert the input fluxes from the 2.5-m
;                    natural system to AB fluxes
;   extinction     - If set, then apply these extinction corrections [5,NOBJ]
;   abfudge        - Additional AB "fudge factors"; default to adding
;                    [0,0,0,0,0] mag to input magnitudes, where a positive
;                    value makes that flux fainter
;   ageburst       - Age of the Universe at the single-burst; default
;                    to value from previous call, or 2.5 Gyr
;   zmetal         - Metallicity at the single-burst; default
;                    to value from previous call, or [Fe/H] = 0.025 dex.
;   filterlist     - List of filter indices to use in fits; default to
;                    using all five filters [0,1,2,3,4]
;   adderr         - Fractional error to add in quadrature; default to 0.03
;
; OUTPUTS:
;   zfit           - Best-fit redshift [NOBJ]
;
; OPTIONAL OUTPUTS:
;   z_err          - Redshift error, or a negative value if an error
;                    occurred in the quadratic minimization estimate [NOBJ]
;   chi2           - Best-fit chi^2 [NOBJ]
;
; COMMENTS:
;   The fluxes should be AB fluxes, or SDSS 2.5-m natural system fluxes
;   if /ABCORRECT is set.
;   The fluxes should already be extinction-corrected, unless
;   the EXTINCTION keyword is passed.
;
; EXAMPLES:
;
; BUGS:
;   The LRG template is not quite correct.  I have modified the spectrum
;   on large scales by multiplying by an empirically-determined function.
;   Also, I've extrapolated the two ends of the spectrum
;   to be essentially flat in f_lambda.
;   The 3-sigma-clipped scatter appears to be 0.023 at z>0.1.
;
; PROCEDURES CALLED:
;   computechi2()
;   filter_thru()
;   find_nminima()
;   readcol
;   sdssflux2ab
;   sxpar()
;
; REVISION HISTORY:
;   09-Dec-2003  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function lrgmodel_photoz, pflux, pflux_ivar, z_err=z_err, $
 abcorrect=abcorrect, extinction=extinction, abfudge=abfudge, $
 ageburst=ageburst1, zmetal=zmetal1, $
 filterlist=filterlist, adderr=adderr, chi2=chi2

   common com_lrgmodel_photoz, zarr, synflux, ageburst_save, zmetal_save

   if (n_elements(filterlist) EQ 0) then filterlist = lindgen(5)
   if (n_elements(adderr) EQ 0) then adderr = 0.03
   if (n_elements(abfudge) EQ 0) then abfudge = [0,0,0,0,0]
   if (n_elements(ageburst_save) EQ 0) then ageburst_save = 2.5
   if (n_elements(zmetal_save) EQ 0) then zmetal_save = 0.025

   if (keyword_set(ageburst1)) then ageburst = ageburst1 $
    else ageburst = 2.5
   if (keyword_set(zmetal1)) then zmetal = zmetal1 $
    else zmetal = 0.016

   ;----------
   ; Read the Bruzual-Charlot model spectra FITS files.
   ; This need only be done the first time this function is called,
   ; then cached for future calls.

   if (NOT keyword_set(zarr)) then begin

      metalstr = ['z008', 'z02', 'z05']
      metalvec = [0.008, 0.02, 0.05]
      agestr = [ '5Myr', '25Myr', '100Myr', '290Myr', '640Myr', '900Myr', $
       '1.4Gyr', '2.5Gyr', '5Gyr', '11Gyr' ]
      agevec = [0.005, 0.025, 0.100, 0.290, 0.640, 0.900, 1.4, 2.5, 5.0, 11.0]

      ; Read in an model LRG spectra, assuming the same wavelengths for all
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
   ; Initialize the colors for the templates as a function of redshift.

   if (ageburst NE ageburst_save OR zmetal NE zmetal_save $
    OR keyword_set(synarr) EQ 0) then begin

      numz = 101
      zarr = 0.01 * findgen(numz)

      synflux = dblarr(5,numz)
      for iz=0L, numz-1 do begin
         print, format='("Z ",i5," of ",i5,a1,$)', $
           iz, numz, string(13b)

         ; Convert redshift to an age, using the WMAP cosmology,
         ; and say that these galaxies formed when the Universe
         ; was AGEBURST Gyr old
         hubble0 = 71. * 3.1558e7 / 3.0856e19 * 1e9 ; units of Gyr^-1
         thisage = lookback(1000., 0.27, 0.73) / hubble0 $
          - lookback(zarr[iz], 0.27, 0.73) / hubble0 - ageburst

         ; Linearly interpolate from the templates in log(age),log(zmetal),
         ; vs. log(flux) space.
         i0 = ((reverse(where(agevec LT thisage)))[0] > 0) < (nage-2)
         j0 = ((reverse(where(metalvec LT zmetal)))[0] > 0) < (nmetal-2)
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

         ; Space equally in log-wavelength
         bigloglam = 3.2000d0 + dindgen(8800) * 1.d-4
         bigwave = 10.d0^bigloglam
         linterp, thiswave, thisflux, bigwave, bigspecflux

         synflux[*,iz] = filter_thru( bigspecflux, $
          waveimg=bigwave, /toair)
;if (iz EQ 0) then begin
; print,synflux[*,iz]/synflux[2,iz]
; stop
;endif
      endfor
      print
   endif

   ;----------
   ; Initialize variables

   ndim = size(pflux, /n_dimen)
   dims = size(pflux, /dimens)
   if (ndim EQ 1) then begin
      nobj = 1
      zfit = 0.
      z_err = 0.
      chi2 = 0.
   endif else begin
      nobj = dims[1]
      zfit = fltarr(nobj)
      z_err = fltarr(nobj)
      chi2 = fltarr(nobj)
   endelse

   ;----------
   ; Loop over each object -- fit redshifts

   numz = n_elements(zarr)

   for iobj=0L, nobj-1 do begin
      print, format='("Object ",i5," of ",i5,a1,$)', $
        iobj, nobj, string(13b)

      chi2arr = dblarr(numz)

      ; Apply AB corrections
      if (keyword_set(abcorrect)) then begin
         thisflux = sdssflux2ab( pflux[*,iobj] )
         thisisig = sqrt( sdssflux2ab(pflux_ivar[*,iobj], /ivar) )
      endif else begin
         thisflux = pflux[*,iobj]
         thisisig = sqrt(pflux_ivar[*,iobj])
      endelse

      ; Apply additional fudge terms to AB corrections
      if (keyword_set(abfudge)) then begin
         thisflux = thisflux * 10.d0^(-abfudge/2.5)
         thisisig = thisisig / 10.d0^(-abfudge/2.5)
      endif

      ; Apply extinction corrections
      if (keyword_set(extinction)) then begin
         thisflux = thisflux * 10.^(0.4*extinction[*,iobj])
         thisisig = thisisig / 10.^(0.4*extinction[*,iobj])
      endif

      ; Insist that we don't use any NaN values
      ibad = where(finite(thisflux) EQ 0 OR finite(thisisig) EQ 0, nbad)
      if (nbad GT 0) then begin
         thisflux[ibad] = 0
         thisisig[ibad] = 0
      endif

      ; Add ADDERR in quadrature
      igood = where(thisisig GT 0, ngood)
      if (ngood GT 0) then begin
         thisisig[igood] = sqrt( 1. / (1./thisisig[igood]^2 $
          + (adderr*(thisflux[igood]>0))^2) )
      endif

      ; Loop over each redshift, and compute the chi^2
      for iz=0L, numz-1 do begin
         chi2arr[iz] = computechi2(thisflux[filterlist], thisisig[filterlist], $
          synflux[filterlist,iz], acoeff=acoeff, dof=dof)
      endfor

      ; Fit a quadratic function to the 3 points surrounding the minimum,
      ; and use that function to estimate the error
      zfit[iobj] = find_nminima(chi2arr, zarr, width=1.5*(zarr[1]-zarr[0]), $
       xerr=xerr1, errcode=errcode, ypeak=ypeak1)
      z_err[iobj] = xerr1 * (errcode EQ 0) + errcode
      chi2[iobj] = ypeak1
   endfor
   print

; Below makes the color-color plots for the photo-z
;gr = 2.5*alog10(synflux[2,*]/synflux[1,*])
;ri = 2.5*alog10(synflux[3,*]/synflux[2,*])
;splot,gr,ri,ps=-4

   return, zfit
end
;------------------------------------------------------------------------------
