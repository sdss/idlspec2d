;+
; NAME:
;   lrgmodel_photoz
;
; PURPOSE:
;   Simple photo-z finder for LRGs using Bruzual-Charlot models.
;
; CALLING SEQUENCE:
;   zfit = lrgmodel_photoz(pflux, pflux_ivar, [ /abcorrect, extinction=, $
;    abfudge=, ageburst=, zmetal=, filterlist=, adderr=, $
;    zsplinearr=, synfluxarr=, z_err=, chi2=, fitflux=, plotfile= ] )
;
; INPUTS:
;   pflux          - Object fluxes in the 5 SDSS filters [5,NOBJ]; if not set,
;                    then simply return with the optional ZSPLINEARR,SYNFLUXARR.
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
;                    to value from previous call, or Z=0.018
;                    (where Z=0.02 is solar)
;   filterlist     - List of filter indices to use in fits; default to
;                    using all five filters [0,1,2,3,4]
;   adderr         - Fractional error to add in quadrature; default to 0.03
;   plotfile       - If set, then send debugging plots to this file
;
; OUTPUTS:
;   zfit           - Best-fit redshift [NOBJ]
;
; OPTIONAL OUTPUTS:
;   zsplinearr     - Array of test redshifts used in the fits [NTEST]
;   synfluxarr     - Array of the fluxes vs. redshift from the models [5,NTEST]
;   z_err          - Redshift error, or a negative value if an error
;                    occurred in the quadratic minimization estimate [NOBJ]
;   chi2           - Best-fit chi^2 [NOBJ]
;   fitflux        - Fluxes for the best-fit model, useful for knowing the
;                    colors of the best-fit model
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
;   Note that I allow the Bruzual-Charlot models to be extrapolated
;   beyond the ages and metallicities that they provide.  This is done
;   via linear interpolation in log(age)-log(metals)-log(flux) spec.
;   This may not be valid.
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
 filterlist=filterlist, adderr=adderr, chi2=chi2, synfluxarr=synfluxarr, $
 zsplinearr=zsplinearr, fitflux=fitflux, plotfile=plotfile

   common com_lrgmodel_photoz, zarr, agevec, metalvec, $
    allwave, allflux, synflux, ageburst_save, zmetal_save

   if (n_elements(filterlist) EQ 0) then filterlist = lindgen(5)
   if (n_elements(adderr) EQ 0) then adderr = 0.03
   if (n_elements(abfudge) EQ 0) then abfudge = [0,0,0,0,0]
   if (n_elements(ageburst_save) EQ 0) then ageburst_save = -999
   if (n_elements(zmetal_save) EQ 0) then zmetal_save = -999

   if (keyword_set(ageburst1)) then ageburst = ageburst1 $
    else ageburst = 2.5
   if (keyword_set(zmetal1)) then zmetal = zmetal1 $
    else zmetal = 0.018

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
    OR keyword_set(synflux) EQ 0) then begin

      dz = 0.01
      zmin = -0.03
      zmax = 1.0
      numz = (zmax - zmin) / dz + 1
      zarr = zmin + dz * findgen(numz)

      synflux = dblarr(5,numz)
      for iz=0L, numz-1 do begin
         print, format='("Z ",i5," of ",i5,a1,$)', $
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

         ; Space equally in log-wavelength
         bigloglam = 3.2000d0 + dindgen(8800) * 1.d-4
         bigwave = 10.d0^bigloglam
         linterp, thiswave, thisflux, bigwave, bigspecflux

         synflux[*,iz] = filter_thru( bigspecflux, $
          waveimg=bigwave, /toair)
      endfor
      print
      ageburst_save = ageburst
      zmetal_save = zmetal
  endif

  if (arg_present(zsplinearr)) then zsplinearr = zarr
  if (arg_present(synfluxarr)) then synfluxarr = synflux
  if (NOT keyword_set(pflux)) then return, 0

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
      print, format='("Object ",i8," of ",i8,a1,$)', $
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
if (finite(chi2arr[iz]) EQ 0) then stop
      endfor

      ; Fit a quadratic function to the 3 points surrounding the minimum,
      ; and use that function to estimate the error
      zfit[iobj] = find_nminima(chi2arr, zarr, width=1.5*(zarr[1]-zarr[0]), $
       xerr=xerr1, errcode=errcode, ypeak=ypeak1)
      z_err[iobj] = xerr1 * (errcode EQ 0) + errcode
      chi2[iobj] = ypeak1

      if (arg_present(fitflux)) then begin
         ; Interpolate the model fluxes to the exact redshift
         if (iobj EQ 0) then fitflux = fltarr(5,nobj)
         for j=0, 4 do begin
            linterp, zarr, transpose(synflux[j,*]), zfit[iobj], flux1
            fitflux[j,iobj] = flux1
         endfor
         ; Re-fit the amplitude of the model fluxes to match the object
         junk = computechi2(thisflux[filterlist], thisisig[filterlist], $
          fitflux[filterlist,iobj], acoeff=acoeff)
         fitflux[filterlist,iobj] = acoeff[0] * fitflux[filterlist,iobj]
       endif

   endfor
   print

   ; Below makes the color-color plots for the photo-z
   if (keyword_set(plotfile)) then begin
      dfpsplot, plotfile, /square, /color
      for iobj=0L, nobj-1 do begin
         for j=0, 2 do begin
            color1 = (2.5*alog10(synflux[j+1,*]/synflux[j,*]))[*]
            color2 = (2.5*alog10(synflux[j+2,*]/synflux[j+1,*]))[*]
            thiscolor1 = 2.5*alog10(pflux[j+1,iobj]/pflux[j,iobj])
            thiscolor2 = 2.5*alog10(pflux[j+2,iobj]/pflux[j+1,iobj])
            thiserr1 = 2.5*alog10(1+1./(pflux[j,iobj]*sqrt(pflux_ivar[j,iobj]))) $
             + 2.5*alog10(1+1./(pflux[j+1,iobj]*sqrt(pflux_ivar[j+1,iobj])))
            thiserr2 = 2.5*alog10(1+1./(pflux[j+1,iobj]*sqrt(pflux_ivar[j+1,iobj]))) $
             + 2.5*alog10(1+1./(pflux[j+2,iobj]*sqrt(pflux_ivar[j+2,iobj])))
            djs_plot, [color1,thiscolor1], [color2,thiscolor2], /nodata, $
             xtitle='('+filtername(j)+'-'+filtername(j+1)+')', $
             ytitle='('+filtername(j+1)+'-'+filtername(j+2)+')', charsize=1.5, $
             title='Object '+strtrim(iobj,2) $
             + ' z(FIT)='+strtrim(string(zfit[iobj],format='(f5.2)'),2)
            djs_oplot, color1, color2, psym=-4
            for k=0, (n_elements(zarr)-1)/10 do $
             djs_xyouts, color1[k*10], color2[k*10], $
              '  z='+strtrim(string(zarr[k*10], format='(f5.2)'),2), $
              charsize=1.5
            djs_oplot, [thiscolor1], [thiscolor2], psym=6, symsize=2, color='red'
            djs_oplot, thiscolor1+[-1,1]*thiserr1, thiscolor2+[0,0], color='red'
            djs_oplot, thiscolor1+[0,0], thiscolor2+[-1,1]*thiserr2, color='red'
         endfor
      endfor
      dfpsclose
   endif

   return, zfit
end
;------------------------------------------------------------------------------
