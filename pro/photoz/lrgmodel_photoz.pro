;+
; NAME:
;   lrgmodel_photoz
;
; PURPOSE:
;   Simple photo-z finder for LRGs using Bruzual-Charlot models.
;
; CALLING SEQUENCE:
;   zfit = lrgmodel_photoz(pflux, pflux_ivar, [ /abcorrect, extinction=, $
;    abfudge=, filterlist=, adderr=, z_err=, chi2= ] )
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
 filterlist=filterlist, adderr=adderr, chi2=chi2

   common com_lrgmodel_photoz, zarr, synflux

   if (n_elements(filterlist) EQ 0) then filterlist = lindgen(5)
   if (n_elements(adderr) EQ 0) then adderr = 0.03
   if (n_elements(abfudge) EQ 0) then abfudge = [0,0.0,0,0,0]

   ;----------
   ; Initialize the template "spectra".
   ; This need only be done the first time this function is called,
   ; then cached for future calls.

   if (NOT keyword_set(zarr)) then begin

      metals = 'z02' ; ???
      eigenfile = [ $
       'ssp_5Myr_'+metals+'.spec', $
       'ssp_25Myr_'+metals+'.spec', $
       'ssp_100Myr_'+metals+'.spec', $
       'ssp_290Myr_'+metals+'.spec', $
       'ssp_640Myr_'+metals+'.spec', $
       'ssp_900Myr_'+metals+'.spec', $
       'ssp_1.4Gyr_'+metals+'.spec', $
       'ssp_2.5Gyr_'+metals+'.spec', $
       'ssp_5Gyr_'+metals+'.spec', $
       'ssp_11Gyr_'+metals+'.spec' ]
      nfile = n_elements(eigenfile)

      ; Age of each template in Gyr
      ages = [0.005, 0.025, 0.100, 0.290, 0.640, 0.900, 1.4, 2.5, 5.0, 11.0]

      ; Read in an model LRG spectra, assuming the same wavelengths for all
      eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
      for ifile=0L, nfile-1 do begin
         readcol, djs_filepath(eigenfile[ifile], root_dir=eigendir), $
          lambda1, flux1, comment='#', format='(F,F)'
         if (ifile EQ 0) then begin
            npix = n_elements(lambda1)
            allwave = lambda1
            allflux = fltarr(npix, nfile)
         endif
         ; Convert to f_nu
         flambda2fnu = lambda1^2 / 2.99792d18
         allflux[*,ifile] = flux1 * flambda2fnu
      endfor

      numz = 101
      zarr = 0.01 * findgen(numz)

      synflux = dblarr(5,numz)
      for iz=0L, numz-1 do begin
         print, format='("Z ",i5," of ",i5,a1,$)', $
           iz, numz, string(13b)

         ; Convert redshift to an age, using the WMAP cosmology,
         ; and say that these galaxies formed when the Universe
         ; was AGEBURST Gyr old
ageburst = 2.5 ; in Gyr ???
         hubble0 = 71. * 3.1558e7 / 3.0856e19 * 1e9 ; units of Gyr^-1
         thisage = lookback(1000., 0.27, 0.73) / hubble0 $
          - lookback(zarr[iz], 0.27, 0.73) / hubble0 - ageburst

         ; Linearly interpolate from the templates in log(age) vs. log(flux)
         ; space.
         j0 = (reverse(where(ages LT thisage)))[0]
         if (j0 EQ nfile-1) then begin
            thisflux = allflux[*,j0]
         endif else begin
            j1 = j0 + 1
            wts = [alog10(ages[j1]/thisage), -alog10(ages[j0]/thisage)]
            wts = wts / total(wts)
            thisflux = 10.^( wts[0] * alog10(allflux[*,j0]) $
             + wts[1] * alog10(allflux[*,j1]))
         endelse

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
