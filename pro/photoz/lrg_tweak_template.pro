;+
; NAME:
;   lrg_tweak_template
;
; PURPOSE:
;   Compute the corrections to the input LRG template used for LRG_PHOTOZ().
;
; CALLING SEQUENCE:
;   lrg_tweak_template, pflux, pflux_ivar, zz, $
;    [ /abcorrect, extinction=, abfudge=, filterlist=, adderr=, $
;    maxiter=, coeff= ]
;
; INPUTS:
;   pflux          - Object fluxes in the 5 SDSS filters [5,NOBJ]
;   pflux_ivar     - Inverse variances for FLUX [5,NOBJ]
;   zz             - Spectroscopic redshifts [NOBJ]
;
; OPTIONAL INPUTS:
;   abcorrect      - If set, then convert the input fluxes from the 2.5-m
;                    natural system to AB fluxes
;   extinction     - If set, then apply these extinction corrections [5,NOBJ]
;   abfudge        - Additional AB "fudge factors"; default to adding
;                    [0,-0.03,0,0,0] mag to input magnitudes, where a positive
;                    value makes that flux fainter
;   filterlist     - List of filter indices to use in fits; default to
;                    using all five filters [0,1,2,3,4]
;   adderr         - Fractional error to add in quadrature; default to 0.03
;   maxiter        - Maximum number of iterations; default to 40
;
; OUTPUTS:
;   coeff          - Best-fit coefficients for tweaking input LRG template
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine computes the modifications to an input LRG spectrum
;   that gives a better total chi^2 when computing photometric redshifts.
;   The form of the transformation is:
;      Flux = Flux * Wave^COEFF[0] * (1 + COEFF[1] * z * Wave)
;
;   The fluxes should be AB fluxes, or SDSS 2.5-m natural system fluxes
;   if /ABCORRECT is set.
;   The fluxes should already be extinction-corrected, unless
;   the EXTINCTION keyword is passed.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   computechi2()
;   filter_thru()
;   mrdfits()
;   sdssflux2ab
;   sxpar()
;
; REVISION HISTORY:
;   21-Oct-2003  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
forward_function mpfit, lrg_tweak_fn

;------------------------------------------------------------------------------
pro lrg_tweak_template_read, retwave, retflux, coeff=coeff

   common com_lrg_tweak_template_read, loglam, fileflux

   if (NOT keyword_set(loglam)) then begin
      ;----------
      ; Initialize the template "spectra".
      ; This need only be done the first time this function is called,
      ; then cached for future calls.

      ; Read in an LRG spectrum
      eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
      eigenfile = 'spLRG-52928.fits'
      fileflux = mrdfits(djs_filepath(eigenfile, root_dir=eigendir), 0, hdr)
      dloglam = sxpar(hdr, 'COEFF1')
      loglam = sxpar(hdr, 'COEFF0') + dindgen(sxpar(hdr, 'NAXIS1')) * dloglam
   endif

   specflux = fileflux * 10.d0^(loglam*coeff[0])

   ; Smooth the end of the spectra
   i1 = where(loglam LT alog10(3200), n1)
   i2 = where(loglam GT alog10(8800), n2)
   val1 = median(specflux[i1])
   val2 = median(specflux[i2])
   w1 = findgen(n1) / (n1-1)
   w2 = (n2-findgen(n2)-1) / (n2-1)
   specflux[i1] = val1 * (1-w1) + specflux[i1] * w1
   specflux[i2] = val2 * (1-w2) + specflux[i2] * w2

   ; Expand the wavelength coverage by extending the ends of the spectrum
   ; to about 1584,12019 Ang.
   bigloglam = 3.2000d0 + dindgen(8800) * 1.d-4
   linterp, loglam, specflux, bigloglam, bigspecflux

   retwave = 10.d0^bigloglam
   retflux = bigspecflux

   return
end
;------------------------------------------------------------------------------
; Return a vector of all the chi values
function lrg_tweak_fn, coeff

   common com_lrg_tweak_fluxes, pflux, pflux_isig, zz, filterlist

   ; Read the tweaked template
   lrg_tweak_template_read, bigwave, bigspecflux, coeff=coeff

   ; Convert from f_lambda to f_nu
   flambda2fnu = bigwave^2 / 2.99792d18
   bigspecflux = bigspecflux * flambda2fnu

   ; Select a binning of 0.005 in redshift, and go out to redshift 0.65
   numz = 131
   deltaz = 0.005
   zarr = deltaz * findgen(numz)

   ; Convolve the spectrum with the filters
   synflux = dblarr(5,numz)
   for iz=0L, numz-1 do begin
      print, format='("Z ",i5," of ",i5,a1,$)', $
        iz, numz, string(13b)
      thiswave = bigwave * (1 + zarr[iz])
      synflux[*,iz] = filter_thru( $
       bigspecflux * (1.d0 + coeff[1] * zarr[iz] * bigwave), $
       waveimg=thiswave, /toair)
   endfor
   print

   ; Compute the chi for each object
   nobj = n_elements(pflux) / 5
   chivec = dblarr(nobj)
   for iobj=0L, nobj-1 do begin
      iz = (round(zz[iobj] / deltaz) > 0) < (numz-1)
      chivec[iobj] = computechi2(pflux[filterlist,iobj], $
       pflux_isig[filterlist,iobj], synflux[filterlist,iz])
   endfor

   return, chivec
end
;------------------------------------------------------------------------------
pro lrg_tweak_template, pflux1, pflux_ivar1, zz1, $
 abcorrect=abcorrect, extinction=extinction, abfudge=abfudge, $
 filterlist=filterlist1, adderr=adderr, maxiter=maxiter

   common com_lrg_tweak_fluxes, pflux, pflux_isig, zz, filterlist

   ; Set defaults
   if (n_elements(abfudge) EQ 0) then abfudge = [0,-0.03,0,0,0]
   if (n_elements(filterlist1) NE 0) then filterlist = filterlist1 $
    else filterlist = lindgen(5)
   if (n_elements(adderr) EQ 0) then adderr = 0.03
   if (NOT keyword_set(maxiter)) then maxiter = 40

   ; Set variables in common blocks
   zz = zz1

   ; Apply AB correction
   if (keyword_set(abcorrect)) then begin
      pflux = sdssflux2ab(pflux1)
      pflux_isig = sqrt( sdssflux2ab(pflux_ivar1, /ivar) )
   endif else begin
      pflux = pflux1
      pflux_isig = sqrt(pflux_ivar1)
   endelse

   ; Apply additional fudge terms to AB corrections
   if (keyword_set(abfudge)) then begin
      for jfilt=0, 4 do begin
         pflux[jfilt,*] = pflux[jfilt,*] * 10.d0^(-abfudge[jfilt]/2.5)
         pflux_isig[jfilt,*] = pflux_isig[jfilt,*] / 10.d0^(-abfudge[jfilt]/2.5)
      endfor
   endif

   ; Apply extinction corrections
   if (keyword_set(extinction)) then begin
      pflux = pflux * 10.^(0.4*extinction)
      pflux_isig = pflux_isig / 10.^(0.4*extinction)
   endif

   ; Insist that we don't use any NaN values
   ibad = where(finite(pflux) EQ 0 OR finite(pflux_isig) EQ 0, nbad)
   if (nbad GT 0) then begin
      pflux[ibad] = 0
      pflux_isig[ibad] = 0
   endif

   ; Add ADDERR in quadrature
   igood = where(pflux_isig GT 0, ngood)
   if (ngood GT 0) then begin
      pflux_isig[igood] = sqrt( 1. / (1./pflux_isig[igood]^2 $
       + (adderr*(pflux[igood]>0))^2) )
   endif

   ; Discard any objects where the baseline photo-z is discrepent by
   ; more than 0.10, and discard any low-redshift objects with z > 0.10.
   zfit = lrg_photoz(pflux, pflux_isig^2, abfudge=0, adderr=0)
   ibad = where(abs(zfit - zz) GT 0.10 OR zz LT 0.10, nbad)
   if (nbad GT 0) then begin
      print, 'Discard ', nbad, ' objects with low and/or discrepent photo-z'
      pflux[*,ibad] = 0
      pflux_isig[*,ibad] = 0
   endif

   ; Call MPFIT to iterate on the solution for the template
   parinfo = {value: 0.D, fixed: 0, limited: [0b,0b], limits: [0.d0,0.d0]}
   parinfo = replicate(parinfo, 2)
   parinfo.value = [0.24d0, -3.d-5]
   ftol = 1d-20
   gtol = 1d-20
   xtol = 1d-20
   coeff = mpfit('lrg_tweak_fn', parinfo=parinfo, perror=perror, $
    maxiter=maxiter, ftol=ftol, gtol=gtol, xtol=xtol, $
    niter=niter, status=status)

   print, 'STATUS = ', status
   print, 'Best-fit coeffs = ', coeff
   print, 'Errors = = ', perror

stop
   return
end
;------------------------------------------------------------------------------
