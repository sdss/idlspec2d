;+ 
; NAME:
;   spflux_read_desi
;
; PURPOSE:
;   Read the Kurucz models at the specified log-lambda and resolution
;
; CALLING SEQUENCE:
;   spflux_read_desi
;   modelflux = spflux_read_desi( loglam, dispimg, [ iselect=, $
;    kindx_return= ] )
;
; INPUTS:
;   loglam     - Log10 wavelengths (vacuum Angstroms) [NPIX]
;   dispimg    - Dispersion image, in units of pixels [NPIX]
;
; OPTIONAL INPUTS:
;   iselect    - If set, then only return these model numbers; default to
;                returning all models
;
; OUTPUTS:
;   modelflux  - Model fluxes [NPIX,NMODEL]
;
; OPTIONAL OUTPUTS:
;   kindx_return- Structure with model parameters for each model
;
; COMMENTS:
;   Note that this function allocates a ridiculous amount of memory
;   for caching the oversampled spectra at many dispersions in GRIDFLUX.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/templates/star_templates_v2.2.fits
;
function spflux_read_desi, loglam, dispimg, iselect=iselect1, $
 kindx_return=kindx_return

   common com_spflux_kurucz, kfile, kflux, kindx, kloglam, nmodel, $
    gridsig, gridlam, gridflux

   if (n_elements(loglam) NE n_elements(dispimg)) then $
    message, 'LOGLAM and DISPIMG must have same number of elements!'

   ;----------
   ; Read the DESI stellar models

   if (NOT keyword_set(kfile)) then begin
      ; Read the file with the desi models as generated
      ; The units should already be erg/cm^2/s/Ang.
      splog, 'Reading DESI stellar models'
      kfile = filepath('desi_stellar_templates_v2.2_trimmed.fits', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
      krawflux = mrdfits(kfile, 0, hdr, /silent)  ; flux
      kindx = mrdfits(kfile, 1, /silent)
      dims = size(krawflux, /dimens)
      nrawpix = dims[0]
      nmodel = dims[1]
      waves = mrdfits(kfile, 2, /silent) 

      ; The models will be sub-sampled relative to the SDSS pix scale of 1d-4:
      subsamp = 5
      kdlog = 1.d-4 / subsamp

      ;-- JEB 2018-05 desi templates are already in vacuum wavelength
      ; These models are sampled linearly in **air** wavelength.
      ; Re-map them to linear in (vacuum) log-wavelength.
      ; airtovac, waves ; Remap wavelengths from air to vacuum
      rawloglam = alog10(waves)
      splog, 'Remapping models to log-wavelengths'
      minlog1 = min(rawloglam, max=maxlog1)
      kloglam = wavevector(minlog1, maxlog1, binsz=kdlog)
      npix = n_elements(kloglam)
      kflux = fltarr(npix,nmodel)
      for imodel=0L, nmodel-1 do begin
       splog, ' rebinning model '+strtrim(imodel,2)
       kflux[*,imodel] = rebin_spectrum(krawflux[*,imodel], rawloglam, kloglam)
      endfor
   
      ;-- JEB 2018-05 removing this step since desi templates are in loglam already
      ; The above rebinning conserves flux, so we are no longer in
      ; units of erg/cm^2/s/Ang.  Divide by one power of wavelength
      ; to recover those units.  (The -4. is an arbitrary scaling.)
      ;for imodel=0L, nmodel-1 do $
      ; kflux[*,imodel] = kflux[*,imodel] / 10.d0^(kloglam-4.)

      
      ;-- JEB 2018-05 removing convolution
      ; Convolve the Kurucz models with a boxcar response
      ; representing the size of the SDSS pixels.
      ;splog, 'Convolving Kurucz models with SDSS pixel size'
      ;if (subsamp GT 1) then begin
      ;   if ((subsamp MOD 2) EQ 0) then $
      ;    kern = ([0.5, fltarr(subsamp-1) + 1.0, 0.5]) / subsamp $
      ;   else $
      ;    kern = (fltarr(subsamp) + 1.0) / subsamp
      ;   for imodel=0L, nmodel-1 do begin
      ;     kflux[*,imodel] = convol(kflux[*,imodel], kern, /center)
      ;   endfor
      ;endif

      ; Compute the spectro-photo fluxes for these spectra
      ; by just using these high-resolution redshift=0 spectra.
      ; Formally, I should re-compute these fluxes after redshifting
      ; and convolving with the dispersion, but I'm choosing to
      ; ignore those tiny effects.
;      wavevec = 10.d0^kloglam
;      flambda2fnu = wavevec^2 / 2.99792e18
;      fthru = filter_thru(kflux * rebin(flambda2fnu,npix,nmodel), $
;       waveimg=wavevec, /toair)
;      mag = - 2.5 * alog10( fthru * 10^((48.6 - 2.5*17.)/2.5) )
;      kindx.mag = transpose(mag)

      ; Construct a grid of these models at various possible resolutions
      splog, 'Convolving Kurucz models with dispersions'
      gridsig = 0.70 + findgen(21) * 0.05 ; span range [0.7,1.7]
      nkpix = 7 * subsamp + 1
      nres = n_elements(gridsig)
      gridlam = kloglam ; Keep the same wavelength mapping
      gridflux = fltarr(npix, nres, nmodel)
      for ires=0L, nres-1 do begin
;         print, format='("Grid ",i5," of ",i5,a1,$)', ires, nres, string(13b)
         kern = exp(-0.5 * (findgen(nkpix*2+1) - nkpix)^2 $
          / (gridsig[ires]*subsamp)^2)
         kern = kern / total(kern)
         for imodel=0L, nmodel-1 do begin
            gridflux[*,ires,imodel] = convol(kflux[*,imodel], kern, /center)
         endfor
      endfor
      print
   endif

   ;----------
   ; Return just with KINDX if LOGLAM,DISPIMG are not set

   if (n_elements(iselect1) GT 0) then iselect = iselect1 $
    else iselect = lindgen(nmodel)
   nselect = n_elements(iselect)
   if (arg_present(kindx_return)) then kindx_return = kindx[iselect]
   if (NOT keyword_set(loglam)) then return, 0

   ;----------
   ; Map the wavelength values and dispersion values onto pixel numbers
   ; on the grid of models.

   xpos = (loglam - gridlam[0]) / (gridlam[1] - gridlam[0])
   xpos = (xpos > 0) < (n_elements(gridlam)-1)
   ypos = (dispimg - gridsig[0]) / (gridsig[1] - gridsig[0])
   ypos = (ypos > 0) < (n_elements(gridsig)-1)

   ;----------
   ; Interpolate within these models

   ndim = size(loglam, /n_dimen)
   dims = size(loglam, /dimens)
   npix = dims[0]
   if (ndim EQ 1) then nobj = 1 $
    else nobj = dims[1]
   modflux = fltarr(npix, nobj, nselect)
   for iobj=0L, nobj-1 do begin
      for j=0L, nselect-1 do begin
         imodel = iselect[j]
         modflux[*,iobj,j] = interpolate(gridflux[*,*,imodel], $
          xpos[*,iobj], ypos[*,iobj])
      endfor
   endfor

   return, modflux
end
;------------------------------------------------------------------------------



