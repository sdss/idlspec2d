;+
; NAME:
;   vdispfit
;
; PURPOSE:
;   Compute velocity dispersions for galaxy spectra.
;
; CALLING SEQUENCE:
;   vdispfit, objflux, objivar, [ objloglam, hdr=, zobj=, npoly=, $
;    sigma=, sigerr= ]
;
; INPUTS:
;   objflux    - Galaxy spectrum (spectra); array of [NPIX,NGALAXY].
;   objivar    - Galaxy inverse variance; array of [NPIX,NGALAXY].
;
; OPTIONAL INPUTS:
;   objloglam  - Log-10 wavelengths; this can be either an NPIX vector
;                if all the galaxy spectra have the same wavelength mapping,
;                or an array with the same dimensions as OBJFLUX.
;                Either OBJLOGLAM or HDR must be specified.
;   hdr        - FITS header from which to read COEFF0, COEFF1 for the
;                wavelength mapping.
;                Either OBJLOGLAM or HDR must be specified.
;   zobj       - Redshift for each galaxy; default to 0.
;   npoly      - Number of polynomial terms to append to eigenspectra;
;                default to 5.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   sigma      - Velocity dispersion in km/sec.
;   sigerr     - Error for SIGMA in km/sec.
;
; COMMENTS:
;   Note that the wavelength spacing in the galaxy and stellar template spectra
;   must be the same.
;
;   We currently mask within +/- 280 km/sec of the following wavelengths
;   that could have emission lines:
;     linelist = [3725.94, 3727.24, 3970.072, 4101.73, 4340.46, $
;      4861.3632, 4958.911, 5006.843, 6300.32, 6548.05, 6562.801, $
;      6583.45, 6716.44, 6730.82]
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/templates/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   airtovac
;   combine1fiber
;   computechi2()
;   djs_filepath()
;   findchi2min
;   mrdfits()
;   poly_array()
;   splog
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;   vdisp_gconv()
;
; REVISION HISTORY:
;   13-Mar-2001  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function vdisp_gconv, x, sigma, _EXTRA=EXTRA

   ; Special case for no smoothing
   if (sigma EQ 0) then return, x

   ksize = round(4*sigma+1) * 2
   xx = findgen(ksize) - ksize/2

   kernel = exp(-xx^2 / (2*sigma^2))
   kernel = kernel / total(kernel)

   return, convol(x, kernel, _EXTRA=EXTRA)
end

;------------------------------------------------------------------------------
pro vdispfit, objflux, objivar, objloglam, hdr=hdr, zobj=zobj, npoly=npoly, $
 sigma=sigma, sigerr=sigerr

   common com_vdispfit, bigflux, bigloglam, bigmask, nsamp, bigsig, nbigpix, nsig

   if (NOT keyword_set(objloglam) AND NOT keyword_set(hdr)) then $
    message, 'Must specify either OBJLOGLAM or HDR!'

   if (n_elements(npoly) EQ 0) then npoly = 5
   dims = size(objflux, /dimens)
   npixobj = dims[0]
   if (size(objflux, /n_dimen) EQ 1) then nobj = 1 $
    else nobj = dims[1]
   if (NOT keyword_set(zobj)) then zobj = fltarr(nobj)

   ;---------------------------------------------------------------------------
   ; If multiple object flux vectors exist, then call this routine recursively.

   if (nobj GT 1) then begin
      sigma = fltarr(nobj)
      sigerr = fltarr(nobj)
      lamdims = size(objloglam, /n_dimens)
      for iobj=0, nobj-1 do begin
         if (lamdims EQ 1) then thisloglam = objloglam $
          else if (lamdims EQ 2) then thisloglam = objloglam[*,iobj]
         vdispfit, objflux[*,iobj], objivar[*,iobj], thisloglam, hdr=hdr, $
          zobj=zobj[iobj], npoly=npoly, sigma=sigma1, sigerr=sigerr1
         sigma[iobj] = sigma1
         sigerr[iobj] = sigerr1
      endfor
      return
   endif

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   if (NOT keyword_set(objloglam)) then begin
      objloglam0 = sxpar(hdr, 'COEFF0')
      objdloglam = sxpar(hdr, 'COEFF1')
      objloglam = objloglam0 + dindgen(npixobj) * objdloglam
   endif else begin
      objdloglam = objloglam[1] - objloglam[0]
   endelse
   restloglam = objloglam - alog10(1 + zobj) ; De-redshift this!

   ;---------------------------------------------------------------------------
   ; Generate the over-sampled eigen-templates for the stellar spectra.
   ; This is saved in a common block between calls.

   if (NOT keyword_set(bigflux)) then begin

      nsamp = 10
      nsig = 20
      bigsig = findgen(nsig) * 25.0 ; in km/sec

      ;----------
      ; Find the most recent template file matching EIGENFILE

      eigenfile = 'spEigenVdisp*.fits'
      eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
      allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
      if (ct EQ 0) then $
       message, 'Unable to find EIGENFILE matching '+eigenfile
      thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]
      splog, 'Selecting EIGENFILE=' + thisfile

      ;----------
      ; Read the stellar templates

      eflux = mrdfits(thisfile, 0, ehdr)
      naxis1 = sxpar(ehdr,'NAXIS1')
      nstar = sxpar(ehdr,'NAXIS2') > 1
      eloglam0 = sxpar(ehdr, 'COEFF0')
      edloglam = sxpar(ehdr, 'COEFF1')
      eloglam = eloglam0 + dindgen(naxis1) * edloglam

      ; Pixel size in km/sec for these oversampled (smaller) pixels
      cspeed = 2.99792458e5
      pixsz = (10.^(edloglam)-1) * cspeed / nsamp

      ;----------
      ; Re-samples to higher resolution by a factor of NSAMP

      nbigpix = (naxis1 - 1) * nsamp + 1
      bigloglam = eloglam0 + dindgen(nbigpix) * edloglam / nsamp
      bigflux = fltarr(nbigpix, nstar, nsig)

      for istar=0, nstar-1 do begin
         combine1fiber, eloglam, eflux[*,istar], $
          newloglam=bigloglam, newflux=tmpflux, maxiter=0
         bigflux[*,istar,0] = tmpflux
         if (istar EQ 0) then bigmask = tmpflux NE 0
      endfor

      ;----------
      ; Generate array of broadened templates

      for isig=1, nsig-1 do begin
         for istar=0, nstar-1 do begin
            bigflux[*,istar,isig] = $
             vdisp_gconv(bigflux[*,istar,0], bigsig[isig]/pixsz, /edge_truncate)
         endfor
      endfor

      ;----------
      ; Mask out the first few and last few pixels, since those would not
      ; have been smoothed properly.  The masked region is 350 km/sec
      ; for a native binning of 70 km/sec.

      bigmask[0:4*nsamp-1] = 0
      bigmask[nbigpix-4*nsamp:nbigpix-1] = 0

      ;----------
      ; Mask out emission lines, setting bigmask=0 near these wavelengths.

      linelist = [3725.94, 3727.24, 3970.072, 4101.73, 4340.46, $
       4861.3632, 4958.911, 5006.843, 6300.32, 6548.05, 6562.801, $
       6583.45, 6716.44, 6730.82]
      vaclist = linelist
      airtovac, vaclist
      vaclist = alog10(vaclist)
      mwidth = 4.e-4 ; Mask out any pixels within +/- 280 km/s
      for iline=0, n_elements(vaclist)-1 do $
       bigmask = bigmask AND (bigloglam LT vaclist[iline] - mwidth $
        OR bigloglam GT vaclist[iline] + mwidth)

   endif

   ;----------
   ; Find the pixel numbers to use from the object and the templates

   ; Find the sub-pixel shifts in the object
   subshift = round(((bigloglam[0]-restloglam[0]) / objdloglam MOD 1) * nsamp)
   indx = subshift + nsamp * lindgen(nbigpix/nsamp)

   if (max(restloglam) LT min(bigloglam[indx]) $
    OR min(restloglam) GT max(bigloglam[indx])) then begin
;      splog, 'No wavelength overlap with template'
      sigma = 0.0
      sigerr = 9999.
      return
   endif

   if (restloglam[0] LT bigloglam[indx[0]]) then begin
      ipixt0 = 0L
      junk = min(abs(restloglam - bigloglam[indx[0]]), ipixo0)
   endif else begin
      ipixo0 = 0L
      junk = min(abs(bigloglam[indx] - restloglam[0]), ipixt0)
   endelse

   npixcomp = (npixobj - ipixo0 + 1) < (n_elements(indx) - ipixt0)

   indxo = ipixo0 + lindgen(npixcomp) ; Indices for object spectrum
   indxt = indx[ipixt0 + lindgen(npixcomp)] ; Indices for template spectra

   ;----------
   ; Add more eigen-templates that represent polynomial terms.

   if (keyword_set(npoly)) then $
    polyflux = poly_array(npixcomp,npoly)

   ;----------
   ; Fit for chi^2 at each possible velocity dispersion

   chi2arr = fltarr(nsig)

   objsmall = objflux[indxo]
   sqivar = sqrt( objivar[indxo] ) * bigmask[indxt]

   for isig=0, nsig-1 do begin

      eigenflux = bigflux[indxt,*,isig]
      if (keyword_set(npoly)) then eigenflux = [[eigenflux], [polyflux]]

      chi2arr[isig] = computechi2(objsmall, sqivar, eigenflux)

   endfor

   ;----------
   ; Fit for the dispersion value at the minimum in chi^2

   findchi2min, bigsig, chi2arr, minchi2, sigma, sigerr

   ;----------
   ; If the best-fit value is at the maximum dispersion value tested,
   ; then we don't really know the answer and should set the error
   ; to a large value.

   if (sigma GE max(bigsig)) then begin
      sigerr = 9999.
   endif

   return
end
;------------------------------------------------------------------------------
