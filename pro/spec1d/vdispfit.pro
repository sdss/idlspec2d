; OPTIONAL KEYWORDS:
;   npoly      - Number of polynomial terms to append to eigenspectra;
;                default to none.

; DATA FILES:
;   $IDLSPEC2D_DIR/templates/TEMPLATEFILES
;------------------------------------------------------------------------------
function vdisp_gconv, x, sigma, _EXTRA=EXTRA

   ; Special case for no smoothing
   if sigma eq 0 then return, x

   ksize = round(4*sigma+1) * 2
   xx = findgen(ksize) - ksize/2

   kernel = exp(-xx^2 / (2*sigma^2))
   kernel = kernel / total(kernel)

   sm = convol(x, kernel, _EXTRA=EXTRA)

   return, sm
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

   ;---------------------------------------------------------------------------
   ; If multiple object flux vectors exist, then call this routine recursively.

   if (nobj GT 1) then begin
      sigma = fltarr(nobj)
      sigerr = fltarr(nobj)
      for iobj=0, nobj-1 do begin
         if (keyword_set(objloglam)) then $
          vdispfit, objflux[*,iobj], objivar[*,iobj], objloglam[*,iobj], $
           zobj=zobj[iobj], npoly=npoly, sigma=sigma1, sigerr=sigerr1 $
         else $
          vdispfit, objflux[*,iobj], objivar[*,iobj], hdr=hdr, $
           zobj=zobj[iobj], npoly=npoly, sigma=sigma1, sigerr=sigerr1
         sigma[iobj] = sigma1
         sigerr[iobj] = sigerr1
      endfor
      return
   endif

   ;---------------------------------------------------------------------------

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

   endif

   ;----------
   ; Find the pixel numbers to use from the object and the templates

   ; Find the sub-pixel shifts in the object
   subshift = round(((bigloglam[0]-restloglam[0]) / objdloglam MOD 1) * nsamp)
   indx = subshift + nsamp * lindgen(nbigpix/nsamp)

   if (max(restloglam) LT min(bigloglam[indx]) $
    OR min(restloglam) GT max(bigloglam[indx])) then begin
      splog, 'No wavelength overlap with template'
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

      chi2arr[isig] = computechi2(objsmall, sqivar, eigenflux, $
       acoeff=acoeff, dof=dof, yfit=yfit)

   endfor

   findchi2min, bigsig, chi2arr, minchi2, sigma, sigerr
print,sigma

end
;------------------------------------------------------------------------------
