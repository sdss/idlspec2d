function zqso_nolya, objflux, objivar, objloglam0, objdloglam, $
 npoly=npoly, pwidth=pwidth, eigenfile=eigenfile, zold=zold

   ndim = size(objflux, /n_dimen)
   if (ndim EQ 1) then nobj = 1 $
    else nobj = (size(objflux, /dimens))[1]

   if (n_elements(zold) GT 0 AND keyword_set(pwidth)) then begin
      pmin = floor( alog10(1.0 + zold.z) / objdloglam - 0.5 * pwidth )
      pmax = floor( alog10(1.0 + zold.z) / objdloglam + 0.5 * pwidth )
   endif

   eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates/')
   if NOT keyword_set(eigenfile) then begin
      if n_elements(zold) GT 0 then eigenfile=zold[0].tfile $
      else message, 'need eigenfile or zold'
   endif

   efile = strtrim(eigendir+eigenfile,2)

   starflux = readfits(efile, shdr)
   starloglam0 = sxpar(shdr, 'COEFF0')
   stardloglam = sxpar(shdr, 'COEFF1')

   ndim = size(starflux, /n_dimen)
   dims = size(starflux, /dimens)
   npixstar = dims[0]
   if (ndim EQ 1) then nstar = 1 $
    else nstar = dims[1]

   starlogwave = findgen(npixstar) * stardloglam + starloglam0
   starmask = starlogwave GT alog10(1270.0)

   ;----------
   ; Add more eigen-templates that represent polynomial terms.

   if (keyword_set(npoly)) then $
    starflux = [ [starflux], [poly_array(npixstar,npoly)] ]

   ;----------
   ; Compute the redshift difference between the first pixel of the object
   ; spectra and the template.

   poffset = (objloglam0 - starloglam0) / objdloglam

   ;----------
   ; Compute the redshifts

   zans = zcompute(objflux, objivar, starflux, starmask, poffset=poffset, $
    pmin=pmin, pmax=pmax, pspace=1, width=101, /doplot)

   ;----------
   ; Convert redshift (and error) from pixels to the conventional dimensionless
   ; value.

   indx = where(zans.dof GT 0)
   if (indx[0] NE -1) then begin
      zans[indx].z = 10.^(objdloglam[indx] * zans[indx].z) - 1.
      zans[indx].z_err = $
       alog(10d) * objdloglam[indx] * zans[indx].z_err * (1 + zans[indx].z)
   endif

   return, zans
end
;------------------------------------------------------------------------------
