
; zciv = zqso_nolya(flux, invvar, loglam, pwidth=400, zold=zans)
function zqso_nolya, objflux, objivar, loglam, $
 npoly=npoly, pwidth=pwidth, eigenfile=eigenfile, zold=zold, $
 minrest=minrest, maxrest=maxrest, synthspec=synthspec

  if NOT keyword_set(minrest) then minrest=1170.0
  if NOT keyword_set(maxrest) then maxrest=1600.0

   ndim = size(objflux, /n_dimen)
   if (ndim EQ 1) then begin
     nobj = 1
     npix = n_elements(objflux)
   endif else begin
     nobj = (size(objflux, /dimens))[1]
     npix = (size(objflux, /dimens))[0]
   endelse

   objdloglam = loglam[1,*] - loglam[0,*]

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


   rlam = loglam - replicate(1,npix) # alog10(1.0+zold.z)
   objmask = rlam GT alog10(minrest) AND rlam LT alog10(maxrest)

   ;----------
   ; Add more eigen-templates that represent polynomial terms.

   if (keyword_set(npoly)) then $
    starflux = [ [starflux], [poly_array(npixstar,npoly)] ]

   ;----------
   ; Compute the redshift difference between the first pixel of the object
   ; spectra and the template.

   poffset = ((loglam[0,*])[*] - starloglam0) / objdloglam

   ;----------------
   ; Compute the redshifts

   smoothflux = objflux*0.0
   tempivar = objivar*objmask

;
;	 Need a subroutine for good continuum fitting, again?!?
;
;

   for i=0,nobj-1 do begin  
     qso_continuum_fit, loglam[*,i], objflux[*,i], objivar[*,i], $
       model=yfit
     smoothflux[*,i] = yfit
   endfor

   zans = zcompute(smoothflux, tempivar, starflux, poffset=poffset, $
    pmin=pmin, pmax=pmax, pspace=1, width=21)

   ;----------------
   ;  Need to do at least one round of rejections
   
   ;----------
   ; Convert redshift (and error) from pixels to the conventional dimensionless
   ; value.

   indx = where(zans.dof GT 0)
   if (indx[0] NE -1) then begin
      zans[indx].z = 10.^(objdloglam[indx] * zans[indx].z) - 1.
      zans[indx].z_err = $
       alog(10d) * objdloglam[indx] * zans[indx].z_err * (1 + zans[indx].z)
   endif


   zreturn = zold
   zreturn.z = (zans.z)[*]
   zreturn.z_err = (zans.z_err)[*]
   zreturn.tcolumn = -1
   zreturn.tcolumn[0] = 0
   zreturn.theta = 0
   zreturn.theta[0] = (zans.theta)[*]
   zreturn.dof   = (zans.dof)[*]
   zreturn.rchi2 = (zans.chi2 / zans.dof )[*]

   
   return, zreturn
end
;------------------------------------------------------------------------------
