;+
; NAME:
;   zfind
;
; PURPOSE:
;   1-D reduction of spectra from 1 plate
;
; CALLING SEQUENCE:
;   zfind
;
; INPUTS:
;   platefile  - Plate file from spectro-2D
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Presently, just writes a file 'z-'+PLATEFILE.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   zcompute()
;   copy_struct
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;   sp1d_struct()
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function sp1d_struct

   result = create_struct( $
    name = 'ZANS', $
    'plate'   ,  0L, $
    'mjd'     ,  0L, $
    'fiberid' ,  0L, $
    'class'   ,  '', $
    'subclass',  '', $
    'z'       , 0.0, $
    'z_err'   , 0.0, $
    'chi2'    , 0.0, $
    'dof'     ,  0L, $
    'theta'   , fltarr(10) )

   return, result
end

;------------------------------------------------------------------------------
function zfind, objflux, objivar, hdr=hdr, fiberid=fiberid, $
 eigenfile=eigenfile, columns=columns, npoly=npoly, $
 zmin=zmin, zmax=zmax, pspace=pspace, nfind=nfind, width=width

   ndim = size(objflux, /n_dimen)
   if (ndim EQ 1) then nobj = 1 $
    else nobj = (size(objflux, /dimens))[1]

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')

   pmin = floor( alog10(1.0 + zmin) / objdloglam )
   pmax = ceil( alog10(1.0 + zmax) / objdloglam )

   ;----------
   ; Read the template file, and optionally trim to only those columns
   ; specified by COLUMNS.
   ; Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.

   starflux = readfits(eigenfile, shdr)
   starloglam0 = sxpar(shdr, 'COEFF0')
   stardloglam0 = sxpar(shdr, 'COEFF1')
   npixstar = sxpar(shdr, 'NAXIS1')

   if (keyword_set(columns)) then $
    starflux = starflux[*,columns]

   ;----------
   ; Add more eigentemplates that represent polynomial terms.

   if (keyword_set(npoly)) then $
    starflux = [ [starflux], [poly_array(npixstar,npoly)] ]

   ;----------
   ; Compute the redshift difference between the first pixel of the object
   ; spectra and the template.

   poffset = (objloglam0 - starloglam0) / objdloglam

   ;----------
   ; Compute the redshifts

   zans = zcompute(objflux, objivar, starflux, poffset=poffset, $
    pmin=pmin, pmax=pmax, pspace=pspace, nfind=nfind, width=width)

   ;----------
   ; Convert redshift (and error) from pixels to the conventional dimensionless
   ; value.

   indx = where(zans.dof GT 0)
   if (indx[0] NE -1) then begin
      zans[indx].z = 10.^(objdloglam * zans[indx].z) - 1.
      zans[indx].z_err = $
       alog(10d) * objdloglam * zans[indx].z_err * (1 + zans[indx].z)
   endif

   ;----------
   ; Copy into the output structure

   result = replicate(sp1d_struct(), nfind, nobj)
   result.plate = sxpar(hdr, 'PLATEID')
   result.mjd = sxpar(hdr, 'MJD')
   for ifind=0, nfind-1 do $
    result[ifind,*].fiberid = transpose(fiberid)
   result.z = zans.z
   result.z_err = zans.z_err
   result.chi2 = zans.chi2
   result.dof = zans.dof
   ntheta = n_elements(zans[0].theta)
   result.theta[0:ntheta-1] = zans.theta

   return, result
end
;------------------------------------------------------------------------------
