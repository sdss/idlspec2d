;+
; NAME:
;   zfind
;
; PURPOSE:
;   Find possible redshift matches for a set of spectra using a set of
;   eigen-templates.
;
; CALLING SEQUENCE:
;   result = zfind( objflux, objivar, hdr=hdr, fiberid=fiberid, $
;    eigenfile=, [ columns=, npoly=, $
;     zmin=, zmax=, pspace=, nfind=, width= ]
;
; INPUTS:
;   objflux    - Object fluxes [NPIXOBJ,NOBJ]
;   objivar    - Object inverse variances [NPIXOBJ,NOBJ]
;
; REQUIRED KEYWORDS:
;   hdr        - FITS header for objects, used for the following fields:
;                COEFF0, COEFF1, PLATEID, MJD.
;   fiberid    - Integer fiber ID number for each object.
;   eigenfile  - Input FITS file with an [NPIXSTAR,NSTAR] image with
;                either templates or eigenspectra.
;                The header keywords COEFF0, COEFF1 are used to specify
;                the wavelength mapping in log-10 Angstroms.
;
; OPTIONAL KEYWORDS:
;   columns    - Column numbers of the eigenspectra image to use in the
;                PCA fit; default to all columns.
;   npoly      - Number of polynomial terms to append to eigenspectra;
;                default to none.
;   zmin       - Minimum redshift to consider; default to no lower bound.
;   zmax       - Maximum redshift to consider; default to no upper bound.
;   pspace     - Keyword for ZCOMPUTE().
;   nfind      - Keyword for ZCOMPUTE().
;   width      - Keyword for ZCOMPUTE().
;
; OUTPUTS:
;   result     - Structure with redshift-fit information.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   zcompute()
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
    'plate'      ,  0L, $
    'mjd'        ,  0L, $
    'fiberid'    ,  0L, $
    'class'      ,  '', $
    'subclass'   ,  '', $
    'z'          , 0.0, $
    'z_err'      , 0.0, $
    'chi2'       , 0.0, $
    'dof'        ,  0L, $
    'wcoverage'  , 0.0, $
    'tfile'      ,  '', $
    'tcolumn'    ,  0L, $
    'npoly'      ,  0L, $
    'theta'      , fltarr(10), $
    'sigma2'     , 0.0, $
    'sigma2_err' , 0.0  $
   )

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
   wcoverage = total(objivar NE 0,1) * objdloglam

   if (n_elements(zmin) NE 0) then $
    pmin = floor( alog10(1.0 + zmin) / objdloglam )
   if (n_elements(zmax) NE 0) then $
    pmax = ceil( alog10(1.0 + zmax) / objdloglam )

   ;----------
   ; Read the template file, and optionally trim to only those columns
   ; specified by COLUMNS.
   ; Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.

   starflux = readfits(eigenfile, shdr)
   starloglam0 = sxpar(shdr, 'COEFF0')
   stardloglam0 = sxpar(shdr, 'COEFF1')
   npixstar = (size(starflux, /dimens))[0]

   if (keyword_set(columns)) then $
    starflux = starflux[*,columns]

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
   for iobj=0, nobj-1 do $
    result[*,iobj].fiberid = fiberid[iobj]
   result.z = zans.z
   result.z_err = zans.z_err
   result.chi2 = zans.chi2
   result.dof = zans.dof
   ntheta = n_elements(zans[0].theta)
   result.theta[0:ntheta-1] = zans.theta
   for iobj=0, nobj-1 do $
    result[*,iobj].wcoverage = wcoverage[iobj]
   result.tfile = tfile
   result.tcolumn = columns[0]
   result.npoly = npoly

   return, result
end
;------------------------------------------------------------------------------
