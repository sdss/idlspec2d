;+
; NAME:
;   zfind
;
; PURPOSE:
;   Find possible redshift matches for a set of spectra using a set of
;   eigen-templates.
;
; CALLING SEQUENCE:
;   result = zfind( objflux, objivar, hdr=hdr, $
;    eigenfile=, [eigendir=, columns=, npoly=, $
;     zmin=, zmax=, zguess=, pwidth=, pspace=, nfind=, width= ]
;
; INPUTS:
;   objflux    - Object fluxes [NPIXOBJ,NOBJ]
;   objivar    - Object inverse variances [NPIXOBJ,NOBJ]
;
; REQUIRED KEYWORDS:
;   hdr        - FITS header for objects, used to construct the wavelengths
;                from the following keywords: COEFF0, COEFF1.
;   eigenfile  - Input FITS file with an [NPIXSTAR,NSTAR] image with
;                either templates or eigenspectra.  If a wildcard appears
;                in the file name, then the file that appears last in a sort
;                is used.
;                The header keywords COEFF0, COEFF1 are used to specify
;                the wavelength mapping in log-10 Angstroms.
;
; OPTIONAL KEYWORDS:
;   eigendir   - Directory for EIGENFILE; default to $IDLSPEC2D/templates.
;   columns    - Column numbers of the eigenspectra image to use in the
;                PCA fit; default to all columns.
;   npoly      - Number of polynomial terms to append to eigenspectra;
;                default to none.
;   zmin       - Minimum redshift to consider; default to no lower bound.
;   zmax       - Maximum redshift to consider; default to no upper bound.
;   zguess     - Initial guess for redshift; search for a solution about
;                this value.  If specified with PWIDTH, then ZMIN and ZMAX
;                are ignoreed.
;   pwidth     - Search width in pixels about the intial guess redshift ZGUESS.
;                If specified with ZGUESS, then ZMIN and ZMAX are ignored.
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
;   One can specify a search domain for the redshift with ZMIN and ZMAX, or
;   with ZGUESS and PWIDTH.  If none of those parameters are set, then all
;   possible redshifts that overlap the object and star (template) are tested.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   concat_dir()
;   djs_filepath()
;   fileandpath()
;   readfits()
;   sxpar()
;   zcompute()
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
    'class'      ,  '', $
    'subclass'   ,  '', $
    'z'          , 0.0, $
    'z_err'      , 0.0, $
    'chi2'       , 0.0, $
    'dof'        ,  0L, $
    'rchi2diff'  , 0.0, $
    'tfile'      ,  '', $
    'tcolumn'    , lonarr(10) - 1L, $
    'npoly'      ,  0L, $
    'theta'      , fltarr(10), $
    'sigma2'     , 0.0, $
    'sigma2_err' , 0.0  $
   )

   return, result
end

;------------------------------------------------------------------------------
function zfind, objflux, objivar, hdr=hdr, $
 eigenfile=eigenfile, eigendir=eigendir, columns=columns, npoly=npoly, $
 zmin=zmin, zmax=zmax, zguess=zguess, pwidth=pwidth, $
 pspace=pspace, nfind=nfind, width=width

   if (n_elements(eigendir) EQ 0) then $
    eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')

   ndim = size(objflux, /n_dimen)
   if (ndim EQ 1) then nobj = 1 $
    else nobj = (size(objflux, /dimens))[1]

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')

   if (n_elements(zmin) NE 0) then $
    pmin = floor( alog10(1.0 + zmin) / objdloglam )
   if (n_elements(zmax) NE 0) then $
    pmax = ceil( alog10(1.0 + zmax) / objdloglam )

   if (n_elements(zguess) GT 0 AND keyword_set(pwidth)) then begin
      pmin = floor( alog10(1.0 + zguess) / objdloglam - 0.5 * pwidth )
      pmax = floor( alog10(1.0 + zguess) / objdloglam + 0.5 * pwidth )
   endif

   ;----------
   ; Find the most recent template file matching EIGENFILE

   allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
   if (ct EQ 0) then $
    message, 'Unable to find EIGENFILE matching '+eigenfile
   thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]
   splog, 'Selecting EIGENFILE=', thisfile

   ;----------
   ; Read the template file, and optionally trim to only those columns
   ; specified by COLUMNS.
   ; Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.

   starflux = readfits(thisfile, shdr)
   starloglam0 = sxpar(shdr, 'COEFF0')
   stardloglam0 = sxpar(shdr, 'COEFF1')

   ndim = size(starflux, /n_dimen)
   dims = size(starflux, /dimens)
   npixstar = dims[0]
   if (ndim EQ 1) then nstar = 1 $
    else nstar = dims[1]

   if (n_elements(columns) NE 0) then begin
      starflux = starflux[*,columns]
   endif else begin
      columns = lindgen(nstar)
   endelse

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
   result.z = zans.z
   result.z_err = zans.z_err
   result.chi2 = zans.chi2
   result.dof = zans.dof
   ntheta = n_elements(zans[0].theta)
   result.theta[0:ntheta-1] = zans.theta
   result.tfile = fileandpath(thisfile)
   for icol=0, n_elements(columns)-1 do $
    result.tcolumn[icol] = columns[icol]
   result.npoly = npoly

   return, result
end
;------------------------------------------------------------------------------
