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
;    [ starflux=, starloglam0=, stardloglam=, $
;    eigenfile=, eigendir=, columns=, npoly=, $
;    zmin=, zmax=, zguess=, pwidth=, nfind=, width=, $
;    zans_fixed=, _EXTRA= ]
;
; INPUTS:
;   objflux    - Object fluxes [NPIXOBJ,NOBJ]
;   objivar    - Object inverse variances [NPIXOBJ,NOBJ]
;
; REQUIRED KEYWORDS:
;   hdr        - FITS header for objects, used to construct the wavelengths
;                from the following keywords: COEFF0, COEFF1.
;
; OPTIONAL KEYWORDS:
;   starflux   - Eigenspectra [NPIXSTAR,NSTAR].
;   starloglam0- Zero-point of log-10(Angstrom) wavelength mapping of STARFLUX.
;   stardloglam- Wavelength spacing for STARFLUX in log-10(Angstroms)
;   eigenfile  - Input FITS file with an [NPIXSTAR,NSTAR] image with
;                either templates or eigenspectra.  If a wildcard appears
;                in the file name, then the file that appears last in a sort
;                is used.
;                The header keywords COEFF0, COEFF1 are used to specify
;                the wavelength mapping in log-10 Angstroms.
;                This must be set if STARFLUX,STARLOGLAM0 are not set.
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
;   nfind      - Keyword for ZCOMPUTE().
;   width      - Keyword for ZCOMPUTE().
;   _EXTRA     - Keywords for ZCOMPUTE(), such as PSPACE, DOPLOT, DEBUG, VERBOSE
;
; OUTPUTS:
;   result     - Structure with redshift-fit information.  Structure
;                elements are left blank if fewer than NFIND peaks are found.
;
; OPTIONAL OUTPUTS:
;   zans_fixed - Structure with fit information for those templates
;                fixed in redshift if FIXED_TEMPLATE is passed.
;
; COMMENTS:
;   One can specify a search domain for the redshift with ZMIN and ZMAX, or
;   with ZGUESS and PWIDTH.  If none of those parameters are set, then all
;   possible redshifts that overlap the object and star (template) are tested.
;
;   Mask any pixels on the templates where the first template contains zeros.
;   This is useful, in particular, where the stellar templates have zeros
;   at the beginning or end of the spectral range due lack of wavelength
;   coverage.
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
;   splog
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
    'rchi2'      , 0.0, $
    'dof'        ,  0L, $
    'rchi2diff'  , 0.0, $
    'tfile'      ,  '', $
    'tcolumn'    , lonarr(10) - 1L, $
    'npoly'      ,  0L, $
    'theta'      , fltarr(10), $
    'theta_covar', fltarr(10,10), $
    'vdisp'      , 0.0, $
    'vdisp_err'  , 0.0, $
    'vdispz'     , 0.0, $
    'vdispz_err' , 0.0, $
    'vdispchi2'  , 0.0, $
    'vdispnpix'  , 0.0, $
    'vdispdof'   ,  0L  $
   )

   return, result
end

;------------------------------------------------------------------------------
function zfind, objflux, objivar, hdr=hdr, $
 starflux=starflux, starloglam0=starloglam0, stardloglam=stardloglam, $
 eigenfile=eigenfile, eigendir=eigendir, columns=columns, npoly=npoly, $
 zmin=zmin, zmax=zmax, zguess=zguess, pwidth=pwidth, $
 nfind=nfind, width=width, zans_fixed=zans_fixed, _EXTRA=EXTRA

   if (n_elements(eigendir) EQ 0) then $
    eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')

   ndim = size(objflux, /n_dimen)
   if (ndim EQ 1) then nobj = 1 $
    else nobj = (size(objflux, /dimens))[1]
   npixobj = (size(objflux))[1]

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
      if (keyword_set(width)) then width1 = width $
       else width1 = pwidth
      pmin = floor( alog10(1.0 + zguess) / objdloglam - 0.5 * (pwidth+1+width1) )
      pmax = floor( alog10(1.0 + zguess) / objdloglam + 0.5 * (pwidth+1+width1) )
   endif

   if (keyword_set(eigenfile)) then begin
      ;----------
      ; Find the most recent template file matching EIGENFILE

      allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
      if (ct EQ 0) then $
       message, 'Unable to find EIGENFILE matching '+eigenfile
      thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]
      splog, 'Selecting EIGENFILE=' + thisfile
      if (keyword_set(columns)) then $
       splog, 'Selecting columns=', columns

      ;----------
      ; Read the template file, and optionally trim to only those columns
      ; specified by COLUMNS.
      ; Assume that the wavelength binning is the same as for the objects
      ; in log-wavelength.

      starflux = readfits(thisfile, shdr,/silent)
      if n_elements(starflux) LE 1 then begin
        splog, 'Looking for a bspline structure ', thisfile
        bspline_set = mrdfits(thisfile, 1, shdr, /silent) 
      endif else begin
        starloglam0 = sxpar(shdr, 'COEFF0')
        stardloglam0 = sxpar(shdr, 'COEFF1')
      endelse
   endif

   if (NOT keyword_set(bspline_set) AND NOT keyword_set(starflux) AND $
                                        NOT keyword_set(starloglam0)) then begin
      message, 'Either EIGENFILE or STARFLUX,STARLOGLAM0 must be set'
   endif

   if (keyword_set(stardloglam)) then begin
      message, 'Object wavelength spacing and STARDLOGLAM must be the same'
   endif

   ndim = size(starflux, /n_dimen)
   dims = size(starflux, /dimens)
   npixstar = dims[0]
   if (ndim EQ 0) then nstar = 0 $
   else if (ndim EQ 1) then nstar = 1 $
   else nstar = dims[1]

   ; Trim to specified columns
   if (keyword_set(nstar)) then begin 
      if (n_elements(columns) NE 0) then starflux = starflux[*,columns] $
       else columns = lindgen(nstar)
   endif

   ;----------
   ; Add more eigen-templates that represent polynomial terms.

   if (keyword_set(npoly) AND keyword_set(starflux)) then begin
       starflux = [ [starflux], [poly_array(npixstar,npoly)] ]
   endif else if (keyword_set(npoly)) then begin
      starflux = poly_array(npixobj,npoly)
   endif

   ;----------
   ; Compute the redshift difference between the first pixel of the object
   ; spectra and the template.

   poffset = 0
   if keyword_set(starloglam0) then poffset = (objloglam0 - starloglam0) / objdloglam

   ;----------
   ; Compute the redshifts

   if keyword_set(bspline_set) then begin
      ; call the QSO composite fit with multiplicative polynomials

      starset = bspline_set
      starset.xmin = (starset.xmin - objloglam0) / objdloglam
      starset.xmax = (starset.xmax - objloglam0) / objdloglam
      starset.fullbkpt = (starset.fullbkpt - objloglam0) / objdloglam

      zans = zcompute_qso(objflux, objivar, starset, starflux, poffset=poffset, $
       pmin=pmin, pmax=pmax, nfind=nfind, width=width, $
       plottitle=plottitle, zans_fixed=zans_fixed, _EXTRA=EXTRA)
   endif else begin
      ; Mask any pixels on the templates where the first template contains zeros
      starmask = starflux[*,0] NE 0
      zans = zcompute(objflux, objivar, starflux, starmask, poffset=poffset, $
       pmin=pmin, pmax=pmax, nfind=nfind, width=width, $
       plottitle=plottitle, zans_fixed=zans_fixed, _EXTRA=EXTRA)
   endelse

   ;----------
   ; Convert redshift (and error) from pixels to the conventional dimensionless
   ; value.  Do not modify any errors that are less than zero, since those
   ; can be used as just warning flags from the fit.

; ASB: change to disallow negative-model stellar fits before the ZWARNING stage:
;   indx = where(zans.dof GT 0, npeak)
   indx = where((zans.dof GT 0) and ((n_elements(columns) NE 1) or (zans.theta[0] GT 0.)), npeak)
   if (npeak GT 0) then $
    zans[indx].z = 10.^(objdloglam * zans[indx].z) - 1.

   jndx = where(zans.dof GT 0 and zans.z_err GE 0)
   if (jndx[0] NE -1) then $
    zans[jndx].z_err = $
     alog(10d) * objdloglam * zans[jndx].z_err * (1 + zans[jndx].z)

   ;----------
   ; Copy valid peaks into the output structure

   result = replicate(sp1d_struct(), nfind, nobj)
   if (npeak GT 0) then begin
      result[indx].z = zans[indx].z
      result[indx].z_err = zans[indx].z_err
      result[indx].rchi2 = zans[indx].chi2 / (zans[indx].dof > 1)
      result[indx].dof = zans[indx].dof
      ntheta = n_elements(zans[0].theta)
      result[indx].theta[0:ntheta-1] = zans[indx].theta
      result[indx].theta_covar[0:ntheta-1,0:ntheta-1] = zans[indx].theta_covar
      if (keyword_set(eigenfile)) then $
       result[indx].tfile = fileandpath(thisfile)
      for icol=0, n_elements(columns)-1 do $
       result[indx].tcolumn[icol] = columns[icol]
      if (keyword_set(npoly)) then result.npoly = npoly
   endif

   return, result
end
;------------------------------------------------------------------------------
