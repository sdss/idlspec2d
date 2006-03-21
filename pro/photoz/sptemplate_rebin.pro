;+
; NAME:
;   sptemplate_rebin
;
; PURPOSE:
;   Read and optionally rebin spectral templates
;
; CALLING SEQUENCE:
;   starflux = sptemplate_rebin([ loglam, dloglam=, eigenfile=, eigendir=, $
;    columns=, npoly=, starloglam= ])
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   loglam     - Log10 of wavelengths [vacuum Ang]; need not be uniformly
;                spaced; if not specified, then templates are returned in
;                their native wavelength scale as defined by OUTLOGLAM
;   dloglam    - If one specifies DLOGLAM instead of LOGLAM, then the
;                output wavelengths will span those in the template file
;                but spaced by this amount in log-wavelength
;   eigenfile  - Input FITS file with an [NPIXSTAR,NSTAR] image with
;                either templates or eigenspectra.  If a wildcard appears
;                in the file name, then the file that appears last in a sort
;                is used.
;                The header keywords COEFF0, COEFF1 are used to specify
;                the wavelength mapping in log-10 Angstroms.
;                Default to 'spEigenGal-*.fits'.
;   eigendir   - Directory for EIGENFILE; default to $IDLSPEC2D/templates.
;   columns    - Column numbers of the eigenspectra image to use in the
;                PCA fit; default to all columns.
;   npoly      - Number of polynomial terms to append to eigenspectra;
;                default to none.
;
; OUTPUTS:
;   starflux   - Output template spectra, typically in units of
;                erg/s/cm^2/Ang regardless of the wavelength spacing
;
; OPTIONAL OUTPUTS:
;   starloglam - Log10 of output wavelengths [vacuum Ang]; if LOGLAM is
;                specified, then this is the same
;
; COMMENTS:
;   There are two methods for rebinning the spectra.  If the output
;   wavelength scale is coarser than the input, then the flux is dumped
;   into pixels using CIC assignment.  Otherwise, BSPLINE_ITERFIT is used
;   to resample the spectra.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/templates/*.fits
;
; PROCEDURES CALLED:
;   combine1fiber
;   djs_filepath()
;   readfits()
;   poly_array()
;   populate_image
;   sxpar()
;
; REVISION HISTORY:
;   21-Mar-2006  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
function sptemplate_rebin, loglam1, dloglam=dloglam1, $
 eigenfile=eigenfile1, eigendir=eigendir1, $
 columns=columns, npoly=npoly, starloglam=starloglam

   if (keyword_set(eigenfile1)) then eigenfile = eigenfile1 $
    else eigenfile = 'spEigenGal-*.fits'
   if (size(eigendir1,/tname) EQ 'STRING') then eigendir = eigendir1 $
    else eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
   if (keyword_set(loglam1) AND keyword_set(dloglam1)) then $
    message, 'Cannot specify both LOGLAM and DLOGLAM!'

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

   starflux = readfits(thisfile, shdr, /silent)
   starloglam0 = sxpar(shdr, 'COEFF0')
   stardloglam = sxpar(shdr, 'COEFF1')

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
      nstar = nstar + npoly
   endif else if (keyword_set(npoly)) then begin
      starflux = poly_array(npixobj,npoly)
      nstar = nstar + npoly
   endif

   if (arg_present(starloglam) OR keyword_set(loglam1) $
    OR keyword_set(dloglam1)) then $
    starloglam = starloglam0 + stardloglam * dindgen(npixstar)

   ;----------
   ; If LOGLAM is specified, then remap these templates onto
   ; the new wavelength scale

   if (keyword_set(loglam1) EQ 0 AND keyword_set(dloglam1) EQ 0) $
    then return, starflux

   if (keyword_set(loglam1)) then begin
      loglam = loglam1
      npix = n_elements(loglam)
      if (npix LT 2) then $
       message, 'LOGLAM must have at least two elements!'
      dloglam = abs(loglam[1] - loglam[0])
   endif else begin
      dloglam = dloglam1
      minlog = min(starloglam, max=maxlog)
      npix = floor((maxlog - minlog) / dloglam)
      loglam = minlog + dloglam * dindgen(npix)
   endelse

   newflux = fltarr(npix,nstar)

   if (dloglam GT stardloglam) then begin
      thispix = interpol(dindgen(npix), loglam, starloglam)
      for istar=0, nstar-1 do begin
         newmask1 = fltarr(npix)
         newflux1 = fltarr(npix)
         populate_image, newmask1, thispix, weights=(starflux[*,istar] NE 0), $
          assign='cic'
         populate_image, newflux1, thispix, weights=starflux[*,istar], $
          assign='cic'
         newflux[*,istar] = newflux1 * (newmask1 GT 0) $
          / (newmask1 + (newmask1 LE 0))
      endfor
   endif else begin
      for istar=0, nstar-1 do begin
         combine1fiber, starloglam, starflux[*,istar], $
          newloglam=loglam, newflux=newflux1
         newflux[*,istar] = newflux1
      endfor
   endelse
   starloglam = loglam

   return, newflux
end
;------------------------------------------------------------------------------
