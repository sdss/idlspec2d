;+
; NAME:
;   synthspec
;
; PURPOSE:
;   Construct synthetic spectrum from eigen-templates.
;
; CALLING SEQUENCE:
;   synflux = synthspec(zans, [ loglam=, hdr=, eigndir= ])
;
; INPUTS:
;   zans       - 
;
; OPTIONAL KEYWORDS:
;   loglam     - Log-10 wavelengths at which to synthesize the spectrum.
;   hdr        - If specified, then use this header to construct LOGLAM.
;                Either LOGLAM or HDR must be specified.
;   eigendir   - Directory for EIGENFILE; default to $IDLSPEC2D/templates.
;
; OUTPUTS:
;   synflux    - 
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   combine1fiber
;   concat_dir()
;   poly_array()
;   readfits()
;   sxpar()
;
; REVISION HISTORY:
;   20-Aug-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function synthspec, zans, loglam=objloglam, hdr=hdr, eigendir=eigendir

   if (n_elements(eigendir) EQ 0) then $
    eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
   tfile = strtrim(zans.tfile,2)

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   if (keyword_set(objloglam)) then begin
      naxis1 = n_elements(objloglam)
      objdloglam = objloglam[1] - objloglam[0]
   endif else if (keyword_set(hdr)) then begin
      naxis1 = sxpar(hdr, 'NAXIS1')
      objloglam0 = sxpar(hdr, 'COEFF0')
      objdloglam = sxpar(hdr, 'COEFF1')
      objloglam = objloglam0 + dindgen(naxis1) * objdloglam
   endif else begin
      print, 'Either LOGLAM or HDR must be specified'
      return, -1
   endelse

   ;----------
   ; If no template file, then return all zeros.

   if (tfile EQ '') then $
    return, fltarr(naxis1)

   ;----------
   ; Read the template file, and optionally trim to only those columns
   ; specified by COLUMNS.
   ; Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.

   starflux = readfits(djs_filepath(tfile, root_dir=eigendir), shdr)
   starloglam0 = sxpar(shdr, 'COEFF0')
   stardloglam0 = sxpar(shdr, 'COEFF1')

   ndim = size(starflux, /n_dimen)
   dims = size(starflux, /dimens)
   npixstar = dims[0]
   if (ndim EQ 1) then nstar = 1 $
    else nstar = dims[1]

   icol = where(zans.tcolumn NE -1, ncol)
   starflux = starflux[*,zans.tcolumn[icol]]

   ;----------
   ; Add more eigen-templates that represent polynomial terms.

   if (keyword_set(zans.npoly)) then $
    starflux = [ [starflux], [poly_array(npixstar,zans.npoly)] ]

   ;----------
   ; Construct the synthetic spectrum as a linear combination of templates

   synflux = starflux # zans.theta[0:ncol-1+zans.npoly]

   ;----------
   ; Apply redshift...

   starloglam = starloglam0 + dindgen(npixstar) * objdloglam
   combine1fiber, starloglam + alog10(1+zans.z), synflux, $
    newloglam=objloglam, newflux=newflux

   return, newflux
end
;------------------------------------------------------------------------------
