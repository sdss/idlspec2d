;+
; NAME:
;   multisynthspec
;
; PURPOSE:
;   Construct synthetic spectrum from eigen-templates.
;
; CALLING SEQUENCE:
;   synflux = multisynthspec(zans, [ loglam=, hdr=, eigendir= ])
;
; INPUTS:
;   zans       - Structure(s) with redshift-fit information (from SPREDUCE1D).
;
; OPTIONAL KEYWORDS:
;   loglam     - Log-10 wavelengths at which to synthesize the spectrum;
;                this can either be a single vector if all spectra have
;                the same wavelength mapping, or an array of [NPIX,NOBJ].
;   hdr        - If specified, then use this header to construct LOGLAM.
;                Either LOGLAM or HDR must be specified.
;   eigendir   - Directory for EIGENFILE; default to $IDLSPEC2D/templates.
;
; OUTPUTS:
;   synflux    - Synthetic spectra
;
; COMMENTS:
;   This routine is meant to be faster than SYNTHSPEC when generating
;   many spectra that share the same templates.  This works by generating
;   the B-spline coefficients only once per template, then evaluating
;   these B-splines at the exact wavelengths of each output spectrum.
;   This will NOT be exactly the same as the results from SYNTHSPEC,
;   because the individual eigen-spectra are spline-shifted before adding
;   them up, rather than adding them up then spline-shifting (as SYNTHSPEC
;   does).
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
function multisynthspec, zans, loglam=objloglam, hdr=hdr, eigendir=eigendir

   nobj = n_elements(zans)
   if (n_elements(eigendir) EQ 0) then $
    eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   if (keyword_set(objloglam)) then begin
      dims = size(objloglam, /dimens)
      naxis1 = dims[0]
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
   ; Construct the output spectra

   newflux = fltarr(naxis1, nobj)

   ;----------
   ; Get name(s) of template file, ignoring any that are blank strings

   tfiles = strtrim(zans.tfile,2)
   alltfile = tfiles[ uniq(tfiles, sort(tfiles)) ]
   if (N_elements(alltfile) EQ 1 AND alltfile[0] EQ '') then return, newflux
   alltfile = alltfile[where(alltfile NE '')]

   ;---------------------------------------------------------------------------
   ; Loop through each template file

   for ifile=0, n_elements(alltfile)-1 do begin

      ; List of objects that use this template file
      iobj = where(tfiles EQ alltfile[ifile], niobj)
      itheta = lonarr(niobj)

      ;----------
      ; Read the template file.
      ; (Assume that the wavelength binning is the same as for the objects
      ; in log-wavelength.)

      starflux = readfits(djs_filepath(alltfile[ifile], root_dir=eigendir), $
       shdr)
      starloglam0 = sxpar(shdr, 'COEFF0')
      stardloglam = sxpar(shdr, 'COEFF1')

      ndim = size(starflux, /n_dimen)
      dims = size(starflux, /dimens)
      npixstar = dims[0]
      starloglam = starloglam0 + dindgen(npixstar) * stardloglam
      if (ndim EQ 1) then nstar = 1 $
       else nstar = dims[1]

      ;----------
      ; Add as many polynomial terms as we might need for this template file

      maxnpoly = max([ zans[iobj].npoly ])
      if (keyword_set(maxnpoly)) then $
       starflux = [ [starflux], [poly_array(npixstar,maxnpoly)] ]

      ;----------
      ; Loop through each column in this template file

      for istar=0, nstar+maxnpoly-1 do begin

         ; Identify the coefficients for each object using this template/column
         if (istar LT nstar) then begin
            for i=0, niobj-1 do $
             itheta[i] = (where(zans[iobj[i]].tcolumn EQ istar))[0]
         endif else begin
            for i=0, niobj-1 do $
             if (istar-nstar LT zans[iobj[i]].npoly) then $
              itheta[i] = (where(zans[iobj[i]].tcolumn EQ -1))[0] + istar-nstar $
             else $
              itheta[i] = -1
         endelse

         indx = where(itheta NE -1, nthis)
         if (indx[0] NE -1) then begin
            thisobj = iobj[indx]
            thistheta = fltarr(nthis)
            for i=0, nthis-1 do $
             thistheta[i] = zans[thisobj[i]].theta[itheta[indx[i]]]
            if (size(objloglam, /n_dimen) EQ 1) then $
             thisloglam = objloglam # replicate(1,nthis) $
            else $
             thisloglam = objloglam[*,thisobj]
            for i=0, nthis-1 do $
             thisloglam[*,i] = thisloglam[*,i] - alog10(1+zans[thisobj[i]].z)

            combine1fiber, starloglam, starflux[*,istar], $
             newloglam=thisloglam, newflux=addflux
            addflux = reform(addflux, naxis1, nthis)

            for i=0, nthis-1 do $
             newflux[*,thisobj[i]] = newflux[*,thisobj[i]] $
              + addflux[*,i] * thistheta[i]
         endif
      endfor
   endfor

   return, newflux
end
;------------------------------------------------------------------------------
