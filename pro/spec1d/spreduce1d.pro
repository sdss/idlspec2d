;+
; NAME:
;   spreduce1d
;
; PURPOSE:
;   1-D reduction of spectra from 1 plate
;
; CALLING SEQUENCE:
;   spreduce1d, platefile, templatefile, outfile
;
; INPUTS:
;   platefile  - Plate file from spectro-2D
;   templatefile - Template file(s)
;   outfile    - Ouptut file
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
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
;   spec_append()
;   struct_addtags()
;   sxpar()
;   veldisp
;
; INTERNAL SUPPORT ROUTINES:
;   sp1d_struct()
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function sp1d_struct, nobj

   result = { $
    plate               : 0L, $
    mjd                 : 0L, $
    fiberid             : 0L, $
    primtarget          : 0L, $
    sectarget           : 0L $
   }

   return, replicate(result, nobj)
end

;------------------------------------------------------------------------------
pro spreduce1d, platefile, templatefile, outfile

   ;----------
   ; Read the 2D output file

   objflux = mrdfits(platefile,0,hdr)
   objivar = mrdfits(platefile,1)
   andmask = mrdfits(platefile,2)
;   ormask = mrdfits(platefile,3)
   plugmap = mrdfits(platefile,4)

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   npix = sxpar(hdr, 'NAXIS1')
   nobj = sxpar(hdr, 'NAXIS2')
   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')
   objwave = 10^(objloglam0 + objdloglam * findgen(npix))

   ;----------
   ; Read the template files
   ; Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.

   ntemplate = n_elements(templatefile)
   for it=0, ntemplate-1 do begin
      tfile = filepath(templatefile[it], $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
      djs_readcol, twave, tflux, tivar
      tloglam0 = alog10(twave[0])

      if (it EQ 0) then begin
         templateloglam0 = tloglam0
      endif else begin
         templateloglam0 = [templateloglam0, tloglam0]
      endelse

      templateflux = spec_append(templateflux, tflux)
      templateivar = spec_append(templateivar, tivar)
   endfor

   ;----------
   ; If the number of spectral pixels is larger in either the objects
   ; or the templates, than increase the smaller's size to agree.

   nptemp = (size(templateflux))[0]
   if (nptemp LT npix) then begin
      templateflux = [templateflux, fltarr(npix-nptemp,ntemplate)]
      templateivar = [templateivar, fltarr(npix-nptemp,ntemplate)]
   endif else if (nptemp GT npix) then begin
      objflux = [objflux, fltarr(nptemp-npix,nobj)]
      objivar = [objivar, fltarr(nptemp-npix,nobj)]
   endif

   ;----------
   ; Mask out the object spectra near the 5577 sky line

   linemask = objwave GE 5573 AND objwave LE 5586 ; =1 for bad pixels
   linemask = rebin(linemask, npix, nobj)

   objflux = djs_maskinterp(objflux, linemask, iaxis=0)
   objivar = objivar * (1 - linemask)

   ;----------
   ; Compute the redshift difference between the first pixel of the object
   ; spectra and each template

   zoffset = templateloglam0 - objloglam0

   ;----------
   ; Compute the redshifts and velocity dispersions

   vres = veldisp(objflux, objivar, templateflux, templateivar, $
    zoffset=zoffset)

   result = sp1d_struct(nobj)
   result.plate = sxpar(hdr, 'PLATEID')
   result.mjd = sxpar(hdr, 'MJD')
   result.fiberid = plugmap.fiberid
   result.primtarget = plugmap.primtarget
   result.sectarget = plugmap.sectarget
   result = struct_addtags(result, vres)

   ;----------
   ; Write the output file

   mwrfits, result, outfile, /create

   return
end
;------------------------------------------------------------------------------
