;+
; NAME:
;   platemerge
;
; PURPOSE:
;   Merge all Spectro-1D outputs with tsObj files.
;
; CALLING SEQUENCE:
;   platemerge, [fullplatefile, outfile=]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   fullzfile   - Redshift file(s) from spectro-1D; default to all files
;                   matching '*/spZbest*.fits'
;   outfile     - Output file name; default to 'spAll.fits'
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
;
; PROCEDURES CALLED:
;   copy_struct_inx
;   headfits()
;   mrdfits()
;   mwrfits
;   sxpar()
;
; REVISION HISTORY:
;   30-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platemerge, fullzfile, outfile=outfile

   if (NOT keyword_set(fullzfile)) then $
    fullzfile = findfile('*/spZbest*.fits', count=nfile) $
   else $
    nfile = n_elements(fullzfile)
   fullzfile = fullzfile[ sort(fullzfile) ]
   if (nfile EQ 0) then return

   if (NOT keyword_set(outfile)) then outfile = 'spAll.fits'

   nout = nfile * 640
   print, 'Total number of objects = ', nout

   ;----------
   ; Loop through each file

   for ifile=0, nfile-1 do begin
      print,'File ',ifile, ' of ', nfile,': '+fullzfile[ifile]

      hdr = headfits(fullzfile[ifile])
      plate = sxpar(hdr, 'PLATEID')
      zans = mrdfits(fullzfile[ifile], 1)
      tsobj = plug2tsobj(plate, zans.plug_ra, zans.plug_dec)

      if (ifile EQ 0) then begin
         outdat = replicate( create_struct(zans[0], tsobj[0]), nout)
         struct_assign, {junk:0}, outdat ; Zero-out all elements
      endif

      copy_struct_inx, zans, outdat, index_to=lindgen(640)+640*ifile
      copy_struct_inx, tsobj, outdat, index_to=lindgen(640)+640*ifile

   endfor

   ;----------
   ; Write the output file

   mwrfits, outdat, outfile

   return
end
;------------------------------------------------------------------------------
