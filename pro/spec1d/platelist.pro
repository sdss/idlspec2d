;+
; NAME:
;   platelist
;
; PURPOSE:
;   Make list of reduced plates
;
; CALLING SEQUENCE:
;   platelist, [fullplatefile, outfile=]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   fullplatefile - Plate file(s) from spectro-2D; default to all files
;                   matching '*/spPlate*.fits'
;   outfile     - If set, then write output to this file.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Also look for the Spectro-1D files with names 'spZbest*.fits'.
;   For example, if PLATEFILE='spPlate-0306-51690.fits', then look for
;   the file ZBESTFILE = 'spZbest-0306-51690.fits'.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   djs_filepath()
;   fileandpath()
;   headfits()
;   mrdfits()
;   sxpar()
;
; REVISION HISTORY:
;   29-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platelist, fullplatefile, outfile=outfile

   if (NOT keyword_set(fullplatefile)) then $
    fullplatefile = findfile('*/spPlate*.fits', count=nfile) $
   else $
    nfile = n_elements(fullplatefile)
   fullplatefile = fullplatefile[ sort(fullplatefile) ]
   if (nfile EQ 0) then return

   ;----------
   ; Open output file

   if (keyword_set(outfile)) then $
    openw, olun, outfile, /get_lun $
   else $
    olun = -1L

   ;----------
   ; Determine names of associated files

   platefile = fileandpath(fullplatefile, path=path)
   platemjd = strmid(fileandpath(platefile), 8, 10)
   zbestfile = 'spZbest-' + platemjd + '.fits'

   ;----------
   ; Loop through all files

   printf, olun, 'PLATE  MJD   TILE SN2_G1 SN2_I1 SN2_G2 SN2_I2 ' $
    + 'Ngal Nqso Nsta Nunk Nsky'
   printf, olun, '-----  ----- ---- ------ ------ ------ ------ ' $
    + '---- ---- ---- ---- ----'

   for ifile=0, nfile-1 do begin
      fullzfile = djs_filepath(zbestfile[ifile], root_dir=path[ifile])
      hdr1 = headfits(fullplatefile[ifile])
      hdr2 = headfits(fullzfile)

      plate = sxpar(hdr1, 'PLATEID')
      tile = sxpar(hdr1, 'TILEID')
      mjd = sxpar(hdr1, 'MJD')
      snvec = [ sxpar(hdr1, 'SPEC1_G'), $
                sxpar(hdr1, 'SPEC1_I'), $
                sxpar(hdr1, 'SPEC2_G'), $
                sxpar(hdr1, 'SPEC2_I') ]
;      if (size(hdr2, /tname) EQ 'STRING') then qdone = 'DONE' $
;       else qdone = 'N/A'
      if (size(hdr2, /tname) EQ 'STRING') then begin
         zans = mrdfits(fullzfile, 1, /silent)
         class = strtrim(zans.class,2)
         nums = [ total(class EQ 'GALAXY'), $
                  total(class EQ 'QSO'), $
                  total(class EQ 'STAR'), $
                  total(class EQ 'UNKNOWN'), $
                  total(class EQ 'SKY') ]
      endif else begin
         nums = lonarr(5)
      endelse

      printf, olun, plate, mjd, tile, snvec, nums, $
       format='(i5,i7,i5,4f7.1,5i5)'

   endfor

   ;----------
   ; Close output file

   if (keyword_set(outfile)) then begin
      close, olun
      free_lun, olun
   endif

   return
end
;------------------------------------------------------------------------------
