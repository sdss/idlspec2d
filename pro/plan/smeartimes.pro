;+
; NAME:
;   smeartimes
;
; PURPOSE:
;   List the TAI timestamps of 'smear' exposures in a log file.
;
; CALLING SEQUENCE:
;   smeartimes
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   All raw spectro files (sdR files) have their headers read and
;   their start and end time stamps computed.  The purpose of this
;   is to get an exact list of times and determine from the TPM logs
;   whether the telescope moved in the way that it was supposed to
;   during these 'smear' exposures.
;
;   The file names and their timestamps are written to the following files.
;     timestamps-all.log
;     timestamps-smear.log
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   sdssproc
;   splog
;   struct_print
;
; REVISION HISTORY:
;   10-Dec-2002  Written by David Schlegel, Princeton (not checked in then)
;-
;------------------------------------------------------------------------------
pro smeartimes

   ;----------
   ; Find all the b1 spFrame files, and sort by exposure number

   spawn, 'ls -d ' + getenv('BOSS_SPECTRO_DATA') + '/?????', mjdlist
   for imjd=0L, n_elements(mjdlist)-1 do begin
      splog, 'Searching for files in ' + mjdlist[imjd]
      flist = findfile( $
       filepath('sdR-b1-????????.fit*', root_dir=mjdlist[imjd]), count=ct)
      if (ct GT 0) then begin
         if (NOT keyword_set(files)) then files = flist $
          else files = [files, flist]
      endif
   endfor

   files = files[ sort(fileandpath(files)) ]
   nfile = n_elements(files)

   ;----------
   ; Loop through, identifying smear frames.

   hstruct1 = create_struct( $
    'mjd'      , 0L, $
    'plate'    , 0L, $
    'expnum'   , 0L, $
    'flavor'   , '', $
    'quality'  , '', $
    'tai_beg'  , 0d, $
    'tai_end'  , 0d, $
    'exptime'  , 0d, $
    'obscomm' , ''  )
   hstruct = replicate(hstruct1, nfile)

   for ifile=0L, nfile-1 do begin
      splog, 'Reading file ', ifile, ' of ', nfile, $
       ': ', fileandpath(files[ifile])

      sdssproc, files[ifile], hdr=hdr
      hstruct[ifile].mjd = sxpar(hdr,'MJD')
      hstruct[ifile].plate = sxpar(hdr,'PLATEID')
;      hstruct[ifile].expnum = sxpar(hdr,'EXPOSURE')
      hstruct[ifile].expnum = long( strmid(fileandpath(files[ifile]),7,8) )
      hstruct[ifile].flavor = strtrim(sxpar(hdr,'FLAVOR'),2)
      hstruct[ifile].quality = strtrim(sxpar(hdr,'QUALITY'),2)
      hstruct[ifile].tai_beg = sxpar(hdr,'TAI-BEG')
      hstruct[ifile].tai_end = sxpar(hdr,'TAI-END')
      hstruct[ifile].exptime = sxpar(hdr,'EXPTIME')
      hstruct[ifile].obscomm = strtrim(sxpar(hdr,'OBSCOMM'),2)
   endfor

   ; Select the smear exposures, and exclude the star cluster plates
   ismear = where(hstruct.expnum NE 0 $
    AND (hstruct.plate LT 798 OR hstruct.plate GT 801) $
    AND hstruct.obscomm NE '{BD17 Smear Pattern}' $
    AND hstruct.flavor EQ 'smear', nsmear)

   formatcodes = [ ['TAI_BEG', 'f11.0'], $
                   ['TAI_END', 'f11.0'] ]

   struct_print, hstruct, filename='timestamps-all.log', $
    formatcodes=formatcodes
   struct_print, hstruct[ismear], filename='timestamps-smear.log', $
    formatcodes=formatcodes

   return
end
;------------------------------------------------------------------------------
