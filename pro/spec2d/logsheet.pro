;+
; NAME:
;   logsheet
;
; PURPOSE:
;   Make a summary of header keywords in a directory of raw SDSS spectro files.
;   Only look at files matching "*.fit"
;
; CALLING SEQUENCE:
;   logsheet, [dir, outfile=]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   dir        - Directory name; default to '.'
;   outfile    - Output file
;
; OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   sdsshead()
;   sxpar()
;
; REVISION HISTORY:
;   06-Oct-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro logsheet, dir, outfile=outfile

   if (NOT keyword_set(dir)) then dir = '.'

   fullname = findfile(filepath('*.fit',root_dir=dir), count=nfile)

   if (nfile EQ 0) then begin
      print, 'No files found.'
   endif else begin

      ; Open output file
      if (keyword_set(outfile)) then begin
         get_lun, olun
         openw, olun, outfile
      endif else begin
         olun = -1
      endelse

      ; Remove the path from the file names
      shortname = strarr(nfile)
      for i=0, nfile-1 do begin
         res = str_sep(fullname[i],'/')
         shortname[i] = res[N_elements(res)-1]
      endfor

      ; Sort the files based upon exposure number + camera number
      isort = sort( strmid(shortname,7,8) + strmid(shortname,4,2) )
      fullname = fullname[isort]
      shortname = shortname[isort]

      printf, olun, $
       'DATEOBS    TAIHMS      FILENAME            PLATE CA EXPTIME  FLAVOR'
      printf, olun, $
       '---------- ----------- ------------------- ----- -- -------- ------'

      for i=0, nfile-1 do begin
         hdr = sdsshead(fullname[i])
         DATEOBS = sxpar(hdr, 'DATE-OBS')
         TAIHMS = sxpar(hdr, 'TAIHMS')
         PLATEID = sxpar(hdr, 'PLATEID')
         CAMERAS = sxpar(hdr, 'CAMERAS')
         EXPTIME = sxpar(hdr, 'EXPTIME')
         FLAVOR = sxpar(hdr, 'FLAVOR') + '               '

         ; Put a blank line before a new plate
         if (i EQ 0) then lastplate = PLATEID
         if (PLATEID NE lastplate) then begin
            printf, olun, ' '
            lastplate = PLATEID
         endif

         printf, olun, DATEOBS, TAIHMS, shortname[i], $
          PLATEID, CAMERAS, EXPTIME, FLAVOR, $
          format='(a10, " ", a11, " ", a19, i6, " ", a2, f9.2, " ", a15)'
      endfor

      ; Close output file
      if (keyword_set(outfile)) then begin
         close, olun
         free_lun, olun
      endif

   endelse

   return
end
