;------------------------------------------------------------------------------
;+
; NAME:
;   plugmap_replace
;
; PURPOSE:
;   Replace plPlugMapM files in the speclog product with new files
;   if any of the FIBERID's differ.
;
; CALLING SEQUENCE:
;   plugmap_replace, [ mjd=, mjstart=, mjend= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd          - MJDs in the speclog product to inspect+replace
;   mjstart      - Start MJD in the speclog product to inspect+replace
;   mjend        - End MJD in the speclog product to inspect+replace
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Consider all files $SPECLOG_DIR/$MJD/plPlugMapM-*.par .
;   Replace if pluggings differ using files in $ASTROLOG_DIR/* .
;
; EXAMPLES:
;   Start with modified versions of the plPlugMapM files in a temp
;   directory (scan2), and copy into the speclog product:
;     setenv,'ASTROLOG_DIR=/home/schlegel/scan2'
;     setenv,'SPECLOG_DIR=/home/schlegel/products/NULL/speclog/trunk'
;     plugmap_replace
;   After the above, the changed files need to be "svn commit" from $SPECLOG_DIR
;
; BUGS:
;
; PROCEDURES CALLED:
;   yanny_free
;   yanny_par()
;   yanny_readone()
;   yanny_write
;
; REVISION HISTORY:
;   31-Mar-2011  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro plugmap_replace, mjd=mjd, mjstart=mjstart, mjend=mjend

   speclog_dir = getenv('SPECLOG_DIR')
   astrolog_dir = getenv('ASTROLOG_DIR')
   if (keyword_set(speclog_dir) * keyword_set(astrolog_dir) EQ 0) then $
    message, 'Must set SPECLOG_DIR and ASTROLOG_DIR'
   splog,file='replace.log'
   mjdlist = get_mjd_dir(speclog_dir, mjd=mjd, mjstart=mjstart, mjend=mjend)

   if (NOT keyword_set(mjdlist)) then begin
      splog, 'No MJDs found'
      return
   endif

   i = 0L
   for imjd=0L, n_elements(mjdlist)-1L do begin
      files = findfile(filepath('plPlugMapM-*.par', root_dir=speclog_dir, $
       subdir=mjdlist[imjd]), count=nfile)
      for ifile=0L, nfile-1L do begin
         file1 = (findfile(filepath(fileandpath(files[ifile]), $
          root_dir=astrolog_dir, subdir='*'), count=ct1))[0]
         if (ct1 GT 0) then begin
            print, i++, ' '+file1, files[ifile]
            dat1 = yanny_readone(file1)
            dat2 = yanny_readone(files[ifile])
            ndiff = long(total(dat1.fiberid NE dat2.fiberid))
            if (ndiff GT 0) then begin
               splog, 'Replacing '+files[ifile], ndiff
               spawn, '\cp -f '+file1+' '+files[ifile]
            endif
         endif
      endfor
   endfor
   splog, /close

   return
end
;------------------------------------------------------------------------------
