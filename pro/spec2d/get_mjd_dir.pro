;+
; NAME:
;   get_mjd_dir
;
; PURPOSE:
;   Get directory list of matching MJD directories.
;
; CALLING SEQUENCE:
;   mjdlist = get_mjd_dir( [topdir, mjd=, mjstart=, mjend= ])
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   topdir     - Search for MJD directories under this top level directory.
;   mjd        - Look for raw data files in TOPINDIR/MJD; default to all
;                subdirectories.  Note that this need not be integer-valued,
;                but could be for example '51441_test'.  Wildcards are allowed,
;                e.g. '514*'.
;   mjstart    - Starting MJD.
;   mjend      - Ending MJD.
;
; OUTPUT:
;   mjdlist    - List of matching direcotries (string array).
;
; COMMENTS:
;   Do not include in the output list any empty directories.
;
; EXAMPLES:
;   If the path '/data/rawdata' contains the subdirectories '0123', '0124',
;   '0200', 'test0123' and the file '0125', then:
;   IDL> print, get_mjd_dir('/data/rawdata')
;        0123 0124 0200 test
;   IDL> print, get_mjd_dir('/data/rawdata', mjstart=124)
;        0124 0200
;   IDL> print, get_mjd_dir('/data/rawdata', mjd='012*')
;        0123 0124
;
; BUGS:
;   Note we assume Unix directory names and we spawn the Unix 'ls' command.
;
; PROCEDURES CALLED:
;   fileandpath()
;
; REVISION HISTORY:
;   30-May-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------

function get_mjd_dir, topdir, mjd=mjd, mjstart=mjstart, mjend=mjend

   ;----------
   ; Save the current directory name, and change directories to TOPDIR

   if (keyword_set(topdir)) then begin
      cd, topdir, current=olddir
   endif else begin
      cd, current=olddir
      cd, olddir
      topdir = olddir
   endelse

   ;----------
   ; Find all directories, or all directories matching MJD

   if (keyword_set(mjd)) then begin
      tmpstring = ''
      if (NOT keyword_set(topdir)) then topdir = '.'
      for i=0, N_elements(mjd)-1 do $
       tmpstring = tmpstring + ' ' $
        + filepath('*'+strtrim(string(mjd[i]),2), root_dir=topdir)
      spawn, '\ls -d '+tmpstring, mjdlist
   endif else begin
      mjdlist = findfile()
   endelse
   if (NOT keyword_set(mjdlist)) then return, ''

   ;----------
   ; Strip leading directory names from MJDLIST.

   mjdonly = strarr(n_elements(mjdlist))
   for imjd=0, N_elements(mjdlist)-1 do begin
      mjdonly[imjd] = fileandpath(mjdlist[imjd], path=newtopdir)

      ; Remove any directories not >= MJDSTART or <= MJDEND, if those
      ; keywords are set.
      if (keyword_set(mjstart) OR keyword_set(mjend)) then begin
         cc = strmid(mjdonly[imjd],0,1) ; first character of this dir name
         if (cc GE '0' AND cc LE '9') then begin
            fixmjd = long(mjdonly[imjd])
            if (keyword_set(mjstart)) then $
             if (fixmjd LT long(mjstart)) then $
              mjdonly[imjd] = ''
            if (keyword_set(mjend)) then $
             if (fixmjd GT long(mjend)) then $
              mjdonly[imjd] = ''
         endif else begin
            ; Reject this directory name because it does not start w/ a number
            mjdonly[imjd] = ''
         endelse
      endif

      ; Remove any directories that do not contain any files (or that are
      ; not directories at all).
      if (keyword_set(mjdonly[imjd])) then begin
         junk = findfile(filepath('', root_dir=mjdlist[imjd]), count=ct)
         if (ct EQ 0) then mjdonly[imjd] = ''
      endif
   endfor

   ii = where(mjdonly NE '', ct)
   if (ct EQ 0) then begin
      return, ''
   endif else begin
      mjdlist = mjdlist[ii]
      mjdonly = mjdonly[ii]
   endelse

   cd, olddir
   if (keyword_set(newtopdir)) then topdir = newtopdir
   return, mjdonly
end
;------------------------------------------------------------------------------
