;+
; NAME:
;   remove2redo
;
; PURPOSE:
;   Removes raw data files which are missing in the current logfile-$MJD.fits
;
; CALLING SEQUENCE:
;   remove2redo, [ logfile=logfile, outdir=outdir, plate=]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   logfile    - FITS log file to which Son-of-Spectro results are ouput;
;                default to the file 'logfile*fits' in the OUTDIR directory.
;   outdir     - The spectrolog directory which should be checked; default
;                to the largest MJD in /data/spectro/spectrologs/$MJD.
;   plate      - Specify a single plate if only those entries should be checked
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   03-April-2001  Written by Scott Burles, FNAL
;-
;------------------------------------------------------------------------------
pro remove2redo, logfile=logfile, outdir=outdir, plate=plate

   if NOT keyword_set(outdir) then begin
     cd, '/data/spectro/spectrologs'
     dirs = findfile() 
     dirs = dirs[where(dirs GT '50000' AND dirs LT '90000')]
     cd, strtrim(string(max(dirs)),2) 
   endif else cd, outdir 

   if NOT keyword_set(logfile) then begin
     logfile = findfile('logfile*fits')
     logfile = logfile[0]
   endif

   print, 'Working on ', logfile   
   for i=1,4 do begin
     list = mrdfits(logfile, i, /silent)


     if size(list, /tname) NE 'INT' then begin
       if keyword_set(plate) then begin
         good = where(list.plate EQ plate)
         if good[0] EQ -1 then list = 0 $
         else list = list[good]
       endif
     endif

     if size(list, /tname) NE 'INT' then begin
       if NOT keyword_set(expnum) then begin
         expnum = list.expnum
         camera = list.camera
       endif else begin
         expnum = [expnum, list.expnum]
         camera = [camera, list.camera]
       endelse
     endif
   endfor

   if NOT keyword_set(expnum) then begin
     print, ''
     print, 'Did not find any exposures to redo!'
     print, ''
     return
   endif

   expmin = min(expnum + 0L)
   expmax = max(expnum + 0L)
   nexp = expmax - expmin 
   if nexp LE 0 then begin
     print, ''
     print, 'Did not find any exposures to redo!'
     print, ''
     return
   endif

   cam = ['b1', 'b2', 'r1', 'r2']
   missing = ' '


   for i=expmin, expmax - 1 do begin 
     for icam = 0,3 do begin
       gotcha = where(cam[icam] EQ camera AND i EQ expnum + 0L)
       if gotcha[0] EQ -1 then $
         missing = missing + string(format='(a,a,a,i8.8,a)', $
                      ' sdR-',cam[icam],'-',i,'.fit')
     endfor
   endfor

   print, missing
   if missing EQ ' ' then return

   mjd = strmid(logfile, 8, 5)
   cd, '/data/spectro/'+mjd

   print, 'rm -f '+missing
   spawn, 'rm -f '+missing

   print, ''
   print, 'Successfully removed files, sos should get these on the next pass!'
   print, ''
   
   return
end
;------------------------------------------------------------------------------
