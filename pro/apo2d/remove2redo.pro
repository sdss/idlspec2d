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
;   logfile    - Fits file where apo results are stored 
;   outdir     - The spectrolog directory which should be checked 
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
     spawn, 'ls
     dirs = findfile() + 0L
     cd, strtrim(string(max(dirs)),2) 
   endif else cd, outdir 

   if NOT keyword_set(logfile) then begin
     logfile = findfile('logfile*fits')
     logfile = logfile[0]
   endif
   
   for i=1,4 do begin
     list = mrdfits(logfile, i)

     if keyword_set(plate) then begin
       good = where(list.plate EQ plate)
       if good[0] EQ - 1 then list = 0 $
       else list = list[good]
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

   if NOT keyword_set(expnum) then return

   expmin = min(expnum + 0L)
   expmax = max(expnum + 0L)
   nexp = expmax - expmin - 1
   if nexp LE 0 then return

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

   if missing EQ ' ' then return

   spawn, 'rm -f '+missing

   return
end


      
         
        
   

   

