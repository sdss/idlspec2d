;+
; NAME:
;   platelinks
;
; PURPOSE:
;   Build links from another directory to the spectroscopic output files.
;
; CALLING SEQUENCE:
;   platelinks, outdir, [ public= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   outdir      - Directory in which to build all links.  These links are
;                 built to the files in $BOSS_SPECTRO_REDUX.
;   public      - If set with /PUBLIC, then build links to any files with
;                 anything in the PUBLIC field of the plate list.
;                 If set to a string, then select those plates that contain
;                 the substring PUBLIC within their PUBLIC field.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   IDL> platelinks, '/u/schlegel/spectro_EDR', public='EDR'
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   djs_filepath()
;   fileandpath()
;   platelist
;   splog
;
; REVISION HISTORY:
;   21-Jan-2003  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platelinks, outdir, public=public

   indir = getenv('BOSS_SPECTRO_REDUX')
   if (NOT keyword_set(outdir)) then begin
      splog, 'OUTDIR must be specified'
      return
   endif

   ;----------
   ; Get the list of plates

   platelist, plist=plist
   if (NOT keyword_set(plist)) then begin
      splog, 'No plates in platelist file'
      return
   endif

   if (keyword_set(public)) then begin
      if (size(public,/tname) EQ 'STRING') then begin
         itrim = where(strmatch(plist.public,'*'+public+'*'), nplate)
         sumfiles = ['spAll-'+public+'.fits', 'spAll-'+public+'.dat', $
          'spAllLine-'+public+'.fits']
      endif else begin
         itrim = where(strtrim(plist.public) NE '', nplate)
         sumfiles = ['spAll-public.fits', 'spAll-public.dat', $
          'spAllLine-public.fits']
      endelse
      if (nplate EQ 0) then begin
         splog, 'No plates matching PUBLIC keyword'
         return
      endif
      plist = plist[itrim]
   endif else begin
      sumfiles = ['spAll.fits', 'spAll.dat', 'spAllLines.dat']
   endelse

   ;----------
   ; Build links to the summary files

   for ifile=0, n_elements(sumfiles)-1 do begin
      file1 = djs_filepath(sumfiles[ifile], root_dir=indir)
      file2 = djs_filepath(sumfiles[ifile], root_dir=outdir)
      junk = findfile(file1, count=ct)
      if (ct GT 0) then begin
         junk = fileandpath(file2, path=thisdir)
         spawn, ['mkdir', '-p', thisdir]
         spawn, ['ln', '-s', file1, file2]
      endif
   endfor

   ;----------
   ; Loop over each plate, and build links

   nplate = n_elements(plist)
   for iplate=0, nplate-1 do begin
      platestr = string(plist[iplate].plate, format='(i4.4)')
      mjdstr = string(plist[iplate].mjd, format='(i5.5)')
      files = [ $
       'spPlate-'+platestr+'-'+mjdstr+'.fits', $
       'spZbest-'+platestr+'-'+mjdstr+'.fits', $
       'spZline-'+platestr+'-'+mjdstr+'.fits', $
       'spZall-'+platestr+'-'+mjdstr+'.fits' ]
      for ifile=0, n_elements(files)-1 do begin
         file1 = djs_filepath(files[ifile], root_dir=indir, subdir=platestr)
         file2 = djs_filepath(files[ifile], root_dir=outdir, subdir=platestr)
         junk = findfile(file1, count=ct)
         if (ct GT 0) then begin
            junk = fileandpath(file2, path=thisdir)
            spawn, ['mkdir', '-p', thisdir], /noshell
            spawn, ['ln', '-s', file1, file2], /noshell
         endif
      endfor
   endfor

   return
end
;------------------------------------------------------------------------------
