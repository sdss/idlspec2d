;+
; NAME:
;   specheaders
;
; PURPOSE:
;   Consolodate all of the raw spectro files headers into a single FITS file
;
; CALLING SEQUENCE:
;   specheaders
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
;   consolodates into a single FITS binary table file 'specheaders.fits'.
;   A specific file is chosen as the template file, and only header cards
;   that match that file are used.
;   Header keywords with '-' are replaced with '_'.
;   Repeat instances of the same keyword in the same file are ignored.
;   The 'END' card is ignored.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   fileandpath()
;   headfits()
;   sxpar()
;
; REVISION HISTORY:
;   02-Sep-2009  Written by David Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro specheaders

   ;----------
   ; Find all the b1 spFrame files, and sort by exposure number

   spawn, 'ls -d ' + getenv('BOSS_SPECTRO_DATA') + '/?????', mjdlist
;mjdlist = mjdlist[2500:2550] ; ???
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
   ; Construct a structure with all the header info
   ; Use a specific file as a template

   file_template = files[(where(strmatch(files,'*sdR-b1-00061798.fit*')))[0]]
   hdr = headfits(file_template)
   struct1 = create_struct('SIMPLE', 'T')
   hname = 'SIMPLE'
   for j=1L, n_elements(hdr)-1L do begin
      hname1 = strtrim(strmid(hdr[j],0,8),2)
      thisname = repstr(hname1,'-','_')
      if (thisname EQ 'NE') then thisname = 'NEON'
      if (hname1 NE '' AND hname1 NE 'END' AND $
       total(hname1 EQ hname) EQ 0) then begin
          struct1 = create_struct(struct1, thisname, (sxpar(hdr,thisname))[0])
         hname = [hname, hname1] ; Header names
      endif
   endfor

   structs = replicate(struct1, nfile)
   tags = tag_names(struct1) ; Tag names

   for ifile=0L, nfile-1L do begin
      splog, 'Reading file ', ifile, ' of ', nfile, $
       ': ', fileandpath(files[ifile])
      hdr = headfits(files[ifile])
      for j=0L, n_elements(hname)-1L do $
       structs[ifile].(j) = sxpar(hdr, hname[j])
   endfor

   mwrfits, structs, 'specheaders.fits' ,/create

   return
end
;------------------------------------------------------------------------------
