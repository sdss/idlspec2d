;+
; NAME:
;   update_plate_release
;
; PURPOSE:
;   Update the plate list file found in this idlspec2d product
;   with new data release notes.
;
; CALLING SEQUENCE:
;   update_plate_release, [ release, rfile, outfile ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   release    - Data release name to put into plate list file;
;                default to 'DR6'
;   rfile      - Yanny param file with PLATEID,MJD,RELEASE columns,
;                where the release names should appear in the RELEASE entries;
;                default to "$IDLSPEC2D_DIR/etc/dr2spectro.par",
;   outfile    - Name of output file; default to 'spPlateList.par' in the
;                current directory.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The file in $SPECLOG_DIR/opfiles/spPlateList.par is read, and we
;   merge in the list of plates returned by the PLATELIST procedure.
;   Manual comments from the first file are retained.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $SPECLOG_DIR/opfiles/spPlateList.par
;
; PROCEDURES CALLED:
;   yanny_readone()
;   yanny_write
;
; REVISION HISTORY:
;   04-Feb-2004  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro update_plate_release, release, rfile, outfile

   if (NOT keyword_set(release)) then release = 'DR6'
   if (NOT keyword_set(rfile)) then $
    rfile = filepath('dr2spectro.par', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdir='etc')
   if (NOT keyword_set(outfile)) then outfile = 'spPlateList.par'

   plist = yanny_readone(filepath('spPlateList.par', $
    root_dir=getenv('SPECLOG_DIR'), subdir='opfiles'), $
    hdr=hdr, enums=enums, structs=structs, stnames=stnames)
   drlist = yanny_readone(rfile)

;   indx = where(strmatch(strupcase(drlist.release), release), ct)
;   if (ct EQ 0) then begin
;      print, 'No matches found to RELEASE=', release
;   endif
;   drlist = drlist[indx]

   for i=0L, n_elements(drlist)-1 do begin
      j = where(plist.plate EQ drlist[i].plateid $
       AND plist.mjd EQ drlist[i].mjd, ct)
      if (ct NE 1) then $
       message, 'Not found: ' + string(drlist[i].plateid) $
        + ' ' + string(drlist[i].mjd)
      if (NOT keyword_set(plist[j].public)) then $
       plist[j].public = strtrim(plist[j].public+' '+release, 2)
   endfor

   yanny_write, outfile, ptr_new(plist), hdr=hdr, enums=enums, $
    structs=structs, stnames=stnames

   return
end
;------------------------------------------------------------------------------
