;+
; NAME:
;   update_platelist
;
; PURPOSE:
;   Update the plate list file found in the speclog product.
;
; CALLING SEQUENCE:
;   update_platelist, [ outfile ]
;
; INPUTS:
;   outfile    - Name of output file; default to 'spPlateList.par' in the
;                current directory.
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The file in $SPECLOG_DIR/opfiles/spPlateList.par is read, and we
;   merge in the list of plates returned by the PLATELIST procedure.
;   Manual comments from the first file are retained.
;   The auto-generated "platequality" is used for plates being added
;   if they didn't already have a quality in the file.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $SPECLOG_DIR/opfiles/spPlateList.par
;
; PROCEDURES CALLED:
;   copy_struct_inx
;   platelist
;   yanny_read
;   yanny_write
;
; REVISION HISTORY:
;   04-Jun-2002  Written by D. Schlegel, Princeton (checked in later)
;------------------------------------------------------------------------------
pro update_platelist, outfile

   if (NOT keyword_set(outfile)) then outfile = 'spPlateList.par'

   ;----------
   ; Read the existing list of plates

   thisfile = filepath('spPlateList.par', root_dir=getenv('SPECLOG_DIR'), $
    subdir='opfiles')
   yanny_read, thisfile, pdat, structs=structs, /anonymous
   thislist = *pdat[0]
   yanny_free, pdat

   ;----------
   ; Find out which plates are reduced

   platelist, plist=plist

   ;----------
   ; Merge these two lists

   platemjd1 = string(thislist.plate, format='(i4.4)') $
    + ' ' + string(thislist.mjd, format='(i5.5)')
   platemjd2 = string(plist.plate, format='(i4.4)') $
    + ' ' + string(plist.mjd, format='(i5.5)')
   platemjd = [platemjd1, platemjd2]
   platemjd = platemjd[ uniq(platemjd, sort(platemjd)) ]

   newlist = thislist[0]
   struct_assign, {junk:0}, newlist
   nplate = n_elements(platemjd)
   newlist = replicate(newlist, nplate)
   for i=0, nplate-1 do $
    newlist[i].plate = long( (strsplit(platemjd[i],/extract))[0] )
   for i=0, nplate-1 do $
    newlist[i].mjd = long( (strsplit(platemjd[i],/extract))[1] )

   ; Retain any values from the original spPlateList.par file
   for i=0L, n_elements(thislist)-1L do begin
      j = where(newlist.plate EQ thislist[i].plate $
       AND newlist.mjd EQ thislist[i].mjd)
      copy_struct_inx, thislist[i], newlist, index_to=j
   endfor

   ; Copy "platequality" from the auto-generated platelist file
   ; if there isn't already a value for that field
   for i=0L, n_elements(plist)-1L do begin
      j = (where(newlist.plate EQ plist[i].plate $
       AND newlist.mjd EQ plist[i].mjd, ct))[0]
      if (ct GT 0) then $
       if (strtrim(newlist[j].platequality) EQ '') then $
        newlist[j].platequality = plist[i].platequality
   endfor

   yanny_write, outfile, ptr_new(newlist), $
    structs=structs, stnames='spplatelist'

   return
end
;------------------------------------------------------------------------------
