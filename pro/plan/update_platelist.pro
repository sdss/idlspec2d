;+
; NAME:
;   update_platelist
;
; PURPOSE:
;   Update the plate list file found in this idlspec2d product.
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

   for i=0, n_elements(thislist)-1 do begin
      j = where(newlist.plate EQ thislist[i].plate $
       AND newlist.mjd EQ thislist[i].mjd)
      copy_struct_inx, thislist[i], newlist, index_to=j
   endfor

   yanny_write, outfile, ptr_new(newlist), $
    structs=structs, stnames='spplatelist'

   return
end
;------------------------------------------------------------------------------
