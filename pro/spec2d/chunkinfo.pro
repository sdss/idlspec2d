;+
; NAME:
;   chunkinfo
;
; PURPOSE:
;   Return chunk information for a given plate number.
;
; CALLING SEQUENCE:
;   cinfo = chunkinfo(plateid)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   plateid     - Plate ID number; scalar or array of integers
;
; OUTPUTS:
;   cinfo       - Structure with chunk information for each plate
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/spChunkList.par
;
; PROCEDURES CALLED:
;   yanny_free
;   yanny_read
;
; REVISION HISTORY:
;   08-Feb-2001  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function chunkinfo, plateid

   nplate = n_elements(plateid)

   ;----------
   ; Read the data file with the chunk information

   chunkfile = filepath('spChunkList.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   yanny_read, chunkfile, pdata
   chunkdata = *pdata[0]
   yanny_free, pdata
   nchunk = n_elements(chunkdata)

   ;----------
   ; Create output structure, setting all information to blanks.
   ; For plates with no chunk info, return this blank information.

   retval = chunkdata[0]
   struct_assign, {junk:0}, retval
   retval = replicate(retval, nplate)

   ;----------
   ; Find the chunk data for each plate

   for iplate=0, nplate-1 do begin
      for ichunk=0, nchunk-1 do begin
         if (plateid[iplate] GE chunkdata[ichunk].plate_start $
          AND plateid[iplate] LE chunkdata[ichunk].plate_end) then $
           retval[iplate] = chunkdata[ichunk]
      endfor
   endfor

   return, retval
end
;------------------------------------------------------------------------------
