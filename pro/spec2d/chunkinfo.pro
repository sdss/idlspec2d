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
;   $PLATELIST_DIR/platePlans.par
;
; PROCEDURES CALLED:
;   yanny_readone()
;
; REVISION HISTORY:
;   08-Feb-2001  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function chunkinfo, plateid,legacy=legacy,plates=plates

   common com_chunkinfo, chunkdata

   nplate = n_elements(plateid)

   ;----------
   ; Read the data file with the chunk information

   if (NOT keyword_set(chunkdata)) then begin
      if keyword_set(plates) or keyword_set(legacy) then begin
         chunkfile = filepath('platePlans.par', root_dir=getenv('PLATELIST_DIR'));PLATELIST_DIR
      endif else begin
         chunkfile = filepath('configPlans.par', root_dir=getenv('CONFIGLIST_DIR'));PLATELIST_DIR
      endelse
      chunkdata = yanny_readone(chunkfile)
   endif
   if (NOT keyword_set(chunkdata)) then begin
      if keyword_set(plates) or keyword_set(legacy) then begin
         splog, 'Empty or missing platePlans.par file'
      endif else begin
         splog, 'Empty or missing configPlans.par file'
      endelse
      return, 0
   endif

   ;----------
   ; Create output structure, setting all information to blanks.
   ; For plates with no chunk info, return this blank information.

   retval = chunkdata[0]
   struct_assign, {junk:0}, retval
   retval = replicate(retval, nplate)

   ;----------
   ; Find the chunk data for each plate

   for iplate=0L, nplate-1L do begin
      if keyword_set(plates) or  keyword_set(legacy) then begin
         indx = where(chunkdata.plateid EQ plateid[iplate], ct)
      endif else begin
         indx = where(chunkdata.field EQ plateid[iplate], ct)
      endelse
      if (ct GT 0) then retval[iplate] = chunkdata[indx[0]]
   endfor

   return, retval
end
;------------------------------------------------------------------------------
