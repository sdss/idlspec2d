;+
; NAME:
;   apo_readlog
;
; PURPOSE:
;   Read a FITS logfile as written by APOREDUCE.
;
; CALLING SEQUENCE:
;   pstruct = apo_readlog, logfile, [ plate=, flavor=, camera= ]
;
; INPUTS:
;   logfile    - Logfile as written by APOREDUCE.  This is a FITS file
;                with an HDU of information for each reduced frame.
;
; OPTIONAL INPUTS:
;   plate      - Return only data for frames matching this plate.
;   flavor     - Return only data for frames matching this flavor.
;   camera     - Return only data for frames matching this camera.
;
; OUTPUT:
;   pstruct    - Array of pointer to structures for data in each reduced frame;
;                reference the first structure with (*pstruct[0]).
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_lockfile()
;   djs_unlockfile
;   mrdfits()
;
; REVISION HISTORY:
;   02-May-2000  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
function apo_readlog, logfile, plate=plate, flavor=flavor, camera=camera

   ; Lock the file to do this - otherwise we might read a partially
   ; written file.
   while(djs_lockfile(logfile, lun=html_lun) EQ 0) do wait, 1

   ; Read in all the HDU's in the log file as structures.
   ; Only include plate numbers and flavors matching PLATE, FLAVOR, CAMERA
   ; if those keywords are specified.

   pstruct = 0
   ihdu = 1
   pp = 1
   while (keyword_set(pp)) do begin
      pp = mrdfits(logfile, ihdu)
      if (keyword_set(pp)) then begin
         append = 1
         if (keyword_set(plate)) then if (pp.plate NE plate) then append = 0
         if (keyword_set(flavor)) then if (pp.flavor NE flavor) then append = 0
         if (keyword_set(camera)) then if (pp.camera NE camera) then append = 0
         if (append) then begin
            if (NOT keyword_set(pstruct)) then pstruct = ptr_new(pp) $
             else pstruct = [pstruct, ptr_new(pp)]
         endif
      endif
      ihdu = ihdu + 1
   endwhile
   nstruct = n_elements(pstruct)

   ; Now unlock the log file.
   djs_unlockfile, logfile, lun=html_lun

   return, pstruct
end
