;+
; NAME:
;   apo_appendlog
;
; PURPOSE:
;   Append to logfile as written by APOREDUCE.
;
; CALLING SEQUENCE:
;   apo_appendlog, logfile, rstruct
;
; INPUTS:
;   logfile    - FITS logfile as written by APOREDUCE.
;   rstruct    - Structure to append to the log file.
;
; OPTIONAL INPUTS:
;
; OUTPUT:
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
;   djs_modfits
;   djs_unlockfile
;   headfits()
;   idlspec2d_version()
;   modfits
;   mrdfits()
;   mwrfits
;   splog
;   sxaddpar
;   struct_append()
;
; REVISION HISTORY:
;   02-Dec-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro apo_appendlog, logfile, rstruct

   ;----------
   ; Determine which HDU in the log file this structure will be appended.

   case rstruct.flavor of
      'bias'    : thishdu = 1
      'dark'    : thishdu = 1
      'flat'    : thishdu = 2
      'arc'     : thishdu = 3
      'science' : thishdu = 4
      'smear'   : thishdu = 4
      else: begin
         splog, 'Unknown structure ', stname
         return
      end
   endcase

   ;----------
   ; Lock the file to do this - otherwise we might read/write to a partially
   ; written file.

   while(djs_lockfile(logfile) EQ 0) do wait, 1

   ;----------
   ; If the log file does not yet exist, then create it.  Otherwise,
   ; append this structure to an existing structure, if it already exists.

   pp = mrdfits(logfile, thishdu, hdr)
   if (NOT keyword_set(hdr)) then begin
      ; Create a new FITS file
      for ihdu=1, 4 do begin
         if (ihdu EQ thishdu) then $
          mwrfits, rstruct, logfile, create=(ihdu EQ 1) $
         else $
          mwrfits, dummy, logfile, create=(ihdu EQ 1)
      endfor
      ; Write the version of the code into the first header
      hdr = headfits(logfile)
      sxaddpar, hdr, 'VERS2D', idlspec2d_version()
      modfits, logfile, 0, hdr, exten_no=0
   endif else begin
      ; Modify to an existing FITS file
      pp = struct_append(pp, rstruct)
      djs_modfits, logfile, pp, exten_no=thishdu
   endelse

   ;----------
   ; Now unlock the log file.

   djs_unlockfile, logfile
   return
end
;------------------------------------------------------------------------------
