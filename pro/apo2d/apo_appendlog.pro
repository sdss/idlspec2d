;+
; NAME:
;   apo_appendlog
;
; PURPOSE:
;   Append to logfile as written by APOREDUCE.
;
; CALLING SEQUENCE:
;   apo_appendlog, logfile, rstruct, tstruct
;
; INPUTS:
;   logfile    - FITS logfile as written by APOREDUCE.
;   rstruct    - Structure to append to the log file with reduced data info
;   tstruct    - Structure to append to the log file with WARNING/ABORT strings
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
pro apo_appendlog, logfile, rstruct, tstruct

   ;----------
   ; Determine which HDU in the log file this structure will be appended.

   if (keyword_set(rstruct)) then begin
      case rstruct.flavor of
         'bias'    : thishdu = 1
         'dark'    : thishdu = 1
         'flat'    : thishdu = 2
         'arc'     : thishdu = 3
         'science' : thishdu = 4
         'smear'   : thishdu = 4
         else      : thishdu = -1
      endcase
   endif else begin
      thishdu = -1
   endelse

   ;----------
   ; Lock the file to do this - otherwise we might read/write to a partially
   ; written file.

   while(djs_lockfile(logfile) EQ 0) do wait, 1

   ;----------
   ; If the log file does not yet exist, then create it.  Otherwise,
   ; append this structure to an existing structure, if it already exists.

   hdr = headfits(logfile)
   if (NOT keyword_set(hdr) OR size(hdr,/tname) NE 'STRING') then begin
      ; Create a new FITS file

      ; Write HDU#0, which  is just a header with the version of the code.
      mkhdr, newhdr, '', /extend
      sxaddpar, newhdr, 'VERS2D', idlspec2d_version()
      writefits, logfile, 0, newhdr

      ; Write HDU numbers 1 through 5
      for ihdu=1, 5 do begin
         if (ihdu EQ thishdu) then $
          mwrfits, rstruct, logfile $
         else $
          mwrfits, dummy, logfile
      endfor
   endif else if (thishdu GT 0) then begin
      ; Modify an existing FITS file
      pp = mrdfits(logfile, thishdu)
      pp = struct_append(pp, rstruct)
      djs_modfits, logfile, pp, exten_no=thishdu
   endif

   ;----------
   ; If TSTRUCT is set, then append that to HDU #5

   if (keyword_set(tstruct)) then begin
      pp = mrdfits(logfile, 5)
      pp = struct_append(pp, tstruct)
      djs_modfits, logfile, pp, exten_no=5
   endif

   ;----------
   ; Now unlock the log file.

   djs_unlockfile, logfile
   return
end
;------------------------------------------------------------------------------
