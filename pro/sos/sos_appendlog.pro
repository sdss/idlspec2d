;+
; NAME:
;   sos_appendlog
;
; PURPOSE:
;   Append to logfile as written by APOREDUCE.
;
; CALLING SEQUENCE:
;   sos_appendlog, logfile, rstruct, tstruct
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
;   idlutils_version()
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
pro sos_appendlog, logfile, rstruct, tstruct
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
    names = ['BIASDARK','FLAT','ARC','SCIENCE','TEXT']

   ;----------
   ; Lock the file to do this - otherwise we might read/write to a partially
   ; written file.

   ; If unable to get lock, then timeout after 100x2 sec
   qdone = 0
   iiter = 0L
   while(qdone EQ 0) do begin
      qdone = djs_lockfile(logfile)
      if (qdone EQ 0) then begin
         wait, 2
         if (iiter GE 100) then begin
            splog, 'ABORT: Timeout on lock file for '+logfile
            return
         endif
      endif
      iiter++
   end

   ;----------
   ; If the log file does not yet exist, then create it.  Otherwise,
   ; append this structure to an existing structure, if it already exists.

   hdr = headfits(logfile)
   if (NOT keyword_set(hdr) OR size(hdr,/tname) NE 'STRING') then begin
      ; Create a new FITS file

      ; Write HDU#0, which  is just a header with the version of the code.
      mkhdr, newhdr, '', /extend
      sxaddpar, newhdr, 'VERSIDL', !version.release, ' Version of IDL'
      sxaddpar, newhdr, 'VERSUTIL', idlutils_version()
      sxaddpar, newhdr, 'VERS2D', idlspec2d_version()
      sxaddpar, newhdr, 'RUN2D', getenv('IDLSPEC2D_VER')
      writefits, logfile, 0, newhdr

      ; Write HDU numbers 1 through 5
      for ihdu=1, 5 do begin
         if (ihdu EQ thishdu) then $
          mwrfits_named, rstruct, logfile, name=names[ihdu-1]$
         else $
          mwrfits_named, dummy, logfile, name=names[ihdu-1]
      endfor
   endif else if (thishdu GT 0) then begin
      ; Modify an existing FITS file
      pp = mrdfits(logfile, thishdu)

      ; check to see if this entry exists, if so overwrite:
      exists = -1
      if keyword_set(pp) then $
        exists = where(pp.expnum EQ rstruct.expnum AND $
                       pp.camera EQ rstruct.camera)

      if exists[0] EQ -1 then pp = struct_append(pp, rstruct) $
       else copy_struct_inx, rstruct, pp, index_to=exists[0]

      modfits_named, logfile, pp, exten_no=thishdu, name=names[thishdu-1]
   endif

   ;----------
   ; If TSTRUCT is set, then append that to HDU #5

   if (keyword_set(tstruct)) then begin
      pp = mrdfits(logfile, 5)
;      print, pp
;      help, pp
;      splog, typename(pp)

      ;-----------------------------------------------
      ; check to see if this entry exists, if so overwrite:
      exists = -1
      if keyword_set(pp) then $
        exists = where(strtrim(pp.filename,2) NE strtrim(tstruct[0].filename,2))

      splog, n_elements(tstruct)
      if (n_elements(tstruct) eq 1) and (strlen(strtrim(tstruct[0].text)) eq 0) then begin
          if exists[0] NE -1 then begin
              pp = pp[exists]
              modfits_named, logfile, pp, exten_no=5, name='text'
          endif else modfits_named, logfile, 0, exten_no=5, /delete_data
      endif else begin 
          if exists[0] NE -1 then begin
              pp = pp[exists] 
              pp = struct_append(pp, tstruct) 
          endif else pp = tstruct
          modfits_named, logfile, pp, exten_no=5, name='text'
      endelse

   endif

   ;----------
   ; Now unlock the log file.

   djs_unlockfile, logfile
   return
end
;------------------------------------------------------------------------------
