;+
; NAME:
;   splog
;
; PURPOSE:
;   Logging routine for spectroscopic pipeline.
;
; CALLING SEQUENCE:
;   splog, v1, v2 ..., [, _EXTRA=extra, /noname, prelog=, filename=, $
;    /append, /close ]
;
; INPUTS:
;   v1, v2 ... - The expressions to be passed to PRINT or PRINTF
;
; OPTIONAL KEYWORDS:
;   _EXTRA     - Any keywords for PRINT or PRINTF, such as FORMAT
;   noname     - If set, then suppress printing name of calling routine
;   prelog     - If set, then print this string immediately after the
;                name of the calling routine on each line, i.e. 'b1'
;   filename   - If specified, then open this file for text output
;   append     - If set at the same time as FILENAME, then append to this file;
;                default is to overwrite file
;   close      - If set, then close the output file
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The output is formatted just like with the IDL PRINT command, except
;   that extraneous whitespace is removed unless a FORMAT keyword is used.
;
; EXAMPLES:
;   Open a file for text output, write to it, then close it:
;     splog, filename='test.log'
;     splog, 'This is a test', 123.4, format='(a,f8.2)'
;     splog, /close
;   Alternatively, this can all be done on one line
;     splog, filename='test.log', /close, $
;      'This is a test', 123.4, format='(a,f8.2)'
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Nov-1999  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro splog, noname=noname, prelog=prelog, $
 filename=filename, append=append, close=close, $
 v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, _EXTRA=extra

   ; Declare LOGLUN in a common block so that it is remembered between calls.
   common com_splog, loglun, fullprelog

   if (keyword_set(filename)) then begin
      ; First close a file if one is already open
      if (keyword_set(loglun)) then splog, /close

      ; Now open the file
      get_lun, loglun
      openw, loglun, filename, append=append
   endif

   if (N_elements(prelog) EQ 1) then fullprelog = prelog

   ; Determine the name of the calling routine
   help, calls=calls

   fname = (str_sep(calls[1], ' '))[0] + ': '

   ; Add spaces for depth of routine
   for i=0,n_elements(calls)-4 do fname = ' ' + fname 

   ;----------
   ; Construct the output text string

   nv = N_params()
;   if (nv GT 0 OR keyword_set(extra)) then begin
   if (nv GT 0) then begin
      case nv of
;      0: textstring = string('', _EXTRA=extra) ; Does not work! IDL bug?
      1: textstring = string(v1, _EXTRA=extra)
      2: textstring = string(v1, v2, _EXTRA=extra)
      3: textstring = string(v1, v2, v3, _EXTRA=extra)
      4: textstring = string(v1, v2, v3, v4, _EXTRA=extra)
      5: textstring = string(v1, v2, v3, v4, v5, _EXTRA=extra)
      6: textstring = string(v1, v2, v3, v4, v5, v6, _EXTRA=extra)
      7: textstring = string(v1, v2, v3, v4, v5, v6, v7, _EXTRA=extra)
      8: textstring = string(v1, v2, v3, v4, v5, v6, v7, v8, _EXTRA=extra)
      9: textstring = string(v1, v2, v3, v4, v5, v6, v7, v8, v9, _EXTRA=extra)
      10: textstring = string(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, _EXTRA=extra)
      11: textstring = string(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, _EXTRA=extra)
      else: textstring = string(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, _EXTRA=extra)
      endcase
   endif

   ;----------
   ; If there is no FORMAT keyword specified, then compress the output string
   ; (remove extraneous whitespace).

   qcompress = 1
   if (keyword_set(extra)) then begin
      tags = tag_names(extra)
      if ( (where(tags EQ 'FORMAT'))[0] NE -1) then qcompress = 0
   endif
   if (qcompress AND keyword_set(textstring)) then $
    textstring = strcompress(textstring)

   ;----------
   ; Write to standard out

   if (nv GT 0 OR keyword_set(extra)) then begin
      if (NOT keyword_set(noname)) then print, fname, format='(a,$)'
      if (keyword_set(fullprelog)) then $
       print, fullprelog+': ', format='(a,$)'
      if (nv EQ 0) then print, _EXTRA=extra $
       else print, textstring
   endif

   ;----------
   ; Write to the output file (if it exists)

   if ((nv GT 0 OR keyword_set(extra)) AND keyword_set(loglun)) then begin
      if (NOT keyword_set(noname)) then printf, loglun, fname, format='(a,$)'
      if (keyword_set(fullprelog)) then $
       printf, loglun, fullprelog+': ', format='(a,$)'
      if (nv EQ 0) then printf, loglun, _EXTRA=extra $
       else printf, loglun, textstring
      flush, loglun
   endif

   ;----------
   ; Close the output file

   if (keyword_set(close) AND keyword_set(loglun)) then begin
      close, loglun
      free_lun, loglun
      loglun = 0
      fullprelog = ''
   endif

   return
end
