;+
; NAME:
;   splog
;
; PURPOSE:
;   Logging routine for spectroscopic pipeline.
;
; CALLING SEQUENCE:
;   splog, v1, v2 ..., [, _EXTRA=extra, /noname, filename=filename, $
;    /append, /close ]
;
; INPUTS:
;   v1, v2 ... - The expressions to be passed to PRINT or PRINTF
;
; OPTIONAL KEYWORDS:
;   _EXTRA     - Any keywords for PRINT or PRINTF, such as FORMAT
;   noname     - If set, then suppress printing name of calling routine
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

pro splog, noname=noname, filename=filename, append=append, close=close, $
 v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, _EXTRA=extra

   ; Declare LOGLUN in a common block so that it is remembered between calls.
   common com_splog, loglun

   if (keyword_set(filename)) then begin
      ; First close a file if one is already open
      if (keyword_set(loglun)) then splog, /close

      ; Now open the file
      get_lun, loglun
      openw, loglun, filename, append=append
   endif

   ; Determine the name of the calling routine
   help, calls=calls
   fname = (str_sep(calls[1], ' '))[0] + ': '

   nv = N_params()
   if (nv GT 0 OR keyword_set(extra)) then begin
      if (NOT keyword_set(noname)) then print, fname, format='(a,$)'

      case nv of
      0: print, _EXTRA=extra
      1: print, v1, _EXTRA=extra
      2: print, v1, v2, _EXTRA=extra
      3: print, v1, v2, v3, _EXTRA=extra
      4: print, v1, v2, v3, v4, _EXTRA=extra
      5: print, v1, v2, v3, v4, v5, _EXTRA=extra
      6: print, v1, v2, v3, v4, v5, v6, _EXTRA=extra
      7: print, v1, v2, v3, v4, v5, v6, v7, _EXTRA=extra
      8: print, v1, v2, v3, v4, v5, v6, v7, v8, _EXTRA=extra
      9: print, v1, v2, v3, v4, v5, v6, v7, v8, v9, _EXTRA=extra
      10: print, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, _EXTRA=extra
      11: print, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, _EXTRA=extra
      else: print, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, _EXTRA=extra
      endcase
   endif

   if ((nv GT 0 OR keyword_set(extra)) AND keyword_set(loglun)) then begin
      if (NOT keyword_set(noname)) then printf, loglun, fname, format='(a,$)'

      case nv of
      0: printf, loglun, _EXTRA=extra
      1: printf, loglun, v1, _EXTRA=extra
      2: printf, loglun, v1, v2, _EXTRA=extra
      3: printf, loglun, v1, v2, v3, _EXTRA=extra
      4: printf, loglun, v1, v2, v3, v4, _EXTRA=extra
      5: printf, loglun, v1, v2, v3, v4, v5, _EXTRA=extra
      6: printf, loglun, v1, v2, v3, v4, v5, v6, _EXTRA=extra
      7: printf, loglun, v1, v2, v3, v4, v5, v6, v7, _EXTRA=extra
      8: printf, loglun, v1, v2, v3, v4, v5, v6, v7, v8, _EXTRA=extra
      9: printf, loglun, v1, v2, v3, v4, v5, v6, v7, v8, v9, _EXTRA=extra
      10: printf, loglun, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, _EXTRA=extra
      11: printf, loglun, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, _EXTRA=extra
      else: printf, loglun, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, _EXTRA=extra
      endcase
      flush, loglun
   endif

   if (keyword_set(close) AND keyword_set(loglun)) then begin
      close, loglun
      free_lun, loglun
      loglun = 0
   endif

   return
end
