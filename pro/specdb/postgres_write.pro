;+
; NAME:
;   postgres_write
;
; PURPOSE:
;   Write an IDL structure to a Postgres-readable file
;
; CALLING SEQUENCE:
;   postgres_write, filename, pdata, [table=, delim=, /append, $
;    scriptfile= ]
;
; INPUTS:
;   filename   - Output file name
;   pdata      - Structure to write
;
; OPTIONAL INPUTS:
;   table      - Name of table; default to the structure name or 'anonymous'
;   delim      - Postgres delimiter; default to the tab character
;   append     - If set, then append to an existing file
;   scriptfile - Name of file for writing Postgres importing script
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Unsupported IDL variable types: UINT, ULONG, ULONG64, COMPLEX, DCOMPLEX.
;   IDL type BYTE is converted to a 2-byte integer.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   06-Apr-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro postgres_write, filename, pdata, table=table, delim=delim, append=append, $
 scriptfile=scriptfile

   if (NOT keyword_set(delim)) then delim = string(9B) ; tab character

   tname = ['TEXT', 'INT2', 'INT2', 'INT4', 'INT8', 'FLOAT4', 'FLOAT8']
   idlname = ['STRING', 'BYTE', 'INT', 'LONG', 'LONG64', 'FLOAT', 'DOUBLE']

   ntag = N_tags( pdata )
   tags = strlowcase( tag_names( pdata ) )
   if (keyword_set(table)) then stname = table $
    else stname = strlowcase( tag_names( pdata, /structure_name) )
   if (stname EQ '') then stname = 'anonymous'

   ;---------------------------------------------------------------------------
   ; Write the table definition
   ;---------------------------------------------------------------------------

   if (keyword_set(scriptfile)) then begin
      get_lun, olun
      openw, olun, scriptfile
   endif else begin
      olun = -1
   endelse

   printf, olun, 'create table ' + stname
   printf, olun, format='("(",$)'

   for itag=0L, ntag-1 do begin          ; Loop through each variable
      tt = size( pdata[0].(itag), /tname )
      dims = size( pdata[0].(itag), /dimens )
      ndim = size( pdata[0].(itag), /n_dimen )
      tagname = tname[(where(idlname EQ tt))[0]]

      if (itag EQ 0) then sline = '' $
       else sline = ' '
      sline = sline + tags[itag] + ' ' + tagname
      for j=0L, ndim-1 do begin
         sline = sline + '[' + strtrim(string(dims[j]),2) + ']'
      endfor
      if (itag EQ ntag-1) then sline = sline + ');' $
       else sline = sline + ','

      printf, olun, sline
   endfor

   printf, olun, 'copy ' + stname + ' from ' + string(39B) + filename $
    + string(39B) + ';'

   if (keyword_set(scriptfile)) then begin
      close, olun
      free_lun, olun
   endif

   ;---------------------------------------------------------------------------
   ; Write the data to a Postgres-importable file
   ;---------------------------------------------------------------------------

   get_lun, olun
   openw, olun, filename, append=append

   for iel=0L, N_elements(pdata)-1 do begin ; Loop thru each row

      for itag=0L, ntag-1 do begin          ; Loop through each variable
         words = pdata[iel].(itag)
         nword = N_elements(words)

         ; If WORDS is type STRING, then check for white-space
         ; in any of its elements.  If there is white space, then
         ; put double-quotes around that element.
         if (size(words,/tname) EQ 'STRING') then begin
            for iw=0L, nword-1 do $
             if (strpos(words[iw],' ') NE -1) then $
              words[iw] = '"' + words[iw] + '"'
         endif else begin
            words = strtrim(string(words),2)
         endelse

         if (itag EQ 0) then sline = '' $
          else sline = sline + delim

         if (nword EQ 1) then begin
            sline = sline + words
         endif else begin
            sline = sline + '{' + words[0]
            for i=1L, N_elements( (pdata)[iel].(itag) )-1 do $
             sline = sline + ',' + words[i]
            sline = sline + '}'
         endelse
      endfor
      printf, olun, sline
   endfor

   close, olun
   free_lun, olun

   return
end
;------------------------------------------------------------------------------
