;+
; NAME:
;   mysql_write
;
; PURPOSE:
;   Write an IDL structure to a MySQL-readable file
;
; CALLING SEQUENCE:
;   mysql_write, filename, pdata, [table=, delim=, modifiers=, /append, $
;    scriptfile= ]
;
; INPUTS:
;   filename   - Output file name
;   pdata      - Structure to write
;
; OPTIONAL INPUTS:
;   table      - Name of table; default to the structure name or 'anonymous'
;   delim      - MySQL delimiter; default to the tab character
;   modifiers  - Modifiers for all variables, such as 'not null';
;                default to ''
;   append     - If set, then append to an existing file
;   scriptfile - Name of file for writing MySQL importing script
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Unsupported IDL variable types: COMPLEX, DCOMPLEX.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   22-Apr-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro mysql_write, filename, pdata, table=table, delim=delim, $
 modifiers=modifiers, append=append, scriptfile=scriptfile

   if (NOT keyword_set(delim)) then delim = string(9B) ; tab character
   if (n_elements(modifiers) EQ 0) then modifiers = '';

   tname = ['MEDIUMTEXT', 'BINARY(1)', $
    'INT2', 'INT4', 'INT8', $
    'INT2 UNSIGNED', 'INT4 UNSIGNED', 'INT8 UNSIGNED', $
    'FLOAT4', 'FLOAT8']
   idlname = ['STRING', 'BYTE', $
    'INT', 'LONG', 'LONG64', $
    'UINT', 'ULONG', 'ULONG64', $
    'FLOAT', 'DOUBLE']

   reserved_names = ['dec']
   reserved_replace = ['decc']

   ntag = N_tags( pdata )
   tags = strlowcase( tag_names( pdata ) )
   if (keyword_set(table)) then stname = table $
    else stname = strlowcase( tag_names( pdata, /structure_name) )
   if (stname EQ '') then stname = 'anonymous'

   ; Replace variable names that are reserved by MySQL to non-conflicting names
   for itag=0L, ntag-1 do begin
      ii = where(tags[itag] EQ reserved_names)
      if (ii[0] NE -1) then tags[itag] = reserved_replace[ii[0]]
   endfor

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
      nel = n_elements( pdata[0].(itag) )
      tagname = tname[(where(idlname EQ tt))[0]]

      if (itag EQ 0) then sline = '' $
       else sline = ' '

      if (nel EQ 1) then begin
         sline = sline + tags[itag] + ' ' + tagname + ' ' + modifiers
      endif else begin
         varname = tags[itag]
         for j=0L, ndim-1 do $
          varname = varname + '_' + $
           strtrim(string( djs_laxisgen(dims, iaxis=j) ), 2)
         for j=0L, nel-1 do begin
            sline = sline + varname[j] + ' ' + tagname + ' ' + modifiers
            if (j NE nel-1) then sline = sline + ', '
         endfor
      endelse

      if (itag EQ ntag-1) then sline = sline + ');' $
       else sline = sline + ','

      printf, olun, sline
   endfor

   printf, olun, 'load data local infile ' + string(39B) + filename $
    + string(39B) +' into table ' + stname + ';'

   if (keyword_set(scriptfile)) then begin
      close, olun
      free_lun, olun
   endif

   ;---------------------------------------------------------------------------
   ; Write the data to a MySQL-importable file
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
            sline = sline + words[0]
            for i=1L, N_elements( (pdata)[iel].(itag) )-1 do $
             sline = sline + delim + words[i]
         endelse
      endfor
      printf, olun, sline
   endfor

   close, olun
   free_lun, olun

   return
end
;------------------------------------------------------------------------------
