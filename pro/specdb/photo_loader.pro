;+
; NAME:
;   photo_loader
;
; PURPOSE:
;   Write PHOTO outputs from SDSS_READOBJ() to an ASCII stream for direct input
;   into Postgres
;
; CALLING SEQUENCE:
;   photo_loader, runnum, camcol, rerun=, [ table=, /create, /clobber ]
;
; INPUTS:
;   runnum     - Run number
;   camcol     - Camera column
;   rerun      - Rerun name
;
; OPTIONAL INPUTS:
;   table      - Database table name; default to "photo"
;   create     - If set, then assume that the database definition
;                does not exist and create it
;   clobber    - If set, then clobber any existing table by this name
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
;   echo "photo_loader,3511,1,rerun=131,/create" | idl | psql foo
;
; BUGS:
;
; PROCEDURES CALLED:
;   sdss_readobj()
;
; REVISION HISTORY:
;   02-May-2003  Written by D. Schlegel and D. Hogg, Princeton.
;-
;------------------------------------------------------------------------------
pro photo_loader, runnum, camcol, rerun=rerun, table=stname, create=create, $
 clobber=clobber

   if (NOT keyword_set(stname)) then stname = 'photo'
   delim = string(9B) ; tab character
   select_tags = ['RUN', 'RERUN', 'CAMCOL', 'FIELD', 'ID', 'PARENT', $
    'RA', 'DEC', $
    'OBJC_TYPE', 'OBJC_FLAGS*', 'FLAGS*', $
    'PSP_STATUS', 'PSF_FWHM', 'AIRMASS', 'EXTINCTION', 'CLOUDCAM', $
    'PETRORAD*', 'PETROR50*', 'PETROR90*', $
    'MODELFLUX*', 'PSFFLUX*', 'PETROFLUX*', $
    'RUN_GLOBAL_STATUS', 'RUN_LOCAL_STATUS']

   tname = ['TEXT', 'INT2', 'INT2', 'INT4', 'INT8', 'FLOAT4', 'FLOAT8']
   idlname = ['STRING', 'BYTE', 'INT', 'LONG', 'LONG64', 'FLOAT', 'DOUBLE']

   ;---------------------------------------------------------------------------
   ; Read the imaging data
   ;---------------------------------------------------------------------------

   pdata = sdss_readobj(runnum, camcol, $
    rerun=rerun, select_tags=select_tags, /silent)
   ntag = N_tags( pdata )
   tags = strlowcase( tag_names( pdata ) )

   ;---------------------------------------------------------------------------
   ; Write the table definition
   ;---------------------------------------------------------------------------

   if (keyword_set(create) OR keyword_set(clobber)) then begin

      if (keyword_set(clobber)) then print, 'drop table ' + stname + ';'

      print, 'create table ' + stname
      print, format='("(",$)'

      for itag=0, ntag-1 do begin          ; Loop through each variable
         tt = size( pdata[0].(itag), /tname )
         dims = size( pdata[0].(itag), /dimens )
         ndim = size( pdata[0].(itag), /n_dimen )
         tagname = tname[(where(idlname EQ tt))[0]]
   
         if (itag EQ 0) then sline = '' $
          else sline = ' '
         sline = sline + tags[itag] + ' ' + tagname
         for j=0, ndim-1 do begin
            sline = sline + '[' + strtrim(string(dims[j]),2) + ']'
         endfor
         if (itag EQ ntag-1) then sline = sline + ');' $
          else sline = sline + ','

         print, sline
      endfor
   endif

   print, 'copy ' + stname + ' from stdin;'

   ;---------------------------------------------------------------------------
   ; Write the data as a Postgres-importable ASCII stream
   ;---------------------------------------------------------------------------

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
      print, sline
   endfor

   print, '\.' ; End of data-input

   flush, -1 ; Bug in IDL that does not flush devices on exit

   return
end
;------------------------------------------------------------------------------
