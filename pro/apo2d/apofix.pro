;+
; NAME:
;   apofix
;
; PURPOSE:
;   Add line to sdHdrFix file to denote change in FITS header for sdR files.
;
; CALLING SEQUENCE:
;   apofix, expnum, [ card, value, camera=, /bad ]
;
; INPUTS:
;   expnum     - Exposure number
;
; OPTIONAL INPUTS:
;   card       - FITS header keyword to change; this is case-insensitive,
;                so that 'exposure' is the same as 'EXPOSURE'.
;   value      - New value for FITS header keyword.
;   camera     - Camera name in which to change values, e.g. 'b1', 'r1',
;                'b2' or 'r2'.  A '?' can be used as a wildcard, for example
;                '?2' to denote a change to both 'b2' and 'r2'.  Default to
;                '??' to denote a change to sdR files for all 4 cameras.
;   bad        - If set, then declare the specified exposure number to be bad.
;                Do this by setting FLAVOR='unknown'.
;                If set, then overwrite and values passed for CARD and VALUE.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   Fix the exposure time for exposure #1234 to be 900 sec:
;     IDL> apofix, 1234, 'exposure', 900
;
;   Fix the TAI time for exposure #1234 to be 4.443852968d+09
;   (use the "d" notation for double-precision, even though it
;   will appear in the sdHdrFix file with an "e"):
;     IDL> apofix, 1234, 'TAI', 4.443852968d+09
;
;   Declare exposure #1234 as bad:
;     IDL> apofix, 1234, /bad
;   or equivalently:
;     IDL> apofix, 1234, 'flavor', 'unknown'
;
; BUGS:
;
; PROCEDURES CALLED:
;   fileandpath()
;   fits_wait
;   headfits()
;   struct_append
;   sxpar()
;   yanny_read
;
; REVISION HISTORY:
;   22-Apr-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro apofix, expnum, card, newval, camera=camera, bad=bad

   common apofix_com, apo_uname

   if (n_params() LT 1) then begin
      print, 'Syntax - apofix, expnum, [ card, value, camera= ]'
      return
   end

   quiet = !quiet
   !quiet = 1

   ;----------
   ; Insist that this proc only run on machines named "sos"

   if (NOT keyword_set(apo_uname)) then begin
      spawn, 'uname -n', uname_string
      apo_uname = (strsplit(uname_string[0], '.', /extract))[0]
   endif

   if (apo_uname NE 'sos') then begin
;      print, 'This procedure can only be run on the machine sos.apo.nmsu.edu'
;      !quiet = quiet
;      return
   endif

   ;----------
   ; Set input directory for sdR files and sdHdrFix file

   rawdata_dir = getenv('RAWDATA_DIR')
   if (NOT keyword_set(rawdata_dir)) then $
    rawdata_dir = '/data/spectro'

   astrolog_dir = getenv('ASTROLOG_DIR')
   if (NOT keyword_set(astrolog_dir)) then $
    astrolog_dir = '/data/spectro/astrolog'

   ;----------
   ; Create output structure

   fstruct = create_struct( $
    name = 'OPHDRFIX', $
    'fileroot', '', $
    'keyword' , '', $
    'value'   , '' )

   ;----------
   ; Sanity checks on EXPNUM, CAMERA, CARD, VALUE

   expnum = long(expnum)
   if (expnum LE 0 OR expnum GT 99999999L OR n_elements(expnum) NE 1) then begin
      print, 'EXPNUM must be a number between 1 and 99999999'
      !quiet = quiet
      return
   endif

   if (NOT keyword_set(camera)) then camera = '??'
   c1 = strmid(camera,0,1)
   c2 = strmid(camera,1,1)
   if (size(camera, /tname) NE 'STRING' $
    OR strlen(camera) NE 2 $
    OR n_elements(camera) NE 1 $
    OR (c1 NE '?' AND c1 NE 'b' AND c1 NE 'r') $
    OR (c2 NE '?' AND c2 NE '1' AND c2 NE '2') ) then begin
      print, 'CAMERA must be a 2-character string'
      !quiet = quiet
      return
   endif

   if (keyword_set(bad)) then begin
      card = 'flavor'
      value = 'unknown'
   endif

   if (size(card, /tname) NE 'STRING' $
    OR n_elements(card) NE 1 $
    OR strlen(card) EQ 0 OR strlen(card) GT 8) then begin
      print, 'CARD must be a string of 1 to 8 characters'
      !quiet = quiet
      return
   endif

   if (n_elements(newval) NE 1) then begin
      print, 'VALUE must be specified (and a scalar)'
      !quiet = quiet
      return
   endif
   if (size(newval, /tname) EQ 'DOUBLE') then format='(e17.10)' $
    else format=''
   strval = strtrim(string(newval,format=format),2)
   if (strpos(strval,'"') NE -1 OR strpos(strval,"'") NE -1) then begin
      print, 'VALUE cannot contain single or double-quotes'
      !quiet = quiet
      return
   endif

   ;----------
   ; Test that sdR files exist that correspond to the exposure number
   ; and camera(s) specified.

   fileroot = string(camera, expnum, format='("sdR-",a2,"-",i8.8)')

   filename = findfile(filepath(fileroot+'.fit*', $
    root_dir=rawdata_dir, subdir='*'), count=nfile)
   if (nfile EQ 0) then begin
      print, 'File=' + fileroot + '.fit not found in dir=' + rawdata_dir
      !quiet = quiet
      return
   endif

   ;----------
   ; Read the values of the cards in these sdR files
   ; Also, compare values in sdR headers to requested values.

   for ifile=0, nfile-1 do begin
      qdone = fits_wait(filename[ifile], deltat=2, tmax=10, /header_only)
      if (qdone) then begin
         thishdr = headfits(filename[ifile])
         thismjd = sxpar(thishdr, 'MJD')
         oldval = sxpar(thishdr, card)
         print, strmid(fileandpath(filename[ifile]),0,15), strupcase(card), $
          strtrim(string(oldval,format=format),2), strval, $
          format='(a, 1x, a8, "=[", a, "] -> [", a, "]")'
      endif
   endfor

   ;----------
   ; Construct the output structure

   fstruct.fileroot = fileroot
   fstruct.keyword = strupcase(card)
   fstruct.value = strval

   ;----------
   ; Read the sdHdrFix file

   if (NOT keyword_set(thismjd)) then begin
      print, 'MJD could not be determined from the FITS headers'
      !quiet = quiet
      return
   endif

   mjdstr = string(thismjd, format='(i5.5)')
   sdfixname = filepath('sdHdrFix-' + mjdstr + '.par', root_dir=astrolog_dir, $
    subdir=mjdstr)

   while (djs_lockfile(sdfixname) EQ 0) do wait, 2
   yanny_read, sdfixname, pdata, hdr=hdr, enums=enums, structs=structs, $
    stnames=stnames, /anonymous
   djs_unlockfile, sdfixname

   ;----------
   ; Use an explicitly-defined structs such that the fTCL Yanny-readers
   ; will still work.

   structs = ['typedef struct {', $
              '  char fileroot[20]; # Root of file name, without any ".fit" suffix', $
              '  char keyword[9]; # Keyword name', $
              '  char value[80]; # Keyword value (as a string)', $
              '} OPHDRFIX;']

   ;----------
   ; Append to the existing data structures in the sdHdrFix file,
   ; or create a new one if the file does not yet exist.

   if (NOT keyword_set(pdata)) then begin
      while (djs_lockfile(sdfixname) EQ 0) do wait, 2
      yanny_write, sdfixname, ptr_new(fstruct), structs=structs
      djs_unlockfile, sdfixname
      ncorr = 1
   endif else begin
      ; Append data to the relevant data structure.
      i = (where(stnames EQ tag_names(fstruct, /structure_name)))[0]
      if (i[0] EQ -1) then begin
         print, 'The file ' + sdfixname + ' appears to be invalid'
         !quiet = quiet
         return
      endif
      pdata[i] = ptr_new(struct_append(*pdata[i], fstruct))

      while (djs_lockfile(sdfixname) EQ 0) do wait, 2
      yanny_write, sdfixname, pdata, hdr=hdr, enums=enums, structs=structs, $
       stnames=stnames
      djs_unlockfile, sdfixname
      ncorr = n_elements(*pdata[i])
   endelse

   yanny_free, pdata
   !quiet = quiet

   print, 'File ' + fileandpath(sdfixname) + ' contains ' $
    + strtrim(string(ncorr),2) + ' declared changes.'

   return
end
;------------------------------------------------------------------------------
