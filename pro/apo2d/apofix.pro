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
;   card       - FITS header keyword to change.
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
;
; BUGS:
;
; PROCEDURES CALLED:
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
;------------------------------------------------------------------------------
; Option for /good, /bad (but not both), setting all needed keywords.
; CARD,VALUE as arrays???
; Report warning if not all 4 files (cameras) exist.
;    Compare values in files to those requested, and add line to op file
;    only if different.
; Another proc that prints all required/optional header keywords.
; Verify that inputs are scalars when they need to be.
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
; ???
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
   if (expnum LE 0 OR expnum GT 99999999L) then begin
      print, 'EXPNUM must be a number between 1 and 99999999'
      !quiet = quiet
      return
   endif

   if (NOT keyword_set(camera)) then camera = '??'
   c1 = strmid(camera,0,1)
   c2 = strmid(camera,1,1)
   if (size(camera, /tname) NE 'STRING' $
    OR strlen(camera) NE 2 $
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
    OR strlen(card) EQ 0 OR strlen(card) GT 8) then begin
      print, 'CARD must be a string of 1 to 8 characters'
      !quiet = quiet
      return
   endif

   if (n_elements(newval) EQ 0) then begin
      print, 'VALUE must be specified'
      !quiet = quiet
      return
   endif
   strval = strtrim(string(newval),2)

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
      print, 'Reading FITS header for ' + filename[ifile]
      qdone = fits_wait(filename[ifile], deltat=2, tmax=10, /header_only)
      if (qdone) then begin
         thishdr = headfits(filename[ifile])
         thismjd = sxpar(thishdr, 'MJD')
         oldval = sxpar(thishdr, card)
; ???
;         if (size(oldval, /tname) EQ 'STRING') then begin
;            qdiff = strtrim(oldval,2) NE strval
;         endif else begin
;            qdiff = oldval NE newval
;         endelse
;print,oldval,newval,qdiff
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
   ; Append to the existing data structures in the sdHdrFix file,
   ; or create a new one if the file does not yet exist.

   if (NOT keyword_set(pdata)) then begin
      while (djs_lockfile(sdfixname) EQ 0) do wait, 2
      yanny_write, sdfixname, ptr_new(fstruct)
      djs_unlockfile, sdfixname
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
   endelse

   yanny_free, pdata
   !quiet = quiet

   return
end
;------------------------------------------------------------------------------
