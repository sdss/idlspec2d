;+
; NAME:
;   inspectfiles
;
; PURPOSE:
;   Create empty spInspect file(s) for manual inspection of classifications
;
; CALLING SEQUENCE:
;   inspectfiles, [ plate, mjd=, inspector=, /clobber ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   plate      - Plate number(s); if not set, then select all plates that
;                are reduced according to the platelist file.
;   mjd        - MJD number(s); if not set, then select the most recent
;                data for each plate (largest MJD).
;   inspector  - Name of inspector, which should be a string without
;                any whitespace.  If set, then output file is
;                "spInspect-$PLATE-$MJD-$INSPECTOR.par", otherwise it is
;                "spInspect-$PLATE-$MJD.par" with inspector='yourname'
;                in the inspector fields of this file.
;   clobber    - If set, then overwrite any files with the same name
;                in the current directory.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   Create spInspect files for all reduced plates and write them in
;   the current directory:
;     IDL> inspectfiles
;
;   Create the file spInspect-0400-51820-knapp.par file for plate 400 observed
;   on MJD=51820, filling in the "inspector" fields with the name "knapp":
;     IDL> inspectfiles, 400, mjd=51820, inspector='knapp'
;
; BUGS:
;
; PROCEDURES CALLED:
;   inspectgen
;   platelist
;   readspec
;   splog
;   yanny_write
;
; REVISION HISTORY:
;   20-May-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro inspectfiles, plate, mjd=mjd, inspector=inspector, clobber=clobber

   ;----------
   ; Check that INSPECTOR does not contain whitespace.

   if (keyword_set(inspector)) then begin
      if (inspector NE strcompress(inspector,/remove_all)) then begin
         print, 'INSPECTOR cannot contain whitespace.'
         return
      endif
   endif

   ;----------
   ; If PLATE is not set, then select all plates

   if (NOT keyword_set(plate)) then begin
      platelist, plist=plist
      indx = where(strtrim(plist.status1d,2) EQ 'Done')
      if (indx[0] EQ -1) then begin
         print, 'No reduced plates found'
         return
      endif
      plate = plist[indx].plate
      mjd = plist[indx].mjd
   endif

   ;----------
   ; If MJD is not set, then find the MJD for each plate

   nplate = n_elements(plate)
   if (NOT keyword_set(mjd)) then begin
      mjd = lonarr(nplate)
      for iplate=0, nplate-1 do begin
         mjd1 = 0
         readspec, plate[iplate], mjd=mjd1, topdir=topdir, /silent
         if (NOT keyword_set(mjd1)) then begin
            print, 'No MJD found for plate ', plate[iplate]
            return
         endif
         mjd[iplate] = mjd1
      endfor
   endif else begin
      if (n_elements(mjd) NE nplate) then begin
         print, 'Number of elements in PLATE and MJD do not agree'
         return
      endif
   endelse

   ;----------
   ; Loop over each plate and generate file

   for iplate=0, nplate-1 do begin
      specinspect = 0
      inspectgen, plate[iplate], mjd=mjd[iplate], inspector=inspector, $
       specinspect=specinspect, spectext=spectext, specblend=specblend, $
       hdr=hdr, structs=structs

      if (keyword_set(specinspect)) then begin
         if (keyword_set(inspector)) then $
          outfile = string(plate[iplate], mjd[iplate], inspector, $
           format='("spInspect-",i4.4,"-",i5.5,"-",a,".par")') $
         else $
          outfile = string(plate[iplate], mjd[iplate], $
           format='("spInspect-",i4.4,"-",i5.5,".par")')
         qexist = keyword_set(findfile(outfile))
         if (qexist) then begin
            if (keyword_set(clobber)) then $
             splog, 'WARNING: Over-writing file: ' + outfile $
            else $
             splog, 'WARNING: Will not over-write file: ' + outfile
         endif
         if ((NOT qexist) OR keyword_set(clobber)) then begin
            splog, 'Writing file ' + outfile
            yanny_write, outfile, [ptr_new(specinspect), ptr_new(spectext), $
             ptr_new(specblend)], hdr=hdr, structs=structs
         endif
      endif else begin
         splog, 'Skipping plate=', plate[iplate], ' MJD=', mjd[iplate]
      endelse
   endfor

   return
end
;------------------------------------------------------------------------------
