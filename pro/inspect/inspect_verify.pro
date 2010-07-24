;+
; NAME:
;   inspect_verify
;
; PURPOSE:
;   Verify the data format of spInspect files
;
; CALLING SEQUENCE:
;   inspect_verify, [ path= ]
;
; INPUTS:
;   path       - Directory name; default to $SPINSPECT_DIR/data/*
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Issue warnings on the following:
;     * File name not matchings spInspect-XXXX-YYYYY.par, where XXXX and YYYYY
;       are integers
;     * Error reading the file as a Yanny-formatted file
;     * Empty data file
;     * BOSSOBJECT structure missing
;     * Any structure name that does not match a structure name in the template
;       file
;     * For any structure that does exist, if any of the elements in the
;       template are missing (those are all considered required)
;     * For any structure that does exist, if any of the elements are a
;       different variable type from the template file
;     * BOSSOBJECT.FIBERID not in the range [1,1000]
;     * BOSSOBJECT.FIBERID has duplicate values
;     * BOSSOBJECT.CLASS_PERSON not in the range [0,4]
;     * BOSSOBJECT.Z_CONF_PERSON not in the range [0,4]
;   The template file is $IDLSPEC2D_DIR/opfiles/spInspect-0000-00000.par
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/opfiles/spInspect-0000-00000.par
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   23-Jul-2010  Written by David Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro inspect_verify1, filename

   common com_inspect_verify, tstruct, tstnames

   if (NOT keyword_set(tstruct)) then begin
      tfile = djs_filepath('spInspect-0000-00000.par', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdir='opfiles')
      yanny_read, tfile, tstructs, stnames=tstnames
   endif

   splog, prelog=filename
   shortname = fileandpath(filename)
   if (strmatch(shortname,'spInspect-[0-9][0-9][0-9][0-9]-[0-9][0-9][0-9][0-9][0-9].par') EQ 0) then $
    splog, 'WARNING: Filename not standard'
   thisplate = strmid(shortname,10,4)
   thismjd = strmid(shortname,15,5)

   yanny_read, filename, pdata, hdr=hdr, stnames=stnames, $
    /anonymous, errcode=errcode
   if (errcode NE 0) then $
    splog,' WARNING: Error reading Yanny file '+string(errcode)
   if (keyword_set(pdata) EQ 0) then $
    splog, 'WARNING: Empty file'
   i = (where(stnames EQ 'BOSSOBJECT', ct))[0]
   if (ct EQ 0) then begin
      splog, 'WARNING: Missing BOSSOBJECT structure'
   endif else begin
      thisdat = *(pdata[i])
      if (tag_exist(thisdat,'PLATE')) then begin
         indx = where(thisdat.plate NE thisplate, nbad)
         if (nbad GT 0) then $
          splog, 'WARNING: PLATE values do not match filename'
      endif
      if (tag_exist(thisdat,'MJD')) then begin
         indx = where(thisdat.mjd NE thismjd, nbad)
         if (nbad GT 0) then $
          splog, 'WARNING: MJD values do not match filename'
      endif
      if (tag_exist(thisdat,'FIBERID')) then begin
         indx = where(thisdat.fiberid LT 1 OR thisdat.fiberid GT 1000, nbad)
         if (nbad GT 0) then splog, 'WARNING: FIBERID out of range'
         fid = thisdat.fiberid
         if (n_elements(uniq(fid, sort(fid))) NE n_elements(fid)) then $
          splog, 'WARNING: Duplicate values of FIBERID'
      endif
      if (tag_exist(thisdat,'CLASS_PERSON')) then begin
         indx = where(thisdat.class_person LT 0 $
          OR thisdat.class_person GT 4, nbad)
         if (nbad GT 0) then splog, 'WARNING: CLASS_PERSON out of range'
      endif
      if (tag_exist(thisdat,'Z_CONF_PERSON')) then begin
         indx = where(thisdat.z_conf_person LT 0 $
          OR thisdat.z_conf_person GT 4, nbad)
         if (nbad GT 0) then splog, 'WARNING: Z_CONF_PERSON out of range'
      endif
   endelse

   ; For any structure in the file, it must contain all of the required
   ; elements in the template file, and be of the same data type
   if (keyword_set(pdata)) then begin
      for i=0, n_elements(stnames)-1 do begin
         j = (where(stnames[i] EQ tstnames, ct1))[0]
         if (ct1 EQ 0) then begin
            splog, 'Unknown structure ' + stnames[i]
         endif else begin
            tags1 = tag_names(*pdata[i])
            tags2 = tag_names(*(tstructs[j]))
            for k=0, n_elements(tags2)-1 do begin
               l = (where(tags1 EQ tags2[k], ct2))[0]
               if (ct2 EQ 0) then begin
                  splog, 'Missing variable '+tags2[k]
               endif else begin
                  if (size((*(pdata[i])).(l), /tname) $
                   NE size((*(tstructs[j])).(k), /tname)) then $
                     splog, 'Wrong data type for '+tags2[k]
               endelse
            endfor
         endelse
      endfor
   endif

   return
end

;------------------------------------------------------------------------------
pro inspect_verify, path=path1

   if (keyword_set(path1)) then path = path1 $
    else path = getenv('SPINSPECT_DIR')+'/data/*'

   alldir = file_search(path, /test_directory, count=ndir)
   if (ndir EQ 0) then begin
      print, 'No directories found'
      return
   endif

   for idir=0L, ndir-1L do begin
      files = file_search(djs_filepath('*', root_dir=alldir[idir]), count=nfile)
      for ifile=0L, nfile-1L do begin
         inspect_verify1, files[ifile]
      endfor
   endfor

   return
end
;------------------------------------------------------------------------------
