;+
; NAME:
;   plugmap_verify
;
; PURPOSE:
;   Sanity checks for BOSS plPlugMapP files
;
; CALLING SEQUENCE:
;   plugmap_verify, files
;
; INPUTS:
;   files     - File name(s), which can contain wildcards
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Do a number of sanity checks on the format and objects in
;   plPlugMapP files for BOSS.
;
; EXAMPLES:
;
; BUGS:
;   No tests are performed on whether the objects are good SDSS targets.
;
; PROCEDURES CALLED:
;   djs_filepath()
;   splog
;   yanny_par()
;   yanny_read
;   yanny_readone()
;
; INTERNAL PROCEDURES:
;   plugmap_warn
;   plugmap_verify1
;
; REVISION HISTORY:
;   09-Oct-2009  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro plugmap_warn, output
   splog, 'WARNING: ', output, /noname
end
;------------------------------------------------------------------------------
pro plugmap_verify1, filename

;   splog, 'Verifying '+filename
   splog, prelog=fileandpath(filename)
   thisplate = long(strmid(fileandpath(filename),11,4))

   ; Read the file and check the structure names
   yanny_read, filename, pp, hdr=hdr, stnames=stnames, /anonymous, $
    errcode=errcode
   if (errcode NE 0) then $
    plugmap_warn, 'Yanny format invalid'
   if (n_elements(pp) NE 1) then $
    plugmap_warn, 'Too many structures '+stnames
   if (stnames[0] NE 'PLUGMAPOBJ') then $
    plugmap_warn, 'Structure name is not PLUGMAPOBJ'
   plugmap = *pp[0]
   yanny_free, pp
   nkey = n_elements(hdr)
   allkey = strarr(nkey)
   for i=0, nkey-1 do allkey[i] = (str_sep(hdr[i],' '))[0]

   ; Test the header
   template_file = djs_filepath('plPlugMapP-verify.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdir='examples')
   junk = yanny_readone(template_file, hdr=thdr, /anonymous)
   if (NOT keyword_set(thdr)) then $
    message, 'Unable to read plPlugMapP template file'
   thdr = thdr[where(strtrim(thdr))] ; Get rid of blank lines
   for i=0, n_elements(thdr)-1 do begin
      thiskey = (str_sep(thdr[i],' '))[0]
      j = where(thiskey EQ allkey, ct)
      if (ct EQ 0) then begin
         plugmap_warn, 'Missing header card '+thiskey
      endif else if (ct GT 1) then begin
         plugmap_warn, 'Multiple copies of header card '+thiskey
      endif else begin
         val1 = yanny_par(thiskey, thdr)
         val2 = yanny_par(thiskey, hdr)
         if (n_elements(val1) NE n_elements(val2)) then begin
            plugmap_warn, 'Incorrect number of values for header card '+thiskey
         endif else begin
            nbad = total(val1 NE val2)
            if (nbad GT 0) then $
             plugmap_warn, 'Incorrect values for header card '+thiskey
         endelse
      endelse
   endfor
   plateid = yanny_par(hdr,'plateId')
   if (n_elements(plateid) NE 1) then $
    plugmap_warn, 'Incorrect number of values for header card plateId'
   if (n_elements(plateid) GE 1) then $
    if (plateid[0] NE thisplate) then $
     plugmap_warn, 'Header plateId inconsistent with file name'
   redden = yanny_par(hdr,'reddeningMed')
   if (n_elements(redden) NE 5) then $
    plugmap_warn, 'Wrong number of elements in reddeningMed'

   if (n_elements(plugmap) EQ 0) then begin
      plugmap_warn, 'No entries in plugmap!'
      splog, prelog=''
      return
   endif

   ; Check the objects
   qobj = strmatch(plugmap.holetype,'OBJECT*')
   iobj = where(qobj, nobj)
   if (nobj NE 1000) then $
    plugmap_warn, 'Wrong number of object fibers '+strtrim(nobj,2)

   qsky = qobj AND strmatch(plugmap.objtype,'SKY*')
   isky = where(qsky, nsky)
   if (nsky LT 80) then $
    plugmap_warn, 'Too few skies '+strtrim(nsky,2)
   if (nsky GT 120) then $
    plugmap_warn, 'Too many skies '+strtrim(nsky,2)

   isky1 = where(qsky AND plugmap.yfocal GE 0, nsky1)
   isky2 = where(qsky AND plugmap.yfocal LT 0, nsky2)
   if (nsky1 LT 8) then $
    plugmap_warn, 'Only '+strtrim(nsky1,2)+' skies on North half'
   if (nsky2 LT 8) then $
    plugmap_warn, 'Only '+strtrim(nsky2,2)+' skies on South half'

   qstd = qobj AND strmatch(plugmap.objtype,'SPECTROPHOTO_STD*')
   istd = where(qstd, nstd)
   if (nstd LT 20) then $
    plugmap_warn, 'Too few std stars '+strtrim(nstd,2)
   if (nstd GT 30) then $
    plugmap_warn, 'Too many std stars '+strtrim(nstd,2)

   istd1 = where(qstd AND plugmap.yfocal GE 0, nstd1)
   istd2 = where(qstd AND plugmap.yfocal LT 0, nstd2)
   if (nstd1 LT 8) then $
    plugmap_warn, 'Only '+strtrim(nstd1,2)+' std stars on North half'
   if (nstd2 LT 8) then $
    plugmap_warn, 'Only '+strtrim(nstd2,2)+' std stars on South half'

   ; Check for the other hole types
   qqual = strmatch(plugmap.holetype,'QUALITY*')
   iqual = where(qqual, nqual)
   if (nqual NE 1) then $
    plugmap_warn, 'Expected 1 QUALITY (center) hole, got '+strtrim(nqual,2)
   icenter = where(qqual AND plugmap.xfocal EQ 0 AND plugmap.yfocal EQ 0, $
    ncenter)
   if (ncenter EQ 0) then $
    plugmap_warn, 'No QUALITY hole in center'

   qguide = strmatch(plugmap.holetype,'GUIDE*')
   iguide = where(qguide, nguide)
   if (nguide NE 16) then $
    plugmap_warn, 'Expected 16 GUIDE stars, got '+strtrim(nguide,2)

   iguide1 = where(qguide AND plugmap.yfocal GE 0, nguide1)
   iguide2 = where(qguide AND plugmap.yfocal LT 0, nguide2)
   if (nguide1 LT 6) then $
    plugmap_warn, 'Only '+strtrim(nguide1,2)+' GUIDE stars on North half'
   if (nguide2 LT 6) then $
    plugmap_warn, 'Only '+strtrim(nstd2,2)+' GUIDE stars on South half'

   qalign = strmatch(plugmap.holetype,'ALIGNMENT*')
   ialign = where(qalign, nalign)
   if (nalign NE 16) then $
    plugmap_warn, 'Expected 16 ALIGNMENT holes, got '+strtrim(nalign,2)

   splog, prelog=''
   return
end
;------------------------------------------------------------------------------
pro plugmap_verify, filename

   allfile = findfile(filename, count=nfile)
   splog, nfile, ' files found'
   if (nfile EQ 0) then return

   for ifile=0L, nfile-1L do $
    plugmap_verify1, allfile[ifile]

   return
end
;------------------------------------------------------------------------------
