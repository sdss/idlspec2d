;+
; NAME:
;   sphdrfix
;
; PURPOSE:
;   Fix header cards in raw SDSS spectroscopic data.
;
; CALLING SEQUENCE:
;   sphdrfix, filename, hdr
;
; INPUTS:
;   filename   - Name of raw FITS file
;   hdr        - FITS header
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   hdr        - FITS header (modified)
;
; COMMENTS:
;   This routine implements "hand edits" of the raw FITS headers for
;   SDSS spectroscopic images.  The list of edits to make can be stored
;   in two possible Yanny parameter file, the static "opHdrFix.par" file
;   in the IDLSPEC2D product, and in the "sdHdrFix-$MJD.par" file for
;   the relevant night.
;
;   This proc only works in IDL version 5.3 and later, because it
;   uses STRMATCH().
;
; EXAMPLES:
;   filename = 'sdR-b2-00003976.fit'
;   hdr = headfits(filename)
;   sphdrfix, filename, hdr
;
; BUGS:
;   This will fail if the MJD needs hand-editing in the sdHdrFix file,
;   since we need the MJD to find the sdHdrFix file in the first place!
;
; PROCEDURES CALLED:
;   fileandpath()
;   splog
;   sxaddpar
;   yanny_free
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES
;   sphdrfix1
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/opHdrFix.par
;   $SPECLOG_DIR/$MJD/sdHdrFix-$MJD.par
;
; REVISION HISTORY:
;   29-Dec-2009  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro sphdrfix1, filename, hdr, hfixpar

   npar = n_elements(hfixpar)

   ;----------
   ; Find any fixes that are relevant for this file

   fileroot = fileandpath(filename)
   i = strpos(fileroot, '.')
   if (i NE -1) then fileroot = strmid(fileroot, 0, i)
   qmatch = bytarr(npar)
   for ipar=0, npar-1 do $
    qmatch[ipar] = strmatch(fileroot, hfixpar[ipar].fileroot)
   indx = where(qmatch, nmatch)

   ;----------
   ; Loop over each relevant fix and apply it

   hdrcards = strmid(hdr, 0, 8)

   for im=0, nmatch-1 do begin
      thiskey = hfixpar[indx[im]].keyword

      thisvalue = hfixpar[indx[im]].value

      ; This value is a string if it contains a single quote.
      qstring = strmatch(thisvalue, "*'*")
      if (qstring) then begin
         ; Extract the first string between single quotes.
         thisvalue = strmid(thisvalue, strpos(thisvalue, "'") )
;         thisvalue = (strsplit(thisvalue, "'", /extract))[0]
         squote = "'"
         thisvalue = (strsplit(thisvalue, squote, extract=1))[0]
;         thisvalue = (strsplit(thisvalue, "'"))[0]
      endif else begin
         ; This value is either a floating-point or integer value.
         ; Floating-point if contains any of the characters '.de'
         qfloat = strmatch(thisvalue, '*[.de]*')
         fval = double(thisvalue)
         ival = long(thisvalue)
         if (ival NE fval OR qfloat) then thisvalue = fval $
          else thisvalue = ival
      endelse

      splog, 'Changing ' + thiskey + '=' + string(sxpar(hdr,thiskey)) $
       + ' -> ' + string(thisvalue)
      sxaddpar, hdr, thiskey, thisvalue
   endfor

   return
end

;------------------------------------------------------------------------------
pro sphdrfix, filename, hdr

   common spec2d_hfixpar, hfix1

   if (n_params() LT 2) then begin
      print, 'Syntax - sphdrfix, filename, hdr'
      return
   endif

   ;----------
   ; Read the "opHdrFix.par" Yanny file only the first time this routine
   ; is called, then save the values in a common block.

   if (NOT keyword_set(hfix1)) then begin
      parfile = filepath('opHdrFix.par', root_dir=getenv('IDLSPEC2D_DIR'), $
       subdirectory='examples')
      yanny_read, parfile, pdata
      hfix1 = *pdata[0]
      yanny_free, pdata
   endif

   if (NOT keyword_set(hfix1)) then $
    message, 'Parameter file opHdrfix.par not found!'

   ; Apply header fixes from  the "opHdrFix.par" file.
   sphdrfix1, filename, hdr, hfix1

   ;----------
   ; Read the sdHdrFix file for this night to look for more possible header
   ; changes from the APO observers.

   speclog_dir = getenv('SPECLOG_DIR')
   if (NOT keyword_set(speclog_dir)) then $
    message, 'Must set environment variable SPECLOG_DIR'
   mjd = sxpar(hdr, 'MJD')
   mjdstr = string(mjd, format='(i05.5)')
   plugdir = concat_dir(speclog_dir, mjdstr)

   reportfile = filepath('sdHdrFix-'+mjdstr+'.par', root_dir=plugdir)
   yanny_read, reportfile, pdata, stnames=stnames, /anonymous, errcode=errcode

   if (keyword_set(pdata)) then begin
      ; Find the ONEEXP structure
      for i=0, N_elements(pdata)-1 do $
;       if (tag_names(*pdata[i], /structure_name) EQ 'OPHDRFIX') then $
       if (stnames[i] EQ 'OPHDRFIX') then $
        hfix2 = *pdata[i]

      ; Apply header fixes from  the "sdHdrFix-$MJD.par" file.
      if (keyword_set(hfix2)) then $
       sphdrfix1, filename, hdr, hfix2

      yanny_free, pdata
   endif else if (errcode NE 0) then begin
      splog, 'WARNING: Error reading sdHdrFix file '+fileandpath(reportfile)
   endif

   return
end
;------------------------------------------------------------------------------
