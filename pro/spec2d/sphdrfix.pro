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
;   SDSS spectroscopic images.  The list of edits to make are stored
;   in a Yanny parameter file.
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
;
; PROCEDURES CALLED:
;   fileandpath()
;   splog
;   sxaddpar
;   yanny_free
;   yanny_read
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/opHdrFix.par
;
; REVISION HISTORY:
;   29-Dec-2009  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro sphdrfix, filename, hdr

   common spec2d_hfixpar, hfixpar

   if (n_params() LT 2) then begin
      print, 'Syntax - sphdrfix, filename, hdr'
      return
   endif

   ;----------
   ; Read this Yanny file only the first time this routine is called,
   ; then save the values in a common block.

   if (NOT keyword_set(hfixpar)) then begin
      parfile = filepath('opHdrFix.par', root_dir=getenv('IDLSPEC2D_DIR'), $
       subdirectory='examples')
      yanny_read, parfile, pdata
      hfixpar = *pdata[0]
      yanny_free, pdata
   endif
   npar = n_elements(hfixpar)

   if (NOT keyword_set(hfixpar)) then $
    message, 'Parameter file not found!'

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
