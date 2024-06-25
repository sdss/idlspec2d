;+
; NAME:
;   apo_checklimits()
;
; PURPOSE:
;   Convert output FITS file from APOREDUCE to HTML format.
;
; CALLING SEQUENCE:
;   markstring = apo_checklimits(flavor, field, camera, value, [ /html ] )
;
; INPUTS:
;   flavor     - FLAVOR to match in the opLimits file.
;   field      - FIELD to match in the opLimits file.
;   camera     - CAMERA to match in the opLimits file.
;   value      - Value to test in the opLimits file.  If this is a
;                string value, then one matches to STRVAL in that file.
;                Otherwise, test for values within [LOVALUE,HIVALUE].
;
; OPTIONAL INPUTS:
;   html       - If set, then convert the color name in MARKSTRING into
;                an HTML string.  For example, a return value of 'red'
;                becomes '<B><FONT COLOR="#FF0000">'.
;
; OUTPUT:
;   markstring - Return the COLOR from the opLimits file.
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
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/opLimits.par
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
function apo_checklimits, flavor, field, camera, value, html=html

   common apo_limits, numlimits, textlimits

   markstring = ''
   if (n_elements(value) EQ 0) then return, markstring

   ;----------
   ; Read this Yanny file only the first time this routine is called,
   ; then save the limits in a common block.

   if (NOT keyword_set(numlimits)) then begin
      limitfile = filepath('opLimits.par', root_dir=getenv('IDLSPEC2D_DIR'), $
       subdirectory='examples')
      yanny_read, limitfile, pdata, stnames=stnames
      numlimits = *pdata[(where(stnames EQ 'SPECLIMIT'))[0]]
      textlimits = *pdata[(where(stnames EQ 'TEXTLIMIT'))[0]]
      yanny_free, pdata
   endif

   if (size(value,/tname) EQ 'STRING') then begin
      ;----------
      ; Case of text limit

      for ilim=0, n_elements(textlimits)-1 do begin
         if (strmatch(field, textlimits[ilim].field) $
          AND strmatch(camera, textlimits[ilim].camera) $
          AND strmatch(flavor, textlimits[ilim].flavor) $
          AND strmatch(strtrim(value,2), textlimits[ilim].strval)) then begin
            markstring = textlimits[ilim].color
            if (keyword_set(html)) then $
             markstring = '<span style="color:' $
              + apo_color2hex(markstring) + ';font-weight:bold;">'
         endif
      endfor
   endif else begin
      ;----------
      ; Case of floating-point limit

      for ilim=0, n_elements(numlimits)-1 do begin
         if (strmatch(field, numlimits[ilim].field) $
          AND strmatch(camera, numlimits[ilim].camera) $
          AND strmatch(flavor, numlimits[ilim].flavor) $
          AND value GE numlimits[ilim].lovalue $
          AND value LE numlimits[ilim].hivalue) then begin
            markstring = numlimits[ilim].color
            if (keyword_set(html)) then $
             markstring = '<span style="color:' $
              + apo_color2hex(markstring) + ';font-weight:bold;">'
         endif
      endfor
   endelse

   return, markstring
end
;------------------------------------------------------------------------------
