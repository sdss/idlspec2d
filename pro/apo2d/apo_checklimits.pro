;------------------------------------------------------------------------------
function apo_checklimits, flavor, field, camera, value, html=html

   common apo_limits, numlimits, textlimits

   markstring = ''
   if (NOT keyword_set(value)) then return, markstring

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
          AND strmatch(value, textlimits[ilim].strval)) then begin
            markstring = textlimits[ilim].color
            if (keyword_set(html)) then $
             markstring = '<B><FONT COLOR="' $
              + apo_color2hex(markstring) + '">'
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
             markstring = '<B><FONT COLOR="' $
              + apo_color2hex(markstring) + '">'
         endif
      endfor
   endelse

   return, markstring
end
;------------------------------------------------------------------------------
