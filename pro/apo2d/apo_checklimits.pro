;------------------------------------------------------------------------------
function apo_checklimits, flavor, field, camera, value, html=html

   common apo_limits, slimits

   markstring = ''
   if (NOT keyword_set(value)) then return, markstring

   ; Read this Yanny file only the first time this routine is called,
   ; then save the limits in a common block.
   if (NOT keyword_set(slimits)) then begin
      limitfile = filepath('opLimits.par', root_dir=getenv('IDLSPEC2D_DIR'), $
       subdirectory='examples')
      yanny_read, limitfile, pdata
      slimits = *pdata[0]
      yanny_free, pdata
   endif

   indx = where(slimits.field EQ field AND slimits.camera EQ camera $
            AND slimits.flavor EQ flavor, nlim)
   for ilim=0, nlim-1 do begin
      if (value GE slimits[indx[ilim]].lovalue $
       AND value LE slimits[indx[ilim]].hivalue) then begin
          markstring = slimits[indx[ilim]].color
          if (keyword_set(html)) then $
           markstring = '<B><FONT COLOR="' $
            + apo_color2hex(markstring) + '">'
      endif
   endfor

   return, markstring
end
;------------------------------------------------------------------------------
