;------------------------------------------------------------------------------
pro matchplugmap, plug1, plug2, indx1, indx2, mindist=mindist

   if (NOT keyword_set(mindist)) then mindist = 1.0 / 3600

   nfiber1 = n_elements(plug1)
   nfiber2 = n_elements(plug2)
   indx1 = lindgen(nfiber1)
   indx2 = lonarr(nfiber1) - 1L

   for ifiber=0, nfiber1-1 do begin
      adist = djs_diff_angle(plug1[ifiber].ra, plug1[ifiber].dec, $
       plug2.ra, plug2.dec)
      dmin = min(adist, imin)
      if (dmin LE mindist) then indx2[ifiber] = imin
   endfor

   igood = where(indx2 NE -1L)
   if (igood[0] NE -1) then begin
      indx1 = indx1[igood]
      indx2 = indx2[igood]
   endif else begin
      indx1 = -1L
      indx2 = -1L
   endelse

   return
end
;------------------------------------------------------------------------------
