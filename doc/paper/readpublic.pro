pro readpublic, public=public1, columns=columns, zans=zans, plug=plug

   q_zans = arg_present(zans)
   q_plug = arg_present(plug)

   if (keyword_set(public1)) then public = public1 $
    else public = ['EDR','DR1','DR2','DR3','DR4','DR5']
   platelist, plist=plist
   qkeep = bytarr(n_elements(plist))
   for i=0, n_elements(public)-1 do $
    qkeep = qkeep OR strmatch(plist.public,public[i])
   qkeep = qkeep AND strmatch(plist.statuscombine,'Done*')
   plist = plist[where(qkeep, nplate)]

   for iplate=0L, nplate-1L do begin
      print, format='("Plate ",i4," of ",i4,a1,$)', iplate, nplate, string(13b)
      if (q_zans AND q_plug) then $
       readspec, plist[iplate].plate, mjd=plist[iplate].mjd, $
        zans=zans1, plug=plug1, /silent $
      else if (q_zans) then $
       readspec, plist[iplate].plate, mjd=plist[iplate].mjd, $
        zans=zans1, /silent $
      else if (q_plug) then $
       readspec, plist[iplate].plate, mjd=plist[iplate].mjd, $
        plug=plug1, /silent
      if (keyword_set(columns)) then begin
         zans1 = struct_selecttags(zans1, select_tags=columns)
         plug1 = struct_selecttags(plug1, select_tags=columns)
      endif
      if (iplate EQ 0 AND keyword_set(zans1)) then $
       zans = replicate(zans1[0], nplate*640L)
      if (iplate EQ 0 AND keyword_set(plug1)) then $
       plug = replicate(plug1[0], nplate*640L)
      if (keyword_set(zans1)) then $
       copy_struct_inx, zans1, zans, index_to=640L*iplate+lindgen(640)
      if (keyword_set(plug1)) then $
       copy_struct_inx, plug1, plug, index_to=640L*iplate+lindgen(640)
   endfor
   print

   return
end
