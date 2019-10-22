; Find the reddening for a bunch of plates.
pro redplate, legacy=legacy, plates=plates

   ;platelist, plist=plist
   conflist, plist=plist, legacy=legacy, plates=plates
   ;print, plist
   plist = plist[where(plist.qsurvey, nplate)]
   ;print, plist
   ebv = fltarr(nplate)
   for iplate=0, nplate-1 do begin
      readspec, plist[iplate].field, mjd=plist[iplate].mjd, plug=plug, $
        legacy=legacy, plates=plates
      euler, plug.ra, plug.dec, ll, bb, 1
      ebv[iplate] = mean(dust_getval(ll,bb,/interp, /noloop))
      print, 'Field ', iplate, ' of ', nplate,' ebv ',ebv[iplate]
   endfor

stop
end
