; Find the reddening for a bunch of plates.
pro redplate

   platelist, plist=plist
   plist = plist[where(plist.qsurvey, nplate)]

   ebv = fltarr(nplate)
   for iplate=0, nplate-1 do begin
print, 'Plate ', iplate, ' of ', nplate
      readspec, plist[iplate].plate, mjd=plist[iplate].mjd, plug=plug
      euler, plug.ra, plug.dec, ll, bb, 1
      ebv[iplate] = mean(dust_getval(ll,bb,/interp, /noloop))
   endfor

stop
end
