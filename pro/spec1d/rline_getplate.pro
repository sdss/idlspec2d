
function rline_getplate, plate, primtarget=primtarget, noz=noz, rchi2=rchi2

     ; first get plugmap

     readspec, [plate], plug=plug, zans=zans

     if (size(plug,/tname) EQ 'INT') then return,0
    
     mask = lonarr(n_elements(plug)) 
     if keyword_set(primtarget ) then $
        mask = mask + ((plug.primtarget AND primtarget) NE 0)

     ; primtarget was not set
     if total(mask) EQ 0 then mask = mask + 1

     if NOT keyword_set(noz) then begin
        bad = where(zans.z EQ 0)
        if bad[0] NE -1 then mask[bad] = 0
     endif

     if keyword_set(rchi2) then begin
        bad = where(zans.rchi2 GT rchi2)
        if bad[0] NE -1 then mask[bad] = 0
     endif

     good = where(mask)

     if good[0] EQ -1 then return, 0

     readspec, zans[good].plate, zans[good].fiberid, $
         flux=flux, invvar=finv, loglam=loglam, plug=plug, zans=zans

     model = loglam * 0
     for i=0,n_elements(zans)-1 do $
       model[*,i] = synthspec(zans[i], loglam=loglam[*,i]) 

     tt = { flux : flux[*,0],      $
            finv : finv[*,0],      $
            loglam : float(loglam[*,0]) , $
            model : float(model[*,0]), $
            plug : plug[0], $
            zans : zans[0] }

     yy = replicate(tt, n_elements(zans))
     yy.flux = temporary(flux)
     yy.finv = temporary(finv)
     yy.loglam = temporary(loglam)
     yy.model  = temporary(model)
     yy.plug   = temporary(plug)
     yy.zans   = temporary(zans)
    
     return, yy
end 

