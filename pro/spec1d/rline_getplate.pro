
; PRIMTARGET - If set, then select only objects where at least one bit
;              in its PRIMTARGET flag matches one bit in this parameter.
; CLASS      - If set, then only select objects that match this classification,
;              i.e. 'GALAXY'
; MAXRCHI2   - If set, then only select objects whose reduced chi^2 is
;              less than or equal to this value.
; QUICK      - If set, then only return ZANS and PLUG in the return structure.
;------------------------------------------------------------------------------
function rline_getplate, plate, mjd=mjd, $
 primtarget=primtarget, class=class, maxrchi2=maxrchi2, quick=quick

     ; First read the plug-map

     readspec, plate, mjd=mjd, plug=plug, zans=zans
     if (NOT keyword_set(plug)) then return, 0
 
     ; Set MASK=1 for objects to return

     mask = lonarr(n_elements(plug)) + 1

     if (keyword_set(primtarget)) then $
      mask = mask * (plug.primtarget NE 0) $
       * ((plug.primtarget AND primtarget) NE 0)

     if (keyword_set(class)) then $
      mask = mask * (strtrim(zans.class,2) EQ class)

     if (keyword_set(maxrchi2)) then $
      mask = mask * (zans.rchi2 LE maxrchi2)

     igood = where(mask)
     if (igood[0] EQ -1) then return, 0

     if (keyword_set(quick)) then begin
        tt = { plug : plug[0], $
               zans : zans[0] }
        tt = replicate(tt, n_elements(igood))
        tt.zans = zans[igood]
        tt.plug = plug[igood]
        return, tt
     endif

     readspec, zans[igood].plate, zans[igood].fiberid, mjd=zans[igood].mjd, $
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

     tt = replicate(tt, n_elements(zans))
     tt.flux = temporary(flux)
     tt.finv = temporary(finv)
     tt.loglam = temporary(loglam)
     tt.model  = temporary(model)
     tt.plug   = temporary(plug)
     tt.zans   = temporary(zans)

     return, tt
end 
;------------------------------------------------------------------------------
