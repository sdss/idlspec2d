
function qso_createspec, loglam, flux, invvar,  plug, zans, $
       restrange = restrange

;
;  This will include up to CaII lines?!?
;
    if NOT keyword_set(restrange) then restrange=[1220.0,4000.0]

    nqso = n_elements(zans)
    tt = { flux : flux[*,0], finv : invvar[*,0], loglam : float(loglam[*,0]) , $
            model : flux[*,0]*0.0, plug : plug[0], zans : zans[0], $
            cont : flux[*,0]*0.0, mask : long(flux[*,0]*0.0)  }

    spec = replicate(tt,nqso)
    spec.flux = 0.0
    spec.finv = 0.0
    spec.loglam = loglam
    spec.plug = plug
    spec.zans = zans
  
    for i=0, nqso - 1 do begin

      range = alog10((1+zans[i].z)*restrange)
      print, zans[i].plate, zans[i].fiberid, format='(i,i)'
      qso_slowcont, loglam[*,i], flux[*,i], invvar[*,i], $
           range=range[*], model=model, mask=mask

      good = where(model GT 0.0)    

      if good[0] NE -1 then begin
        spec[i].flux[good] = 1.0 - flux[good,i]/model[good]
        spec[i].finv[good] = invvar[good,i]*model[good]^2
        spec[i].mask = mask
        spec[i].cont = model
      endif
    endfor

return, spec 
end

