
;
;
;   Read in QSO from a list of plate, mjds
; IDL> getqso, 410, 51877, zmin=0.3, zmax=2.2, imax=20.0, $
;      flux=flux, invvar=invvar, mask=mask, $
;      disp=disp, plugmap=plugmap, loglam=loglam, zans=zans, model=model
;

pro getqso, plate, mjd,  zmin=zmin, zmax=zmax, imax=imax, dangerous=dangerous, $
      flux=flux, invvar=invvar, mask=mask, mgspec=mgspec, $
      disp=disp, plugmap=plugmap, loglam=loglam, zans=zans, model=model


    if NOT keyword_set(zmin) then zmin=0.0
    if NOT keyword_set(zmax) then zmax=10.0
    if NOT keyword_set(imax) then imax=30.0
    if NOT keyword_set(dangerous) then dangerous = 0

    nplate = n_elements(plate)
    nmjd = n_elements(mjd)

    if nmjd EQ 0 then mjd = lonarr(nplate)
    nmjd = n_elements(mjd)
    if nmjd NE nplate then begin
      message, 'plate does not match mjd'
      return
    endif

    for i=0,nplate - 1 do begin

       readspec, plate[i], mjd=mjd[i], zans=zans, plug=plug, /silent
 
       if keyword_set(zans) then begin
         good = where(strtrim(zans.class,2) EQ 'QSO' AND $
                    (dangerous OR zans.zwarning EQ 0) AND $
                    zans.z GT zmin AND zans.z LT zmax AND $
                    plug.mag[3] LE imax ,ng)


         if good[0] NE -1 then begin
            zg =  zans[good]

            if NOT keyword_set(fullplate) then fullplate=zg.plate $
            else fullplate = [fullplate, zg.plate]

            if NOT keyword_set(fullmjd) then fullmjd =zg.mjd $
            else fullmjd = [fullmjd, zg.mjd]

            if NOT keyword_set(fullfiber) then fullfiber=zg.fiberid $
            else fullfiber= [fullfiber, zg.fiberid]
          endif
        endif
    endfor

    if NOT keyword_set(fullplate) then return

    readspec, fullplate, fullfiber, mjd=fullmjd, $
      flux=flux, invvar=invvar, andmask=mask, $
      disp=disp, plugmap=plugmap, loglam=loglam, zans=zans, /silent

    mgspec = qso_createspec(loglam, flux, invvar, plugmap, zans)

    return
end


                    
      	
