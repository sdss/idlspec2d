;------------------------------------------------------------------------------
pro spframe_read, filename, indx, objflux=objflux, objivar=objivar, $
 mask=mask, wset=wset, loglam=loglam, dispset=dispset, dispimg=dispimg, $
 plugmap=plugmap, skyflux=skyflux, hdr=hdr, adderr=adderr

   qtrim = n_elements(indx) GT 0

   if (NOT keyword_set(filename)) then $
    message, 'Must specify FILENAME'

   if (arg_present(hdr)) then hdr = headfits(filename[0])
   if (arg_present(objflux) $
    OR (arg_present(objivar) AND keyword_set(adderr))) then begin
      objflux = mrdfits(filename[0], 0, /silent)
      if (qtrim) then objflux = objflux[*,indx]
   endif
   if (arg_present(objivar)) then begin
      objivar = mrdfits(filename[0], 1, /silent)
      if (qtrim) then objivar = objivar[*,indx]
      if (keyword_set(adderr)) then begin
         gmask = objivar NE 0 ; =1 for good points
         objivar = 1.0 / ( 1.0/(objivar + (1-gmask)) $
          + (adderr * (objflux>0))^2 ) * gmask
      endif
   endif
   if (arg_present(mask)) then begin
      mask = mrdfits(filename[0], 2, /silent)
      if (qtrim) then mask = mask[*,indx]
   endif
   if (arg_present(wset) OR arg_present(loglam)) then begin
      wset = mrdfits(filename[0], 3, /silent)
      if (qtrim) then wset = traceset_trim(wset, indx)
      if (arg_present(loglam)) then traceset2xy, wset, xtmp, loglam
      xtmp = 0
   endif
   if (arg_present(dispset) OR arg_present(dispimg)) then begin
      dispset = mrdfits(filename[0], 4, /silent)
      if (qtrim) then dispset = traceset_trim(dispset, indx)
      if (arg_present(dispimg)) then traceset2xy, dispset, xtmp, dispimg
      xtmp = 0
   endif
   if (arg_present(plugmap)) then begin
      plugmap = mrdfits(filename[0], 5, /silent)
      if (qtrim) then plugmap = plugmap[indx]
   endif
   if (arg_present(skyflux)) then begin
      skyflux = mrdfits(filename[0], 6, /silent)
      if (qtrim) then skyflux = skyflux[*,indx]
   endif

   return
end
;------------------------------------------------------------------------------
