;------------------------------------------------------------------------------
pro ztweak_star

   snmax = 100

   ;----------
   ; Read the input spectra

   filename = filepath('eigeninput_star.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
   yanny_read, filename, pdat
   slist = *pdat[0]
   yanny_free, pdat

; ???
;ii=where(strtrim(slist.class,2) EQ 'K')
;slist=slist[ii]
;slist=slist[0:9]
;slist.plate = 406
;slist.mjd = 51869
;slist.fiberid = [9,18,21,34,39,53,62,64,80,102]
   readspec, slist.plate, slist.fiberid, mjd=slist.mjd, $
    flux=objflux, invvar=objivar, $
    andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=objloglam, /align
   nobj = n_elements(slist)

   ;----------
   ; Insist that all of the requested spectra exist

   imissing = where(plugmap.fiberid EQ 0, nmissing)
   if (nmissing GT 0) then begin
      for i=0, nmissing-1 do $
       print, 'Missing plate=', slist[imissing[i]].plate, $
        ' mjd=', slist[imissing[i]].mjd, $
        ' fiber=', slist[imissing[i]].fiberid
      message, string(nmissing) + ' missing object(s)'
   endif

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-sub residuals.

   objivar = skymask(objivar, andmask, ormask)
andmask = 0 ; Free memory
ormask = 0 ; Free memory

   if (keyword_set(snmax)) then begin
      ifix = where(objflux^2 * objivar GT snmax^2)
      if (ifix[0] NE -1) then objivar[ifix] = (snmax/objflux[ifix])^2
   endif

   ;----------
   ; Find the best-fit Elodie star for each spectrum

   objdloglam = objloglam[1] - objloglam[0]
   res = elodie_best(objflux, objivar, $
    objloglam0=objloglam[0], objdloglam=objdloglam)
stop


   return
end
;------------------------------------------------------------------------------
