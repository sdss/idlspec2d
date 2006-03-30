
;------------------------------------------------------------------------------
pro plot_thru, plate, mjd=mjd

   platelist, plist=plist
   indx = where(keyword_set(strtrim(plist.public)) $
AND plist.plate GT 400 AND plist.plate LT 410 $ ; ???
    AND strmatch(plist.statuscombine,'Done*'), nplate)
   plist = plist[indx]

   ;----------
   ; Start by counting the number of CCD frames

   print, 'Counting exposures'
   ntot = 0
   for iplate=0L, nplate-1L do begin
      print, format='("Plate ",i4," of ",i4,a1,$)', iplate, nplate, string(13b)
      readspec, plist[iplate].plate, mjd=plist[iplate].mjd, objhdr=objhdr
      expid = sxpar(objhdr, 'EXPID*')
      ntot = ntot + n_elements(expid)
   endfor
   print

   npix = 4000
   loglam = 3.5700 + dindgen(npix)*1d-4
   efficiency = fltarr(npix,ntot)
   explist = strarr(ntot)
   airmass = fltarr(ntot)

   ;----------

   print, 'Reading exposures'
   inum = 0L
   for iplate=0L, nplate-1L do begin
      print, format='("Plate ",i4," of ",i4,a1,$)', iplate, nplate, string(13b)
      readspec, plist[iplate].plate, mjd=plist[iplate].mjd, objhdr=objhdr
      expid = sxpar(objhdr, 'EXPID*')
      camname = strmid(expid,0,2)
      for iexp=0, n_elements(expid)-1 do begin
         foo = spthroughput(plist[iplate].plate, $
          camname=strmid(expid[iexp],0,2), $
          expnum=long(strmid(expid[iexp],3,8)), loglam=loglam, $
          airmass=airmass1, efficiency=efficiency1, /median)
         efficiency[*,inum] = efficiency1
         explist[inum] = expid[iexp]
         airmass[inum] = airmass1
         inum = inum + 1L
      endfor
   endfor
   print

stop

   return
end
;------------------------------------------------------------------------------
