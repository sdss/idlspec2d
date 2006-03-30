
;------------------------------------------------------------------------------
pro plot_thru, plate, mjd=mjd

   ;----------
   ; Get the list of EDR-DR5 plates

   public = ['EDR','DR1','DR2','DR3','DR4','DR5']
   platelist, plist=plist
   qkeep = bytarr(n_elements(plist))
   for i=0, n_elements(public)-1 do $
    qkeep = qkeep OR strmatch(plist.public,public[i])
   qkeep = qkeep AND strmatch(plist.statuscombine,'Done*')
   plist = plist[where(qkeep, nplate)]

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
