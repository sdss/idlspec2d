; Plot the throughput as a function of MJD
;------------------------------------------------------------------------------
pro plot_thru, plate, mjd=mjd

   savefile = 'plot_thru.ss'
   plotfile = 'plot_thru.ps'

   mjd_mirror = [52144,52516,52900,53243,53609] ; realuminizations
   ; Ignore plates where we don't have uber-photometry
   badplate = [1038,1292,2041,2042,2043,2044,2061,2062,2063,2064,$
    2183,2211+lindgen(34), $
    2247,2256,2234,2335,2336,2339,2340,2368+lindgen(9),2426+lindgen(15)]
; Add to this list other special plates...
badplate = [1038,1292,1490,1665,1666, $
 2041,2042,2043,2044,2061,2062,2063,2064,$
 2183,2184,2211+lindgen(34), $
 2247,2256,2234,2335,2336,2339,2340,2368+lindgen(9),2426+lindgen(15)]

   if (findfile(savefile)) then begin
      restore, savefile
   endif else begin
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
         print, format='("Plate ",i4," of ",i4,a1,$)', $
          iplate, nplate, string(13b)
         readspec, plist[iplate].plate, mjd=plist[iplate].mjd, objhdr=objhdr
         expid = sxpar(objhdr, 'EXPID*')
         ntot = ntot + n_elements(expid)
      endfor
      print

      npix = 4000
      loglam = 3.5700 + dindgen(npix)*1d-4
      efficiency = fltarr(npix,ntot)
      explist = strarr(ntot)
      platelist = lonarr(ntot)
      mjdlist = lonarr(ntot)
      proglist = strarr(ntot)
      airmass = fltarr(ntot)

      ;----------

      print, 'Reading exposures'
      inum = 0L
      for iplate=0L, nplate-1L do begin
         print, format='("Plate ",i4," of ",i4,a1,$)', $
          iplate, nplate, string(13b)
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
            platelist[inum] = plist[iplate].plate
            mjdlist[inum] = plist[iplate].mjd
            proglist[inum] = plist[iplate].progname
            airmass[inum] = airmass1
            inum = inum + 1L
         endfor
      endfor
      print
      save, file=savefile, plist, loglam, efficiency, $
       platelist, mjdlist, airmass, explist

   endelse

   ntot = n_elements(platelist)
   qgood = bytarr(ntot)
   for i=0L, ntot-1L do qgood[i] = total(platelist[i] EQ badplate) EQ 0

   ; Convolve the efficiencies with the filter curves
   res = filter_thru(efficiency, waveimg=10.^loglam)

   dfpsplot, plotfile, /square, /encap
   ii = where(strmatch(explist,'r2*') AND qgood)
   djs_plot, mjdlist[ii], res[ii,3], psym=3, charsize=1.3, charthick=2, $
    xrange=[51600,53600], /xstyle, yrange=[0,0.2], /ystyle, $
    xtitle='MJD', ytitle='r2 Efficiency in i-band', xtickformat='(i5)'
   for imjd=0, n_elements(mjd_mirror)-1 do $
    djs_oplot, mjd_mirror[imjd]+[0,0], !y.crange
   dfpsclose

   return
end
;------------------------------------------------------------------------------
