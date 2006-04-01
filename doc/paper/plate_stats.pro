
   plist = mrdfits('platelist-public.fits',1)
   nplate = n_elements(plist)

   readspec, plist[0].plate, mjd=plist[0].mjd, zans=zans
   zall = replicate(zans[0], 640, nplate)
   for i=1,nplate-1 do begin & $
     readspec, plist[i].plate, mjd=plist[i].mjd, zans=zans & $
     zall[*,i] = zans & $
   endfor

   pl = plist.plate
   s = sort(pl)
   u = uniq(pl[s])
  
   u = s[u]

   print, 'Number of unique plates', n_elements(u)
   print, 'Total number of exposures', long(total(plist.nexp_b1))
   print, 'Total on-sky exposure time (hrs)', total(plist.expt_b1)/3600.

   ra = zall.plug_ra
   dec = zall.plug_dec

   ingroup = spheregroup(ra, dec, 3.0/3600., firstgroup=f)
   f = where(f NE -1)
   zf = zall[f]


   print, 'Unique objects ', n_elements(f)
   print, 'Sky fibers', long(total(strtrim(zf.objtype,2) EQ 'SKY'))
   print, 'QSO fibers', long(total(strtrim(zf.class,2) EQ 'QSO'))
   print, 'Galaxy fibers', long(total(strtrim(zf.class,2) EQ 'GALAXY'))
   print, 'Star fibers', long(total(strtrim(zf.class,2) EQ 'STAR'))
  

 
   dfpsplot, 'plate_stats.ps', /square, /color
   !p.thick=2

   g = where(plist.seeing50 GT 0.8)
   djs_plot, plist[g].mjd-51600, plist[g].platesn2/plist[g].expt_b1 * 3600, /ylog, yr=[0.1,100], $
       ps=1, xtitle='MJD - 2451600', ytitle='RMS Guiding("), Median Seeing ("), Normalized Plate S/N^2', $
       charsize=1.5, charthick=4, xthick=4, ythick=4
   oplot, plist[g].mjd-51600, plist[g].seeing50, ps=2   
   oplot, plist[g].mjd-51600, plist[g].rmsoff50, ps=4   

   for i=0,5 do oplot, [140,140] + i*365 + (i GE 4), [0.1,100] , thick=4
   dfpsclose

