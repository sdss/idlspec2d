
pro writespectra, objhdr, plugsort, flux, fluxivar, wset, $
 filebase=filebase

   dims = size(flux, /dimens)
   ntrace = dims[1]
   
   for i=0, ntrace-1 do begin

     ; Copy flux and error for this object
     outflux = flux[*,i]
     outerr = 0*outflux - 1.0
     good = where(fluxivar[*,i] GT 0.0)
     if (good[0] NE -1) then $
      outerr[good] = 1.0/sqrt(fluxivar[good,i])

     ; If this fiber is unplugged, then set the flux equal to zero
     if (plugsort[i].fiberid LE 0) then begin
        outflux = 0 * outflux
        outerr = 0 * outerr - 1.0
     endif

     ; Modify the header
     outhdr = objhdr
     sxaddpar, outhdr, 'OBJID', string(format='(5(i))', $
                plugsort[i].objid)
     sxaddpar, outhdr, 'MAG', string(format='(5(f8.3))', $
                plugsort[i].mag)
     sxaddpar, outhdr, 'RA', plugsort[i].ra
     sxaddpar, outhdr, 'DEC', plugsort[i].dec
     sxaddpar, outhdr, 'OBJTYPE', plugsort[i].objtype
     sxaddpar, outhdr, 'PRIMTARG', plugsort[i].primtarget
     sxaddpar, outhdr, 'SECTARGE', plugsort[i].sectarget
     sxaddpar, outhdr, 'FIBERID', plugsort[i].fiberId

     ; Add wavelength cards...
     nparams = (size(wset,/dimens))[0]
     sxaddpar, outhdr, 'NWORDER', nparams
     sxaddpar, outhdr, 'WFITTYPE', 'LOG-'+wset.func
     sxaddpar, outhdr, 'PIXMIN', wset.xmin
     sxaddpar, outhdr, 'PIXMAX', wset.xmax

     for j=0, nparams-1 do begin
       keyw = 'COEFF'+string(format='(i1)',j) 
       sxaddpar, outhdr, keyw, wset.coeff[j,i]
     endfor

     outname = filebase + '-' + string(format = '(i3.3)', i+1) + '.fit'
     writefits, outname, [[outflux],[outerr]], outhdr
   endfor

   save, filename=filebase+'.ss', objhdr, flux, fluxivar, wset, plugsort

   return
end
