
pro writespectra, objhdr, plugsort, flux, fluxivar, wset, $
 filebase=filebase

   dims = size(flux, /dimens)
   ntrace = dims[1]
   npix = dims[0]
   nparams = (size(wset.coeff,/dimens))[0]
   
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
     sxaddpar, outhdr, 'RAOBJ', plugsort[i].ra, 'RA (deg) of object'
     sxaddpar, outhdr, 'DECOBJ', plugsort[i].dec, 'DEC (deg) of object'
     sxaddpar, outhdr, 'OBJTYPE', plugsort[i].objtype
     sxaddpar, outhdr, 'PRIMTARG', plugsort[i].primtarget
     sxaddpar, outhdr, 'SECTARGE', plugsort[i].sectarget
     sxaddpar, outhdr, 'FIBERID', plugsort[i].fiberId

     ; Add wavelength cards...
     sxaddpar, outhdr, 'NWORDER', nparams
     sxaddpar, outhdr, 'WFITTYPE', 'LOG-'+wset.func
     sxaddpar, outhdr, 'PIXMIN', wset.xmin
     sxaddpar, outhdr, 'PIXMAX', wset.xmax

     for j=0, nparams-1 do begin
       keyw = 'COEFF'+string(format='(i1)',j) 
       sxaddpar, outhdr, keyw, wset.coeff[j,i]
     endfor

     outname = filebase + '-' + string(format = '(i3.3)', i+1) + '.fit'

     sxaddpar, outhdr, 'BITPIX', -32
     sxaddpar, outhdr, 'NAXIS', 2
     sxaddpar, outhdr, 'NAXIS1', npix
     sxaddpar, outhdr, 'NAXIS2', 2
     writefits, outname, [[outflux],[outerr]], outhdr
   endfor

   save, filename=filebase+'.ss', objhdr, flux, fluxivar, wset, plugsort

   return
end
