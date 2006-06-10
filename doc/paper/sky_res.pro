
 outfile = 'sky_res_v5.fits'
 path = 0
 tt = mrdfits('platelist-public.fits',1)
 n = n_elements(tt)
 chi =0 
 chihist = fltarr(241,n)
 chirej = fltarr(241,n)
 chilrg = fltarr(241,n)
 chilrg_rej = fltarr(241,n)
 lrgflag = 2L^5 + 2L^26
 for i=0,n-1 do begin
   readspec, tt[i].plate, mjd=tt[i].mjd, plug=plug
   sky =  where(strtrim(plug.objtype,2) EQ 'SKY', nsky)
   if nsky GT 0 then begin
     readspec, tt[i].plate, sky+1, mjd=tt[i].mjd, flux=flux, invvar=invvar, $
       ormask=ormask
  
     res = flux * sqrt(invvar) 
     good = where(invvar GT 0, ngood)
     goodrej = where((invvar GT 0) AND ((ormask AND 2L^7) EQ 0), ngoodrej)

     chihist[*,i] = histogram(res[good], min=-12.05, max=12.0499, bins=0.1)
     chirej[*,i] = histogram(res[goodrej], min=-12.05, max=12.0499, bins=0.1)
   endif
 
   lrg =  where((plug.primtarget AND lrgflag), nlrg)
   if nlrg EQ 0 then continue

     readspec, tt[i].plate, lrg+1, mjd=tt[i].mjd, flux=flux, invvar=invvar, $
       ormask=ormask, synflux=synflux
     res = (flux-synflux) * sqrt(invvar) 
     good = where(invvar GT 0 , ngood)
     goodrej = where((invvar GT 0) AND ((ormask AND 2L^23) EQ 0), ngoodrej)
     chilrg[*,i] = histogram(res[good], min=-12.05, max=12.0499, bins=0.1)
     chilrg_rej[*,i] = histogram(res[goodrej], $
                min=-12.05, max=12.0499, bins=0.1)

   print, i
endfor

   mwrfits, chihist, outfile, /create
   mwrfits, chirej, outfile
   mwrfits, chilrg, outfile
   mwrfits, chilrg_rej, outfile

end

