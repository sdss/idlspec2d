
 outfile = 'sky_cross.fits'
 path = 0
 tt = mrdfits('platelist-public.fits',1)
 n = n_elements(tt)
 fullstr =0 
 for i=0,n-1 do begin
   readspec, tt[i].plate, mjd=tt[i].mjd, plug=plug
   sky =  where(strtrim(plug.objtype,2) EQ 'SKY', nsky)

   k = where((sky GT 0 AND sky LT 319) OR (sky GT 320 AND sky LT 639),nk)

   if nk EQ 0 then continue

   readspec, tt[i].plate, sky[k]+1, mjd=tt[i].mjd, flux=flux, invvar=invvar, $
         ormask=ormask
   readspec, tt[i].plate, sky[k], mjd=tt[i].mjd, flux=fl, invvar=il, $
         ormask=orl
   readspec, tt[i].plate, sky[k]+2, mjd=tt[i].mjd, flux=fr, invvar=ir, $
         ormask=orr
 
   coeff_str = replicate({coeff1 : 0., coeff2 : 0., $
                          mag1 : fltarr(5), mag2 : fltarr(5)}, nk)
   mask1 = (il GT 0 AND invvar GT 0 AND $
           ((orl AND 2L^23) EQ 0) AND ((ormask AND 2L^23) EQ 0))
   denom1 = total(fl*fl*mask1,1)
   coeff_str.coeff1 = total(flux*fl*mask1,1) / $
             (denom1 + (denom1 EQ 0)) * (denom1 GT 0)
   coeff_str.mag1 = plug[sky[k]-1].mag

   mask2 = (ir GT 0 AND invvar GT 0 AND $
           ((orr AND 2L^23) EQ 0) AND ((ormask AND 2L^23) EQ 0))
   denom2 = total(fr*fr*mask2,1)
   coeff_str.coeff2 = total(flux*fr*mask2,1) / $
             (denom2 + (denom2 EQ 0)) * (denom2 GT 0)
   coeff_str.mag2 = plug[sky[k]+1].mag

   str = struct_addtags(plug[sky[k]], coeff_str)

   fullstr = struct_append(fullstr, str)

 endfor

 mwrfits, fullstr, outfile, /create

 allmag = [fullstr.mag1[2], fullstr.mag2[2]] 
 allcoeff = [fullstr.coeff1, fullstr.coeff2]
 keep = where(allcoeff NE 0 AND allmag GT 14.5 AND allmag LT 20)

 dfpsplot, 'sky_cross.ps', /color, /square
; djs_plot, allmag[keep], allcoeff[keep], ps=3, yr=[-0.005,0.005],$
;       xr=[14.5,20],/xs
 hogg_scatterplot, allmag[keep], allcoeff[keep], cthick=4., $
      /cond, xr=[14.5,20], yr=[-0.006, 0.006], xnpix=11, /noerase,/xs, $
      xtitle='Sky neighbor r magnitude', ytitle='Linear coefficient', $
      ythick=4, xthick=4, charth=4., chars=1.5, levels=levels, ynpix=48, $
      /nogreyscale
 dfpsclose
 
end

