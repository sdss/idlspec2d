pro platesn, platefile, snplot

   finalflux=mrdfits( platefile, 0, hdr)
   finalivar=mrdfits( platefile , 1)
   finalandmask=mrdfits( platefile, 2)
   finalormask=mrdfits( platefile, 3)
   finaldispersion=mrdfits( platefile, 4)
   finalplugmap=mrdfits( platefile, 5)

   nfiber = (size(finalflux))[2]
   npix   = (size(finalflux))[1]
   finalwave = sxpar(hdr,'COEFF0') + sxpar(hdr,'COEFF1')*findgen(npix)
   
   gwave = where(finalwave GT alog10(4000) AND finalwave LT alog10(5500))
   rwave = where(finalwave GT alog10(5600) AND finalwave LT alog10(6900))
   iwave = where(finalwave GT alog10(6910) AND finalwave LT alog10(8500))

   sn = finalflux*sqrt(finalivar)
   snvec = [ transpose(djs_median(sn[gwave,*],1)), $
             transpose(djs_median(sn[rwave,*],1)), $
             transpose(djs_median(sn[iwave,*],1))]

   ;--------------------------------------------------------------------
   ;  spectra are already in 10^-17 flambda
   ;
   waveimg = finalwave # replicate(1,nfiber)
   flambda2fnu = (10^(waveimg))^2 / 2.99792e11

   filter=transpose(filter_thru(finalflux*flambda2fnu, waveimg=waveimg, $
             mask=(finalivar GT 0), /norm))

   synthetic_mags = fltarr(3,nfiber)
   posfilter = where(filter[1:3,*] GT 0)
   if posfilter[0] NE -1 then $
       synthetic_mags[posfilter] = $
                -2.5 * alog10((filter[1:3,*])[posfilter]) - 48.6 + 60.0

   ;---------------------------------------------------------------------------
   ;  Make S/N plot
   ;
 
   plotsn, snvec, finalplugmap, plotfile=snplot, plottitle=plottitle

   ;---------------------------------------------------------------------------
   ; Write combined output file
   ;---------------------------------------------------------------------------

   ; 1st HDU is flux
;   mwrfits, finalflux, platefile, hdr, /create
;
   ; 2nd HDU is inverse variance
;   mwrfits, finalivar, platefile

   ; 3rd HDU is AND-pixelmask
;   mwrfits, finalandmask, platefile

   ; 4th HDU is OR-pixelmask
;   mwrfits, finalormask, platefile

   ; 5th HDU is dispersion map
;   mwrfits, finaldispersion, platefile

   ; 6th HDU is plugmap
;   mwrfits, finalplugmap, platefile

;   mwrfits, snvec, platefile
;   mwrfits, synthetic_mags, platefile

    stop

   return
end
