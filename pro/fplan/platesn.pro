;
;platesn, 'spPlate-0306-51637.fits', 'spSN2d-0306-51637.ps'
;

pro platesn, platefile, snplot, planfile=planfile

   if keyword_set(planfile) then begin
      yanny_read, planfile, pdata, hdr=hdr
      platefile = yanny_par(hdr,'combinefile')
      logfile = yanny_par(hdr,'logfile')
      snplot = yanny_par(hdr,'snplot')
      yanny_free, pdata
   endif

   if keyword_set(logfile) then splog, filename=logfile, /append
   if (findfile(platefile))[0] EQ '' then begin
      splog, 'ABORT: No spPlate file exists?!?'
      return
   endif

   finalflux=mrdfits( platefile, 0, hdr, /silent)
   finalivar=mrdfits( platefile , 1, /silent)
   finalandmask  =mrdfits( platefile , 2, /silent)
   finalplugmap=mrdfits( platefile, 5, /silent)

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
   ;  that's why we add 2.5*17 to the magnitude
   ;
   waveimg = 10^(finalwave) 
   flambda2fnu = (waveimg*waveimg / 2.99792e18) # replicate(1,nfiber)

   filter=transpose(filter_thru(finalflux*flambda2fnu, waveimg=waveimg, $
             mask=(finalivar LE 0), /norm))

   synthetic_mags = fltarr(3,nfiber)
   posfilter = where(filter[1:3,*] GT 0)
   if posfilter[0] NE -1 then $
       synthetic_mags[posfilter] = $
                -2.5 * alog10((filter[1:3,*])[posfilter]) - 48.6 + 2.5*17.0

   ;---------------------------------------------------------------------------
   ;  Make S/N plot
   ;
 

   plotsn, snvec, finalplugmap, plotfile=snplot, plottitle=plottitle, $
      synthmag=synthetic_mags, snplate=snplate

   ;--------------------------------------------------------------------
   ;          Bad Fiber roll call
   ;
  
   rollcall = lonarr(3,8)

   for i=0,7 do $
      rollcall[*,i] = $
          [  total(total((finalandmask[gwave,*] AND 2L^i) GT 0, 1) GT 0), $
             total(total((finalandmask[rwave,*] AND 2L^i) GT 0, 1) GT 0), $
             total(total((finalandmask[iwave,*] AND 2L^i) GT 0, 1) GT 0)]
   
   splog, "Fiber Problem",  "  # in g'",  "  # in r'", "  # in i'" , $
          format='(3x,a20, 3(a9))'
   common com_maskbits, maskbits
   test = pixelmask_bits("NODATA")
   for i=0,7 do begin
     label = maskbits[where(maskbits.bit EQ i)].label
     info = string(format='(a20,3(i9.3))',label,rollcall[*,i])
     splog, label, rollcall[*,i], format='(a20,3(i9))'
   endfor

   if keyword_set(logfile) then splog, /close
   ;---------------------------------------------------------------------------
   ;	!!! This is crazy, but I'm going to write out the S/N 
   ;    !!!   per plate in the hdr !!
  

   if NOT keyword_set(snplate) then return

   bands = ['G','R','I']
   snmag = [20.2, 20.25, 19.9]


   for ispec=1,2 do begin
     for bb=0,n_elements(bands) - 1 do begin

         key1 = 'SPEC'+ strtrim(ispec,2)+'_'+bands[bb]
         comment = string(format='(a,i2,a,f5.2)', '(S/N)^2 for spec ', ispec, ' at mag ', snmag[bb])

         sxaddpar, hdr, key1, snplate[ispec-1,bb], comment, before='LOWREJ'
     endfor
   endfor

   ;---------------------------------------------------------------------------
   ; Write combined output file
   ;---------------------------------------------------------------------------

   finalormask=mrdfits( platefile, 3, /silent)
   finaldispersion=mrdfits( platefile, 4, /silent)

   ; 1st HDU is flux
   mwrfits, finalflux, platefile, hdr, /create

   ; 2nd HDU is inverse variance
   mwrfits, finalivar, platefile

   ; 3rd HDU is AND-pixelmask
   mwrfits, finalandmask, platefile

   ; 4th HDU is OR-pixelmask
   mwrfits, finalormask, platefile

   ; 5th HDU is dispersion map
   mwrfits, finaldispersion, platefile

   ; 6th HDU is plugmap
   mwrfits, finalplugmap, platefile

   mwrfits, snvec, platefile
   mwrfits, synthetic_mags, platefile

   return
end
