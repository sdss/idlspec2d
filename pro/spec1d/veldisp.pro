
;-----------------------------------------------------------
;  This is the IDL way to store FFT frequencies,
;  IDL will do odd pixel FFT's, but it's likely much slower
function getidlfreq, npix

   if (npix mod 2 EQ 0) then $
       return, [findgen(npix/2+1),-npix/2.0 + 1.0 + findgen(npix/2-1)]/npix

   return, [findgen(npix/2+1),-(npix- 1.0)/2.0 + findgen(npix/2)]/npix
end


;
;	Do we need to pass wavelengths here?  
;       Or do we do this before this proc?
;
pro veldisp, objflux, objerr, star, starerr, bestz=bestz, bestsigma=bestsigma, lowcutoff=lowcutoff
       
   if N_params() LT 2 then begin
     print, 'syntax - veldisp, flux, err, star, starerr, lowcutoff=lowcutoff'
     return
   endif

   if (NOT keyword_set(lowcutoff)) then lowcutoff = 1.0/30.0

   ;  We want to loop over objects, and just take star FFT once.
   ;  Let's do stellar FFT first
   if (size(star))[0] NE 1 then error, 'Stellar template is not 1d'

   npixstar = n_elements(star)
   starmean = mean(star)
     
   if (starmean LE 0) then error, 'star has negative mean counts'

   ; Take cosbell of normalized stellar spectra,  0.2 is cosbell fraction
   starbell = (star/starmean - 1.0) * cosbell(star, 0.2)
   starerrbell = (starerr/starmean)* cosbell(star, 0.2)
    
   ;---------------------------------------------------------------------
   pad = 2L^(fix(alog(npixstar)/alog(2) + 2)) - npixstar


   starbell = [starbell,fltarr(pad)]
   starerrbell = [starerrbell,fltarr(pad)]
   npixbell = n_elements(starbell)
   starfft = fft(starbell)  * npixbell
   starvariancefft = fft(starerrbell*starerrbell)  * npixbell
   starvar0 = float(starvariancefft[0])

   ; Now we have our star ready, let's loop out objects

   if (size(objflux))[0] EQ 1 then nobj = 1 $
   else if (size(objflux))[0] EQ 2 then nobj = (size(objflux))[2] $
   else error, 'flux array is neither 1d or 2d'

   bestsigma = fltarr(nobj)
   bestz = fltarr(nobj)

   for iobj=0,nobj-1 do begin

     if (nobj EQ 1) then begin
       tempflux = objflux
       temperr = objerr
     endif else begin
       tempflux = objflux[*,iobj]
       temperr = objerr[*,iobj]
     endelse

     npixflux = n_elements(tempflux)
     npixerr = n_elements(temperr)

     if (npixerr NE npixflux) then error, 'Flux and Error arrays do not match'


     ;
     ;	Here let's interpolate all pixels where temperr LE 0.0
     ;   We could (should) use djs_maskinterp here?
     ;

     nonzero = where(temperr GT 0.0)
     x = findgen(npixerr)
     lessthanzero = where(temperr LE 0.0)
     inrange = where(temperr LE 0.0 AND x GT min(nonzero) AND x LT max(nonzero))

     if (nonzero[0] NE -1 AND inrange[0] NE -1) then begin
       spl0 = spl_init(x[nonzero],tempflux[nonzero])
       tempflux[lessthanzero] = 0.0
       tempflux[inrange] = spl_interp(x[nonzero],tempflux[nonzero], spl0, x[inrange])
     endif

     ; Here we need to match the stellar spectrum, 
     ; after cosbell 

     objbell = cosbell(tempflux, 0.2)
     fluxmean = mean(tempflux)
 
     if (fluxmean LE 0) then begin
       print, i, 'th  galaxy has negative mean counts', galmean     
     endif else begin

       fluxbell = (tempflux/fluxmean - 1.0) * objbell
       errbell = (temperr/fluxmean)*objbell

       if npixbell LE npixflux then error, 'Not enough padding'

       fluxbell = [fluxbell, fltarr(npixbell-npixflux)]
       errbell = [errbell, fltarr(npixbell-npixflux)]
       
       fluxfft = fft(fluxbell) * npixbell
       fluxvariancefft = fft(errbell*errbell)  * npixbell
       fluxvar0 = float(fluxvariancefft[0])

       fftfreq = getidlfreq(npixbell)

       ;
       ;  Do bandpass filter of fft here
       ;

        ;  Use a cosbell for high pass filter  at 15?? pixels

       ; Should write a routine to perform band pass filtering here

        lowpixels = where(abs(fftfreq) LT lowcutoff)
        highpass = fftfreq * 0.0 + 1.0
        if (lowpixels[0] NE -1) then highpass[lowpixels] = $
           0.5 * (1.0 - cos(!Pi*abs(fftfreq[lowpixels])/lowcutoff))

       
        starfilt = starfft * highpass 
        fluxfilt = fluxfft * highpass 

        corr =  float(fft(fluxfilt * conj(starfilt),1))

        corr = shift(corr,pad)


        ;  This loop finds the redshift by searching the 5 highest peaks

        for i=0,5 do begin
           nhalffit = 5
           xtemp = findgen(2*nhalffit + 1)-nhalffit
           peak = max(corr,velcen)
           parabola = poly_fit(xtemp, corr[xtemp+velcen], 2)

;           plot, corr
;           oplot, xtemp+velcen, poly(xtemp, parabola), color=500


           if (parabola[2] GE 0.0) then begin
             print, 'peak is not well fit at ', velcen
             corr[xtemp+velcen] = 0.0
           endif else if (total(corr[xtemp+velcen]) LT 0.0) then begin
             print, 'total corr is less than zero at ', velcen
             corr[xtemp+velcen] = 0.0
           endif else i=5

        endfor

        fitcen = velcen - 0.5 * parabola[1]/parabola[2]  - pad
        bestz[iobj] = fitcen   ; this redshift is in pixels!

        twopiei = 2.0 * !Pi * complex(0.0,1.0)
        phase = exp( - twopiei * fftfreq * fitcen)
        starshift = starfft * phase

;
;	Need to pick lower and upper limits to do comparison
;	Let's try to compare from 80 pixels to 2.2 pixels
;
 
        lowlimit = 1.0/80.0
        highlimit = 1.0/2.2

        inside = where(abs(fftfreq) GT lowlimit AND abs(fftfreq) LT highlimit)

        if (inside[0] EQ -1) then begin
           print, 'no pixels in correct frequency range'
           return
        endif
	
;
;       Tonry and Davis show a simple expression to maximize
;       This is a slow minimizer, stepping through 100 sigmas to
;       find best one.  This method is the fft difference method
;       We need a routine for each method.
;

        chi2 = fltarr(100)
        sigma = fltarr(100)
        alpha = fltarr(100)
        for i=0,99 do begin
           
          sigma[i] = i/10.0;  in pixels
          broad = exp(-(fftfreq*sigma[i] * 2.0 * !Pi)^2/2.0)

          numer = fluxfft * conj(starshift) * broad
          denom = starshift * conj(starshift) * broad^2
          alpha[i] = total(numer[inside])/total(denom[inside])

          diff = fluxfft[inside] - alpha[i] * starshift[inside] * broad[inside]
;          chi2[i] = float(total(diff * conj(diff)/(fluxvar0 + starvar0 * broad[inside]^2)))
          chi2[i] = float(total(diff * conj(diff)/(fluxvar0)))

        endfor

      minchi2 = min(chi2,minplace)
      bestsigma[iobj] = sigma[minplace] 

      chip = poly_fit(sigma/10.0,chi2,10,/double)
      lotsofsigma = findgen(1000)/100.0
      chifit = poly(lotsofsigma/10.0,chip)
      minchifit = min(chifit,fitplace)
     
      bestsigma[iobj] = lotsofsigma[fitplace]
      print,bestz[iobj], bestsigma[iobj], minchi2, alpha[minplace]

      plot,sigma,chi2,ps=1,/yno
      oplot,lotsofsigma,chifit
     endelse
   endfor 

   return
end 
