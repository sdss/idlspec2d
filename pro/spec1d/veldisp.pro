;+
; NAME:
;   veldisp
;
; PURPOSE:
;   Fit a series of galaxy spectrum with a single stellar template.
;    For each object, this procedure will first find the best redshift 
;    between object and template.  The correlation function is formed
;    in fitredshift, and the best redshift and width of the correlation
;    peak is calculated (along with error estimates).  Next perform chi2
;    fitting with fourier_difference and fourier_quotient methods
;
; CALLING SEQUENCE:
;   veldisp, objflux, objerr, star, starerr, result, lowcutoff=lowcutoff, $
;     nloop=nloop, sigmastep=sigmastep, highcutoff=highcutoff, doplot=doplot, $
;     nodiff=nodiff
;
; INPUTS:
;   objflux    - array of spectra [npix, nobj]
;   objerr     - corresponding error [npix, nobj]
;   star       - stellar template [nstarpix]
;   starerr    - corresponding error [nstarpix]
;
; OPTIONAL KEYWORDS:
;   nloop      - number of sigmas to fit in chi2 tests
;   sigmastep  - steps between sigmas to test
;   lowcutoff  - low frequency cutoff for cross-correlation peak finding
;   highcutoff - high frequency cutoff for cross-correlation peak finding
;   doplot     - Show plots of spectra, chi2 and correlation peaks
;   nodiff     - skip fourier_difference (as it's slow right now)
;
; OUTPUTS:
;   result     - structure array containing desired outputs
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;   We assume that objflux and star have the same zeropoint
;    And that all spectra are binned log-linear in wavelength
;    If there is a zeropoint difference between objects and star
;    this needs to be included after veldisp has run.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;  getidlfreq  - quick function to return idl specific frequencies in FFT
;  fft()
;  fitredshift
;  fourier_difference()
;  fourier_quotient()
;
; REVISION HISTORY:
;   25-Mar-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
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
pro veldisp, objflux, objerr, star, starerr, result, lowcutoff=lowcutoff, $
    nloop=nloop, sigmastep=sigmastep, highcutoff=highcutoff, doplot=doplot, $
    nodiff=nodiff
       
   if N_params() LT 2 then begin
     print, 'syntax - veldisp, flux, err, star, starerr, lowcutoff=lowcutoff,'
     print, '           nloop=nloop, sigmastep=sigmastep'
     return
   endif

   if (NOT keyword_set(lowcutoff)) then lowcutoff = 1.0/30.0
   if (NOT keyword_set(highcutoff)) then highcutoff = 1.0/3.0
   if (NOT keyword_set(nloop)) then nloop = 40
   if (NOT keyword_set(sigmastep)) then sigmastep = 0.2


   tempresult = { $
       z                  : 0.0, $
       z_err              : 0.0, $
       sigma_cc           : 0.0, $
       sigma_cc_err       : 0.0, $
       sigma_quotient     : 0.0, $
       sigma_quotient_err : 0.0, $
       sigma_diff         : 0.0, $
       sigma_diff_err     : 0.0 }



   ;  We want to loop over objects, and just take star FFT once.
   ;  Let's do stellar FFT first
   if (size(star))[0] NE 1 then message, 'Stellar template is not 1d'

   npixstar = n_elements(star)
   starmean = mean(star)
     
   if (starmean LE 0) then error, 'star has negative mean counts'

   ; Take cosbell of normalized stellar spectra,  0.2 is cosbell fraction

   stargood = where(starerr GT 0.0)
   starup = max(stargood)
   stardown = min(stargood)
   nfill = starup - stardown + 1
   filled = lindgen(nfill) + stardown
   
   starbell = fltarr(npixstar) 
   starerrbell = fltarr(npixstar) 

   starbell[filled] = ((star/starmean - 1.0))[filled] * cosbell(filled, 0.2)
   starerrbell[filled] = ((starerr/starmean))[filled] * cosbell(filled, 0.2)
    
   ;---------------------------------------------------------------------
   pad = 2L^(fix(alog(npixstar)/alog(2) + 2)) - npixstar


   starbell = [starbell,fltarr(pad)]
   starerrbell = [starerrbell,fltarr(pad)]
   npixbell = n_elements(starbell)
   fftfreq = getidlfreq(npixbell)
   starfft = fft(starbell)  * npixbell
   starvariancefft = fft(starerrbell*starerrbell)  * npixbell
   starvar0 = float(starvariancefft[0])

   ; Should write a routine to perform band pass filtering here

     lowpixels = where(abs(fftfreq) LT lowcutoff)
     highpass = fftfreq * 0.0 + 1.0

     if (lowpixels[0] NE -1) then highpass[lowpixels] = $
        0.5 * (1.0 - cos(!Pi*abs(fftfreq[lowpixels])/lowcutoff))

     highpixels = where(abs(fftfreq) GT highcutoff)
     lowpass = fftfreq * 0.0 + 1.0
     sep = max(fftfreq) - highcutoff
     if (highpixels[0] NE -1) then lowpass[highpixels] = $
        0.5 * (1.0 - cos(!Pi*(max(fftfreq) - abs(fftfreq[highpixels]))/sep))


   starfilt = starfft * highpass * lowpass

   fitredshift, starfilt, starfilt, pad=pad, $
      nsearch=5, zfit=starcen, z_err=starcen_err, $
      veldispfit=starsigma, veldisp_err=starsigma_err, doplot=doplot

   ; Now we have our star ready, let's loop out objects

   if (size(objflux))[0] EQ 1 then nobj = 1 $
   else if (size(objflux))[0] EQ 2 then nobj = (size(objflux))[2] $
   else error, 'flux array is neither 1d or 2d'

   result = replicate(tempresult, nobj)

   print,'    Gal    z      z_err   vel_cc (err)  vel_q  (err)  vel_d  (err)  alpha_d  alpha_q'

   for iobj=0,nobj-1 do begin

     if (nobj EQ 1) then begin
       tempflux = objflux
       temperr = objerr
     endif else begin
       tempflux = objflux[*,iobj]
       temperr = objerr[*,iobj]
     endelse

;     window,0 
;     plot, tempflux
;     oplot, temperr, color=500

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
     if (lessthanzero[0] NE -1) then tempflux[lessthanzero] = 0.0

     galup = max(nonzero)
     galdown = min(nonzero)
     nfill = galup - galdown + 1
     filled = lindgen(nfill) + galdown
    

     inrange = where(temperr LE 0.0 AND x GT galdown AND x LT galup)


     ; Here we need to match the stellar spectrum, 
     ; after cosbell 

     objbell = tempflux*0.0
     objbell = cosbell(filled, 0.2)

     fluxmean = mean(tempflux)
 
     if (fluxmean LE 0) then begin
       print, iobj, 'th  galaxy has negative mean counts', fluxmean     
     endif else begin

       fluxbell = fltarr(npixflux) 
       errbell = fltarr(npixerr)

     if (nonzero[0] NE -1 AND inrange[0] NE -1) then begin
       spl0 = spl_init(x[nonzero],tempflux[nonzero])
;       tempflux[inrange] = spl_interp(x[nonzero],tempflux[nonzero], spl0, x[inrange])
       tempflux[inrange] = fluxmean
     endif

       fluxbell[filled] = (tempflux/fluxmean - 1.0)[filled] * objbell
       errbell[filled] = (temperr/fluxmean)[filled]*objbell

       if npixbell LE npixflux then error, 'Not enough padding'

       fluxbell = [fluxbell, fltarr(npixbell-npixflux)]
       errbell = [errbell, fltarr(npixbell-npixflux)]
       
       fluxfft = fft(fluxbell) * npixbell
       fluxvariancefft = fft(errbell*errbell)  * npixbell
       fluxvar0 = float(fluxvariancefft[0])


       ;
       ;  Do bandpass filter of fft here
       ;


       
        fluxfilt = fluxfft * highpass  * lowpass


       ;
       ;  Need to fill an array of length fluxfilt which records the number
       ;   of good pixels cross-correlated between star and galaxy as a 
       ;   a function of shift

        npixtot = n_elements(fluxbell)
	ngoodpixels = lonarr(npixtot)
        for i=0, npixtot-1 do ngoodpixels[i] = $
             min([starup,galup+pad-i]) - max([stardown,galdown+pad-i])

        fitredshift, fluxfilt, starfilt, pad=pad, ngoodpixels= ngoodpixels, $
          nsearch=5, zfit=fitcen, z_err=fitcen_err, $
          veldispfit=galsigma, veldisp_err=galsigma_err, doplot=doplot

        result[iobj].z = fitcen   ; this redshift is in pixels!
        result[iobj].z_err = fitcen_err 

        if (keyword_set(doplot)) then begin
          window,2, ypos = 50
          plot, star/starmean, xr=[0,2000], yr=[0,2.0], $
            title='Rest frame spectra of template (white) and galaxy (red)'
          oplot, smooth(shift(tempflux/fluxmean,-fitcen),5) + 0.1, color=500
        endif
        if (galsigma GT starsigma AND starsigma GT 0.0) then begin
           result[iobj].sigma_cc = sqrt(galsigma^2 - starsigma^2)
           result[iobj].sigma_cc_err = sqrt((galsigma*galsigma_err)^2 + $
                (starsigma*starsigma_err)^2)/result[iobj].sigma_cc
        endif


        twopiei = 2.0 * !Pi * complex(0.0,1.0)
        phase = exp( - twopiei * fftfreq * fitcen)
        starshift = starfft * phase

	testsigma = findgen(nloop)*sigmastep

;
;	Need to pick lower and upper limits to do comparison
;	Let's try to compare from 80 pixels to 2.2 pixels
;       Don't do fourier_difference if nodiff keyword is set
;
	if (NOT keyword_set(nodiff)) then $
	answer = fourier_difference(fftfreq, fluxfft, starshift, fluxvar0, $
             starvar0, testsigma=testsigma, deltachisq=1.0, $ 
             lowlimit = 1.0/80.0, highlimit=1.0/2.2, doplot=doplot)

;---------------------------------------------------------
;          quotient = float(fluxfft/(starshift))
;          quotient = float(fluxfft/(alpha[i]*starshift))
;          diff = (quotient[inside] - broad[inside])
;          chi2quotient[i] = total(diff^2)/fluxvar0
 

       bestalpha = -9999.0

       if (n_elements(answer) EQ 4) then begin
         result[iobj].sigma_diff = answer[1]
         result[iobj].sigma_diff_err = answer[2]
         bestalpha = answer[3]
       endif

;         print, iobj, result[iobj].z, exp(result[iobj].z*0.000230259)-1.0, $
;              result[iobj].sigma_diff*70.0, result[iobj].sigma_diff_err*70.0, $
;             answer[0], bestalpha

        if (keyword_set(doplot)) then quoplot = 2 else quoplot = 0
	answerq = fourier_quotient(fftfreq, fluxfft, starshift, fluxvar0, $
             starvar0, testsigma=testsigma, deltachisq=1.0, $ 
             lowlimit = 1.0/80.0, highlimit=1.0/2.2, doplot=quoplot)

       if (n_elements(answerq) EQ 4) then begin
         result[iobj].sigma_quotient = answerq[1]
         result[iobj].sigma_quotient_err  = answerq[2]
         bestalpha_q = answerq[3]
       endif

         print,iobj,result[iobj], format='(i,f9.2,7(f7.2),$)' 
         print, bestalpha, bestalpha_q, format='(f10.4,f9.4)'

;      window,1 
;      plot,sigma,chi2diff,ps=1,/yno
;      oplot,lotsofsigma,chifit, color=500
     endelse
   endfor 

   return
end 
