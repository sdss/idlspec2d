;+
; NAME:
;   mveldisp
;
; PURPOSE:
;   Fit a series of galaxy spectrum with a single stellar template.
;    For each object, this procedure will first find the best redshift 
;    between object and template.  The correlation function is formed
;    in mfitredshift, and the best redshift and width of the correlation
;    peak is calculated (along with error estimates).  Next perform chi2
;    fitting with mfourier_difference and mfourier_quotient methods
;
; CALLING SEQUENCE:
;   mveldisp, objflux, objivar, starflux, starivar, result, $
;    klo_cut=, khi_cut=, maxsig=, sigmastep=, /doplot, /nodiff ]
;
; INPUTS:
;   objflux    - Array of object spectra [npix, nobj]
;   objivar    - Array of object inverse variance [npix, nobj]
;   starflux   - Template spectrum [nstarpix]
;   starivar   - Template inverse variance [nstarpix]
;
; OPTIONAL KEYWORDS:
;   klo_cut    - Low frequency cutoff for cross-correlation peak finding;
;                default to 1/30.
;   khi_cut    - High frequency cutoff for cross-correlation peak finding;
;                default to 1/3.
;   maxsig     - Maximum velocity dispersion to search for; default to 2 pix
;   sigmastep  - Steps between each test sigma; default to 0.2 pix
;   doplot     - Show plots of spectra, chi2 and correlation peaks
;   nodiff     - skip mfourier_difference (as it's slow right now)
;
; OUTPUTS:
;   result     - Structure array with outputs
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;   We assume that objflux and star have the same zeropoint
;   And that all spectra are binned log-linear in wavelength
;   If there is a zeropoint difference between objects and star
;   this needs to be included after mveldisp has run.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;  mveldisp_fft
;  djs_maskinterp()
;  djs_mean()
;  fft_apodize
;  mfitredshift
;  mfourier_difference()
;  mfourier_quotient()
;
; REVISION HISTORY:
;   25-Mar-2000  Written by S. Burles, FNAL
;   29-Mar-2000  Modified by D. Finkbeiner & D. Schlegel, APO
;   25-Jun-2000  Cleaned up and commented by D. Finkbeiner, APO
;   2000-Sep-11  Worked-over by the SWAT team
;-
;------------------------------------------------------------------------------
pro mveldisp, objflux, objivar, objwave, starflux, starivar, starwave, $
 result, redshifts= redshifts, czmin=czmin, czmax=czmax, klo_cut=klo_cut, $
 khi_cut=khi_cut, maxsig=maxsig, sigmastep=sigmastep, doplot=doplot, $
 nodiff=nodiff, noquotient=noquotient, nobe=nobe, keep=keep, $
 starnoise=starnoise, continuum=continuum, noxcorr=noxcorr

; set keyword defaults
   if (NOT keyword_set(klo_cut)) then klo_cut = 1.0/128.
   if (NOT keyword_set(khi_cut)) then khi_cut = 1.0/3.0
   if (NOT keyword_set(maxsig)) then maxsig = 6.0
   if (NOT keyword_set(sigmastep)) then sigmastep = 0.2
   IF (keyword_set(nobe)) THEN BEGIN 
       nodiff = 1 &  noquotient=1
   ENDIF 
; prepare plot windows
   IF (keyword_set(doplot)) THEN BEGIN
      window, 0, retain=2 &  window, 1, retain = 2 &  window, 2, retain =2
   ENDIF

   if (size(objflux, /tname) EQ 'DOUBLE') then PI = !dpi $
    else PI = !pi

; check dimensions of everything
   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   if (ndim EQ 1) then begin
      nobj = 1
      npixobj = dims[0]
   endif else if (ndim EQ 2) then begin
      npixobj = dims[0]
      nobj = dims[1]
   endif else begin
      message, 'OBJFLUX is neither 1-D or 2-D'
   endelse

   if total(abs(size(starflux, /dimens)-size(starivar, /dimens))) NE 0 $
    OR size(starflux, /n_dimen) NE size(starivar, /n_dimen) THEN  $
    message, 'Dimensions of STARFLUX and STARIVAR do not match'

   if total(abs(size(objflux, /dimens)-size(objivar, /dimens))) NE 0 $
    OR size(objflux, /n_dimen) NE size(objivar, /n_dimen) THEN  $
    message, 'Dimensions of OBJFLUX and OBJIVAR do not match'

   nstar = n_elements(starflux)/(size(starflux))[1]	
	
   tempresult = {VELDISP,   $
       z                  : fltarr(nstar), $
       z_err              : fltarr(nstar), $
       zconf              : fltarr(nstar), $
       sigma_cc           : intarr(nstar), $
       sigma_cc_err       : intarr(nstar), $
       sigma_quotient     : intarr(nstar), $
       sigma_quotient_err : intarr(nstar), $
       sigma_quot_diff    : intarr(nstar), $
       sigma_quot_diff_err : intarr(nstar), $
       sigma_diff         : intarr(nstar), $
       sigma_diff_err     : intarr(nstar) }
   result = replicate(tempresult, nobj)

   ;---------------------------------------------------------------------------
   ; Decide how large the padded spectra should be, based upon the
   ; large of the size of STARFLUX and OBJFLUX.
   ; Pad to larger (or equal) 2^N, and then doubled for isolated b.c.


;;;;; <- Mariangela's idea of a comment

   print, ' nobj nstar  Deltaredshift    veldisp_cc veldisp_q_d veldisp_q  veldisp_d  redshift'


     khi_cut_inp=khi_cut

   FOR iobj=0, nobj-1 DO BEGIN


    if (keyword_set(continuum)) then begin
      npix=(size(objflux[*, iobj]))[1]
      x=findgen(npix)*0.0001
      galerr = 1.0/objivar[*, iobj]^2
      inf=where(finite(galerr) EQ 0)
      if (inf[0] NE -1) then galerr[inf]=0
      fullbkpt = slatec_splinefit(x,objflux[*, iobj] , coeff, nbkpt=5, invvar=galerr,upper=2.8, lower=0.1)
      continuum = slatec_bvalu(x, fullbkpt, coeff)	
      objflux[*, iobj]=objflux[*, iobj]/continuum
      objflux[*, iobj]=djs_maskinterp(objflux[*, iobj], objivar[*, iobj] LE 0.0, /const)
      objivar[*, iobj]=objivar[*, iobj]/continuum
      normobj = djs_mean(objflux[*, iobj]) > djsig(objflux[*, iobj])
      objflux[*, iobj]=objflux[*, iobj]/ normobj -1
      objivar[*, iobj]=objivar[*, iobj]/ normobj
      endif
;	plot, objflux[*, iobj], ps=0
;stop		

      fluxobj= objflux[*, iobj]	
      fluxivar = objivar[*, iobj]
      waveobj=objwave[*,iobj]

    if (keyword_set(redshifts)) then begin    
	if (NOT keyword_set(czmin)) then czmin=3000
	if (NOT keyword_set(czmax)) then czmax=10000	
;	galsize=(size(objflux))[1]
;        xp=fltarr(galsize)
;        xp=3.57403 + findgen(galsize)*0.0001
        xp=objwave[*,iobj]
        xpp=xp-alog10(1+redshifts[iobj])
        ppp=where(10^xpp GT czmin AND 10^xpp LT czmax)
; cut off last 300 km/s
        czmaxgal=max(10^xpp[ppp]) - 300
        czmingal=min(10^xpp[ppp])
        if (redshifts[iobj] GT 0.16) then begin
        if (czmaxgal GT 6000) then czmaxgal= 6000
        czmingal= 4000
	endif 
        ppp1=where(10^xpp GT czmingal AND 10^xpp LT czmaxgal)

; remove sky lines
        skyremove, objflux[ppp1, iobj], xp[ppp1], $
                   5570, 5590, ppp2, kcor
        ppp1=ppp1[ppp2]
        kcor5570=kcor 
        skyremove, objflux[ppp1, iobj], xp[ppp1], $
                   6295, 6310, ppp3, kcor
        ppp1=ppp1[ppp3]
        kcor6295=kcor

	fluxobj= objflux[ppp1, iobj]	
;        if (redshifts[iobj] GT 0.2) then fluxobj=smooth(fluxobj,8)
	fluxivar = objivar[ppp1, iobj]
        waveobj= 10^xpp[ppp1]
;	plot, waveobj, fluxobj,ps=0
;stop
    endif


;--------------------------------------------------------------------------- 
; Compute FFT for stellar template

;   nstar = n_elements(starflux)/(size(starflux))[1]

   FOR istar=0, nstar-1 DO BEGIN 

t1=systime(1)
print,'BEGIN',systime(1)-t1

    meanstar=mean(starflux[*, istar]) 
    if (keyword_set(starnoise)) then begin $ &
        starflux[*, istar] = starflux[*, istar] + randomu(105,3918,/normal)*meanstar/starnoise
        starivar[*, istar]=sqrt(starivar[*, istar]^2+(meanstar/starnoise)^2)
    endif
    if (NOT keyword_set(starnoise)) then begin
        starflux[*, istar] = starflux[*, istar]
        starivar[*, istar]= starivar[*, istar]
    endif


    if (keyword_set(continuum)) then begin
      npix=(size(starflux[*, istar]))[1]
      x=findgen(npix)*0.0001
      galerr = 1.0/starivar[*, istar]^2
      inf=where(finite(galerr) EQ 0)
      if (inf[0] NE -1) then galerr[inf]=0
      fullbkpt = slatec_splinefit(x,starflux[*, istar] , coeff, nbkpt=5, invvar=galerr,upper=2.8, lower=0.1)
      continuum = slatec_bvalu(x, fullbkpt, coeff)	
      starflux[*, istar]=starflux[*, istar]/continuum
      starflux[*, istar]=djs_maskinterp(starflux[*, istar], starivar[*, istar] LE 0.0, /const)
      starivar[*, istar]=starivar[*, istar]/continuum
      normstar = djs_mean(starflux[*, istar]) > djsig(starflux[*, istar])
      starflux[*, istar]=starflux[*, istar]/ normstar -1
      starivar[*, istar]=starivar[*, istar]/ normstar
      endif	

;	plot, starflux[*, istar], ps=0
		
      fluxstar= starflux[*, istar]	
      fluxstarivar = starivar[*, istar]
      wavestar=starwave[*,istar]

    if (keyword_set(redshifts)) then  begin   
	starsize=(size(starflux))[1]
;        xp=fltarr(starsize)
;        xp=3.57403 + findgen(starsize)*0.0001
        xp=starwave[*,istar]
        pp=where(10^xp GT czmingal AND 10^xp LT czmaxgal)

        if (kcor5570 EQ 1) then begin
        xxp=xp[pp]+alog10(1+redshifts[iobj]) 
	cpp=where(10^xxp LT 5570 OR 10^xxp GT 5590)
        pp=pp[cpp]
	endif

	if (kcor6295 EQ 1) then begin
        xxp=xp[pp]+alog10(1+redshifts[iobj]) 
	cpp=where(10^xxp LT 6295 OR 10^xxp GT 6310)
        pp=pp[cpp]
	endif

	fluxstar= starflux[pp, istar]	
;	if (redshifts[iobj] GT 0.2) then fluxstar=smooth(fluxstar,8)
	fluxstarivar = starivar[pp, istar]
	wavestar=10^xp[pp]
;	plot,wavestar, fluxstar,ps=0
;stop
    endif
;	fluxobj1=fluxobj/mean(fluxobj)
;	if (keyword_set(continuum)) then plot,10^xpp[ppp],fluxobj,ps=0 else $
;        plot,10^xpp[ppp],fluxobj1,ps=0  	
;	fluxstar1=fluxstar/mean(fluxstar)
;	if (keyword_set(continuum)) then djs_oplot,10^xp[pp],fluxstar,ps=0,color='red' $
;        else djs_oplot,10^xp[pp],fluxstar1,ps=0,color='red'
;stop

   npixobj=(size(fluxobj))[1]
   npixstar = (size(fluxstar))[1]
   npixbig = 2L^(fix(alog(npixstar > npixobj)/alog(2) + 1.9999))

        khi_cut = khi_cut_inp

	mveldisp_fft, fluxobj, fluxivar, npixbig,  $
          fluxfft, fluxfilt, fluxvar0, fluxvariancefft, fluxivar_pad,  $
          khicut, klo_cut=klo_cut, khi_cut=khi_cut

        khi_cut=khicut


   mveldisp_fft, fluxstar, fluxstarivar, npixbig, starfft,  $
     starfilt, starvar0, starvariancefft, starivar_pad, $
     khigh, klo_cut=klo_cut, khi_cut=khi_cut, wave=starwave

   mfitredshift, starfilt, starivar_pad, starfilt, starivar_pad, $
     nsearch=5, zfit=starcen, z_err=starcen_err, $
     veldispfit=starsigma, veldisp_err=starsigma_err, doplot=doplot

   ;---------------------------------------------------------------------------
   ; LOOP OVER OBJECT SPECTRA




;   FOR iobj=0, nobj-1 DO BEGIN 


;      fluxivar = objivar[*, iobj]


	mveldisp_fft, fluxobj, fluxivar, npixbig,  $
          fluxfft, fluxfilt, fluxvar0, fluxvariancefft, fluxivar_pad,  $
          khigh, klo_cut=klo_cut, khi_cut=khi_cut
	
      if not keyword_set(noxcorr) then begin 
      mfitredshift, fluxfilt, fluxivar_pad, starfilt, starivar_pad, $
        nsearch=5, zfit=fitcen, z_err=fitcen_err, $
        veldispfit=galsigma, veldisp_err=galsigma_err, zconf=zconf, $
        doplot=doplot

      if (keyword_set(redshifts)) then $
      result[iobj].z[istar] = redshifts[iobj] + 10.^(fitcen/10000)-1. $
	else $
      result[iobj].z[istar]     = 10.^(fitcen/10000)-1. ; dimensionless z
      result[iobj].z_err[istar] = alog(10)*1e-4*fitcen_err* $
        (1+result[iobj].z[istar])

      result[iobj].zconf[istar] = zconf

      if (keyword_set(doplot)) then begin
         window,2,retain=2,XSIZE=800, YSIZE=400
           title='Rest frame spectra of template (white) and galaxy (red)'
	if (NOT keyword_set(redshifts)) then begin
        fluxobj= shift(fluxobj, -fitcen)
        waveobj= shift(waveobj, -fitcen)
        endif
	fluxobj1=fluxobj/mean(fluxobj)
	fluxstar1=fluxstar/mean(fluxstar)
	if (keyword_set(continuum)) then plot,wavestar,fluxstar,ps=0 $
        else plot,wavestar,fluxstar1,ps=0,xr=[5050,5250]
	if (keyword_set(continuum)) then djs_oplot,waveobj,fluxobj,ps=0,color='red'$
        else  oplot,waveobj,fluxobj1,ps=0,color=150

      endif

; Should really store sigma squared, and allow negative values; error
; should reflect it - DPF ???
      if (galsigma GT starsigma AND starsigma GT 0.0) then begin
         result[iobj].sigma_cc[istar] = sqrt(galsigma^2 - starsigma^2)*70

; fix this
;         result[iobj].sigma_cc_err[istar] = sqrt((galsigma*galsigma_err)^2 + $
;          (starsigma*starsigma_err)^2)/(result[iobj].sigma_cc[istar]/70)*70
         result[iobj].sigma_cc_err[istar] = sqrt((galsigma_err)^2 + $
          (starsigma_err)^2)*70
      endif
      ENDIF

      twopiei = 2.0 * PI * complex(0.0,1.0)
      knums = fft_wavenums(npixbig)
; fitcen generated from passed redshift
;      fitcen = alog10(redshifts[iobj]+1)*1E4  ; in pixels

;      phase = exp( - twopiei * knums * fitcen)
;      starshift = starfft * phase
       starshift = starfft

      testsigma = findgen(ceil(float(maxsig)/sigmastep) + 1) * sigmastep

; Need to pick lower and upper limits to do comparison
; Let's try to compare from 80 pixels to 2.2 pixels

print,'BEFORE newdiff  ',systime(1)-t1

      if (NOT keyword_set(nodiff)) then $
       answer = newdiff(fluxfft, starshift, fluxvar0, $
                starvar0, testsigma=testsigma, deltachisq=1.0, $
                lowlimit = 1.0/80.0, highlimit=1.0/5.0, $
		broadarr=broadarr, doplot=doplot)

      bestalpha = -9999.0
print,'END',systime(1)-t1

      if (n_elements(answer) EQ 4) then begin
         result[iobj].sigma_diff[istar] = answer[1]*70
         result[iobj].sigma_diff_err[istar] = answer[2]*70
         bestalpha = answer[3]
      endif

      if (keyword_set(doplot)) then quoplot = 2 else quoplot = 0
      IF NOT keyword_set(noquotient) THEN BEGIN 
         answerq = mfourier_quotient(fluxfft, starshift, fluxvar0, $
                  starvar0, testsigma2=testsigma, deltachisq=1.0, $
                  lowlimit = 1.0/80.0, highlimit=1.0/4., doplot=quoplot, $
                                     broadarr=broadarr)
      ENDIF 

      if (n_elements(answerq) EQ 4) then begin
;         result[iobj].sigma_quotient[istar] = sqrt(answerq[1])*70
;         result[iobj].sigma_quotient_err[istar]  = $
;                         (sqrt(answerq[1]+answerq[2])-sqrt(answerq[1]))*70
          result[iobj].sigma_quotient[istar] = answerq[1]*70
          result[iobj].sigma_quotient_err[istar] = answerq[2]*70
         bestalpha_q = answerq[3]
      endif


        if (keyword_set(doplot)) then qtoplot = 3 else qtoplot = 0
        IF NOT keyword_set(noquotient) THEN BEGIN 
         answerq = mfourier_quot_diff(fluxfft, starshift, fluxvar0, $
                  starvar0, testsigma2=testsigma, deltachisq=1.0, $
                  lowlimit = 1.0/80.0, highlimit=1.0/5., doplot=qtoplot, $
                                     broadarr=broadarrqt)
      if (n_elements(answerq) EQ 4) then begin
         result[iobj].sigma_quot_diff[istar] = answerq[1]*70
         result[iobj].sigma_quot_diff_err[istar]  = answerq[2]*70
          bestalpha_q_d = answerq[3]
      endif
        ENDIF

          r = result[iobj]
          print, iobj, istar, r.z[istar], r.z_err[istar], r.sigma_cc[istar], $
            r.sigma_cc_err[istar], r.sigma_quot_diff[istar], r.sigma_quot_diff_err[istar], $
            r.sigma_quotient[istar], r.sigma_quotient_err[istar], r.sigma_diff[istar], $
            r.sigma_diff_err[istar],redshifts[iobj],$
            format='(2(i4),x,f8.5," +/-",f8.5,4(2x,i3," +-",i3),x,f8.5)'


   endfor
   endfor
   return
end
