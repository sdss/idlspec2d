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
;   veldisp, objflux, objivar, starflux, starivar, result, $
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
;   nodiff     - skip fourier_difference (as it's slow right now)
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
;   this needs to be included after veldisp has run.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;  bandpassfilter()
;  djs_maskinterp()
;  djs_mean()
;  fft_apodize
;  fitredshift
;  fourier_difference()
;  fourier_quotient()
;
; REVISION HISTORY:
;   25-Mar-2000  Written by S. Burles, FNAL
;   29-Mar-2000  Modified by D. Finkbeiner & D. Schlegel, APO
;   25-Jun-2000  Cleaned up and commented by D. Finkbeiner, APO
;-
;------------------------------------------------------------------------------
pro veldisp, objflux, objivar, objwave, starflux, starivar, starwave, result, klo_cut=klo_cut, $
 khi_cut=khi_cut, maxsig=maxsig, sigmastep=sigmastep, doplot=doplot, $
 nodiff=nodiff, noquotient=noquotient, nobe=nobe, keep=keep
       
; set keyword defaults
   if (NOT keyword_set(klo_cut)) then klo_cut = 1.0/128.
   if (NOT keyword_set(khi_cut)) then khi_cut = 1.0/3.0
   if (NOT keyword_set(maxsig)) then maxsig = 2.0
   if (NOT keyword_set(sigmastep)) then sigmastep = 0.2
   IF (keyword_set(nobe)) THEN BEGIN 
       nodiff = 1 &  noquotient=1
   ENDIF 
; prepare plot windows
   IF (keyword_set(doplot)) THEN BEGIN
      window, 0 &  window, 1 &  window, 2
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


   ;---------------------------------------------------------------------------
   ; Decide how large the padded spectra should be, based upon the
   ; large of the size of STARFLUX and OBJFLUX.
   ; Pad to larger (or equal) 2^N, and then doubled for isolated b.c.

   npixstar = (size(starflux))[1]
   npixbig = 2L^(fix(alog(npixstar > npixobj)/alog(2) + 1.9999))

   ;---------------------------------------------------------------------------
   ; Compute FFT for stellar template

   nstar = n_elements(starflux)/(size(starflux))[1]

   FOR istar=0, nstar-1 DO BEGIN 


   veldisp_fft, starflux[*, istar], starivar[*, istar], npixbig, starfft,  $
     starfilt, starvar0, starvariancefft, starivar_pad, $
     klo_cut=klo_cut, khi_cut=khi_cut, wave=starwave, keep=keep

   fitredshift, starfilt, starivar_pad, starfilt, starivar_pad, $
      nsearch=5, zfit=starcen, z_err=starcen_err, $
      veldispfit=starsigma, veldisp_err=starsigma_err, doplot=doplot

   ;---------------------------------------------------------------------------
   ; LOOP OVER OBJECT SPECTRA

   print, '     plt fiber     redshift          veldisp_cc^2     veldisp_q^2'

   FOR iobj=0, nobj-1 DO BEGIN 


      fluxivar = objivar[*, iobj]
      veldisp_fft, objflux[*,iobj], fluxivar, npixbig,  $
        fluxfft, fluxfilt, fluxvar0, fluxvariancefft, fluxivar_pad,  $
        klo_cut=klo_cut, khi_cut=khi_cut

      fitredshift, fluxfilt, fluxivar_pad, starfilt, starivar_pad, $
        nsearch=5, zfit=fitcen, z_err=fitcen_err, $
        veldispfit=galsigma, veldisp_err=galsigma_err, zconf=zconf, $
        doplot=doplot
      
; 2nd try

      IF keyword_set(keep) THEN BEGIN 
          veldisp_fft, objflux[*,iobj], fluxivar, npixbig,  $
            fluxfft, fluxfilt, fluxvar0, fluxvariancefft, fluxivar_pad, $
            keep=keep*10.^(fitcen/10000.)
          
;          fitredshift, fluxfilt, fluxivar_pad, starfilt, starivar_pad, $
;            nsearch=5, zfit=fitcen, z_err=fitcen_err, $
;            veldispfit=galsigma, veldisp_err=galsigma_err, doplot=doplot
      ENDIF 

      result[iobj].z[istar]     = 10.^(fitcen/10000)-1. ; dimensionless z
      result[iobj].z_err[istar] = alog(10)*1e-4*fitcen_err* $
        (1+result[iobj].z[istar])
      result[iobj].zconf[istar] = zconf

      if (keyword_set(doplot)) then begin
         wset,2
         djs_plot, starflux/normstar, xr=[0,2000], yr=[0,2.0], $
           title='Rest frame spectra of template (white) and galaxy (red)'
         djs_oplot, $
          smooth(shift(objflux[*,iobj]/normobj, -fitcen),5), $
           color='red'
      endif

; Should really store sigma squared, and allow negative values; error
; should reflect it - DPF ???
      if (galsigma GT starsigma AND starsigma GT 0.0) then begin
         result[iobj].sigma2_cc[istar] = galsigma^2 - starsigma^2

; fix this
         result[iobj].sigma2_cc_err[istar] = sqrt((galsigma*galsigma_err)^2 + $
          (starsigma*starsigma_err)^2)/result[iobj].sigma2_cc[istar]
      endif

      twopiei = 2.0 * PI * complex(0.0,1.0)
      knums = fft_wavenums(npixbig)
      phase = exp( - twopiei * knums * fitcen)
      starshift = starfft * phase

      testsigma = findgen(ceil(float(maxsig)/sigmastep) + 1) * sigmastep

; Need to pick lower and upper limits to do comparison
; Let's try to compare from 80 pixels to 2.2 pixels

      if (NOT keyword_set(nodiff)) then $
       answer = fourier_difference(fluxfft, starshift, fluxvar0, $
                starvar0, testsigma=testsigma, deltachisq=1.0, $
                lowlimit = 1.0/80.0, highlimit=1.0/2.2, doplot=doplot)

      bestalpha = -9999.0

      if (n_elements(answer) EQ 4) then begin
         result[iobj].sigma2_diff[istar] = answer[1]
         result[iobj].sigma2_diff_err[istar] = answer[2]
         bestalpha = answer[3]
      endif

      if (keyword_set(doplot)) then quoplot = 2 else quoplot = 0
      IF NOT keyword_set(noquotient) THEN BEGIN 
         answerq = fourier_quotient(fluxfft, starshift, fluxvar0, $
                  starvar0, testsigma2=testsigma^2, deltachisq=1.0, $
                  lowlimit = 1.0/250.0, highlimit=1.0/5., doplot=quoplot, $
                                     broadarr=broadarr)
      ENDIF 

      if (n_elements(answerq) EQ 4) then begin
         result[iobj].sigma2_quotient[istar] = answerq[1]
         result[iobj].sigma2_quotient_err[istar]  = answerq[2]
         bestalpha_q = answerq[3]
      endif

      IF NOT keyword_set(noquotient) THEN BEGIN
          r = result[iobj]
          print, iobj, r.plate, r.fiber, r.z[istar], r.z_err[istar], r.sigma2_cc[istar], $
            r.sigma2_cc_err[istar], r.sigma2_quotient[istar], r.sigma2_quotient_err[istar], $
            format='(i4,i5,i4,f8.5," +/-",f8.5,2(f8.3," +-",f6.3))'

      ENDIF ELSE BEGIN 
          print, iobj, result[iobj].z[istar], result[iobj].z_err[istar],  $
            result[iobj].sigma2_cc[istar], result[iobj].sigma2_cc_err[istar], $
            format='(i4,f9.3,3(f8.3))'
      ENDELSE 

   endfor
   endfor
   return
end
