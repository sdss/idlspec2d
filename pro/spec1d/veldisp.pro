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
;    klo_cut=, khi_cut=, maxsig=, sigmastep=, zoffset= ]
;
; INPUTS:
;   objflux    - Array of object spectra [NPIX] or [NPIX,NOBJ]
;   objivar    - Array of object inverse variance [NPIX] or [NPIX,NOBJ]
;   starflux   - Template spectrum [nstarpix]
;   starivar   - Template inverse variance [nstarpix]
;
; OPTIONAL KEYWORDS:
;   klo_cut    - Low frequency cutoff for cross-correlation peak finding;
;                default to 1/128.
;   khi_cut    - High frequency cutoff for cross-correlation peak finding;
;                default to 1/3.
;   maxsig     - Maximum velocity dispersion to search for; default to 2 pix
;   sigmastep  - Steps between each test sigma; default to 0.2 pix
;   zoffset    - ???
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
function veldisp_struc, nn

   struc = { $
            z                   : 0., $
            z_err               : 0., $
            zconf               : 0., $
            sigma2_cc           : 0., $
            sigma2_cc_err       : 0., $
            sigma2_quotient     : 0., $
            sigma2_quotient_err : 0. }

   if (keyword_set(nn)) then return, replicate(struc, nn) $
    else return, struc
end

;------------------------------------------------------------------------------
; this routine:
;   1) interpolates accross missing pixels
;   2) normalizes spectrum to unit flux (or variance)
;   3) apodizes spectrum with cos bell in real space
;   4) pads with zeros up to TWICE the next higher value of 2^N
;   5) computes FFT

pro veldisp_fft, flux_in, err_in, npixbig, fluxfft, fluxfilt, fluxvar0, $
 err, klo_cut=klo_cut, khi_cut=khi_cut

   flux = flux_in-mean(flux_in)
   err = err_in

   ; Interpolate over bad regions
   flux = djs_maskinterp(flux, err LE 0.0, /const)
  
   ; Apodize
   fft_apodize, flux, err
  
   ; Pad
   npad = npixbig - n_elements(flux)
   if (npad GT 0) then begin
      flux = [flux, fltarr(npad)]
      err = [err, fltarr(npad)]
   endif

   ; Take FFT
   fluxfft = fft(flux) * npixbig
   fluxvariancefft = fft(err^2)  * npixbig
   fluxvar0 = float(fluxvariancefft[0])
 
   ; Band-pass filter the object spectrum
   w = where(err NE 0)
   fluxfilt = bandpassfilter(fluxfft, klo_cut=klo_cut, khi_cut=khi_cut)
   norm = stddev((float(fft(fluxfilt, /inv)))[w],/double)

   fluxfilt = fluxfilt / norm

   if (total(finite(fluxfilt) EQ 0) NE 0) then $
    message, 'Infinite value in FFT'

   return
end

;------------------------------------------------------------------------------
function veldisp, objflux, objivar, starflux, starivar, zoffset=zoffset, $
 klo_cut=klo_cut, khi_cut=khi_cut, maxsig=maxsig, sigmastep=sigmastep
       
   if (NOT keyword_set(klo_cut)) then klo_cut = 1./128.
   if (NOT keyword_set(khi_cut)) then khi_cut = 1./3.
   if (NOT keyword_set(maxsig)) then maxsig = 2.0
   if (NOT keyword_set(sigmastep)) then sigmastep = 0.2

   if (size(objflux, /tname) EQ 'DOUBLE') then PI = !dpi $
    else PI = !pi

   ;---------------------------------------------------------------------------
   ; Check dimensions of input arrays

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npixobj = dims[0]
   if (ndim EQ 1) then nobj = 1 $
    else if (ndim EQ 2) then nobj = dims[1] $
    else message, 'OBJFLUX is neither 1-D or 2-D'

   ndim = size(starflux, /n_dimen)
   dims = size(starflux, /dimens)
   npixstar = dims[0]
   if (ndim EQ 1) then nstar = 1 $
    else if (ndim EQ 2) then nstar = dims[1] $
    else message, 'STARFLUX is neither 1-D or 2-D'

   if total(abs(size(starflux, /dimens)-size(starivar, /dimens))) NE 0 $
    OR size(starflux, /n_dimen) NE size(starivar, /n_dimen) THEN  $
    message, 'Dimensions of STARFLUX and STARIVAR do not match'

   if total(abs(size(objflux, /dimens)-size(objivar, /dimens))) NE 0 $
    OR size(objflux, /n_dimen) NE size(objivar, /n_dimen) THEN  $
    message, 'Dimensions of OBJFLUX and OBJIVAR do not match'

   ;---------------------------------------------------------------------------
   ; Decide how large the padded spectra should be, based upon the
   ; larger of the size of STARFLUX and OBJFLUX.
   ; Pad to larger (or equal) 2^N, and then doubled for isolated B.C.

   npixbig = 2L^(fix(alog(npixstar > npixobj)/alog(2) + 1.9999))

   result = veldisp_struc(nobj)

   ;---------------------------------------------------------------------------
   ; Compute FFT for stellar template

   for istar=0, nstar-1 do begin

      veldisp_fft, starflux[*,istar], starivar[*,istar], npixbig, starfft, $
       starfilt, starvar0, starivar_pad, $
       klo_cut=klo_cut, khi_cut=khi_cut

      fitredshift, starfilt, starivar_pad, starfilt, starivar_pad, $
       nsearch=5, zfit=starcen, z_err=starcen_err, $
       veldispfit=starsigma, veldisp_err=starsigma_err

   ;---------------------------------------------------------------------------
   ; LOOP OVER OBJECT SPECTRA

   splog, 'Object#    z [pixels]        veldisp_cc^2     veldisp_q^2'

   for iobj=0, nobj-1 do begin

      veldisp_fft, objflux[*,iobj], objivar[*,iobj], npixbig,  $
        fluxfft, fluxfilt, fluxvar0, fluxivar_pad,  $
        klo_cut=klo_cut, khi_cut=khi_cut

      fitredshift, fluxfilt, fluxivar_pad, starfilt, starivar_pad, $
        nsearch=5, zfit=fitcen, z_err=fitcen_err, $
        veldispfit=galsigma, veldisp_err=galsigma_err, zconf=zconf

      result[iobj].z[istar] = fitcen
      if (keyword_set(zoffset)) then $
       result[iobj].z[istar] = result[iobj].z[istar] + zoffset[istar]
      result[iobj].z_err[istar] = fitcen_err
      result[iobj].zconf[istar] = zconf

;      if (keyword_set(doplot)) then begin
;         wset,2
;         djs_plot, starflux/normstar, xr=[0,2000], yr=[0,2.0], $
;           title='Rest frame spectra of template (white) and galaxy (red)'
;         djs_oplot, $
;          smooth(shift(objflux[*,iobj]/normobj, -fitcen),5), $
;           color='red'
;      endif

; Should really store sigma squared, and allow negative values; error
; should reflect it - DPF ???
      if (galsigma GT starsigma AND starsigma GT 0.0) then begin
         result[iobj].sigma2_cc[istar] = galsigma^2 - starsigma^2

; fix this ???
         result[iobj].sigma2_cc_err[istar] = sqrt((galsigma*galsigma_err)^2 + $
          (starsigma*starsigma_err)^2)/result[iobj].sigma2_cc[istar]
      endif

      twopiei = 2.0 * PI * complex(0.0,1.0)
      knums = fft_wavenums(npixbig)
      phase = exp( - twopiei * knums * fitcen)
      starshift = starfft * phase

      testsigma = findgen(ceil(float(maxsig)/sigmastep) + 1) * sigmastep

      answerq = fourier_quotient(fluxfft, starshift, fluxvar0, $
       starvar0, testsigma2=testsigma^2, deltachisq=1.0, $
       lowlimit = 1.0/250.0, highlimit=1.0/5., broadarr=broadarr)

      if (n_elements(answerq) EQ 4) then begin
         result[iobj].sigma2_quotient[istar] = answerq[1]
         result[iobj].sigma2_quotient_err[istar]  = answerq[2]
         bestalpha_q = answerq[3]
      endif

      r = result[iobj]
      splog, iobj, r.z[istar], r.z_err[istar], $
       r.sigma2_cc[istar], r.sigma2_cc_err[istar], $
       r.sigma2_quotient[istar], r.sigma2_quotient_err[istar], $
       format='(i5,f9.1," +/-",f7.1,2(f9.3," +-",f6.3))'

   endfor
   endfor

   return, result
end
