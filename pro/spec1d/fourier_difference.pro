
;+
; NAME:
;   fourier_difference
;
; PURPOSE:
;   Perform a chi2 fit to the fourier difference between a single
;   galaxy and a broadened stellar template to calculate velocity dispersion
;   and uncertainty on velocity dispersion
;
; CALLING SEQUENCE:
;   answers = fourier_difference(galfft, starfft, galvar0, starvar0, $
;          testsigma=, lowlimit=, highlimit=, $
;          deltachisq=deltachisq, /doplot)
;
; INPUTS:
;   galfft     - Fourier transform of galaxy
;   starfft    - Fourier transform of stellar template
;   galvar0    - error in galaxy fft (0th element of galaxy error FFT)
;   starvar0   - error in stellar fft (0th element of stellar error FFT)
;
; OPTIONAL KEYWORDS:
;   testsigma  - Array of sigma values to calculate chi2
;   lowlimit   - lower boundary of chi2 sum (in knums units)
;   highlimit  - upper boundary of chi2 sum (in knums units)
;   deltachisq - chi2 difference from minimum to set error on velocity dispersion
;   doplot     - Output plots to xwindow
;
; OUTPUTS:
;   answers    - Four element array with:
;                [minchi2, minsigma, errsigma, bestalpha]
;                bestalpha is the normalization constant between galaxy and star
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Currently, this is very slow, as we have to check 11 normalizations
;   for each element in testsigma array
;
;   Need to ensure that confidence level returned as errsigma is proper
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   25-Mar-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function alpha_chisq, alpha, gal, star, galvar, starvar, br, minchi2

   ntry = n_elements(alpha)
   if (ntry EQ 0) then return, -1

   chisq = fltarr(ntry)
   for i=0,ntry -1 do begin
      diff = gal - star*br*alpha[i]
      rediff = float(diff*conj(diff))
      chisq[i] = total(rediff/(galvar + starvar*alpha[i]^2*br^2))
   endfor

   findchi2min, alpha, chisq, minchi2, bestalpha
  
   return, bestalpha
end

;------------------------------------------------------------------------------
function fourier_difference, galfft, starfft, galvar0, starvar0, $
 testsigma=testsigma, lowlimit = lowlimit, highlimit=highlimit, $
 deltachisq=deltachisq, doplot=doplot

   if (NOT keyword_set(lowlimit)) then lowlimit = 1.0/80.0
   if (NOT keyword_set(highlimit)) then highlimit = 1.0/2.2

   if (size(galfft, /tname) EQ 'DOUBLE') then PI = !dpi $
    else PI = !pi

   knums = fft_wavenums(N_elements(galfft))

   inside = where(abs(knums) GT lowlimit AND $
                  abs(knums) LT highlimit, ninside)

   if (inside[0] EQ -1) then begin
      print, 'No pixels in correct frequency range'
      return, -1
   endif

   if (n_elements(testsigma) EQ 0) then testsigma = findgen(30)*0.2

   ; This is a slow minimizer, stepping through 30 sigmas to
   ; find best one.  This method is the fft difference method

   nloop = n_elements(testsigma)
   chi2diff = fltarr(nloop)
   sigma = fltarr(nloop)
   alpha = fltarr(nloop)

   alphatry = findgen(11)*0.2 
     
   for i=0,nloop-1 do begin
      broad = exp(-(knums*testsigma[i] * 2.0 * PI)^2/2.0)

      alpha[i] = alpha_chisq(alphatry, galfft[inside], starfft[inside], $
       galvar0, starvar0, broad[inside], bestchisq)
      chi2diff[i] = bestchisq

;     answer = amoeba(1.0e-2, function_name='alpha_chisq',  $
;      p0 = [1.1], scale=[0.1], nmax = 20, function_value=areturn)

   endfor

   findchi2min, testsigma, chi2diff, minchi2, minsigma, errsigma, $
    deltachisq = deltachisq, doplot=doplot, npts=ninside

;   oplot, testsigma, alpha, ps=2
   minc = min(chi2diff, alphaplace)
   bestalpha =alpha[alphaplace]
   return, [minchi2, minsigma, errsigma, bestalpha]
end
;------------------------------------------------------------------------------
