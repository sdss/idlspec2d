;+
; NAME:
;   fourier_quotient
;
; PURPOSE:
;   Perform a chi2 fit to the fourier quotient of a single
;    galaxy and a broadened stellar template to calculate velocity dispersion
;    and uncertainty on velocity dispersion
;
; CALLING SEQUENCE:
;   answers = fourier_quotient(galfft, starfft, galvar0, starvar0, $
;    testsigma=, lowlimit=, highlimit=, $
;    deltachisq=, /doplot)
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
;	Same inputs and outputs as fourier_difference
;
; EXAMPLES:
;
; BUGS:
;
;	Need to ensure that confidence level returned as errsigma is proper
;
; PROCEDURES CALLED:
;   
;
; REVISION HISTORY:
;   25-Mar-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function fourier_quotient, galfft, starfft, galvar0, starvar0, $
 testsigma=testsigma, lowlimit = lowlimit, highlimit=highlimit, $
 deltachisq=deltachisq, doplot=doplot

   if (NOT keyword_set(lowlimit)) then lowlimit = 1.0/80.0
   if (NOT keyword_set(highlimit)) then highlimit = 1.0/2.2

   if (size(galfft, /tname) EQ 'DOUBLE') then PI = !dpi $
    else PI = !pi

   knums = fft_wavenums(n_elements(galfft))
   inside = where(abs(knums) GT lowlimit AND $
                  abs(knums) LT highlimit, ninside)

   if (inside[0] EQ -1) then begin
      print, 'No pixels in correct frequency range'
      return, -1
   endif

   if (n_elements(testsigma) EQ 0) then testsigma = findgen(30)*0.2

;       Tonry and Davis show a simple expression to maximize
;       This is a slow minimizer, stepping through 30 sigmas to
;       find best one.  This method is the fft difference method
;       We need a routine for each method.

      nloop = n_elements(testsigma)
      chi2 = fltarr(nloop)
      sigma = fltarr(nloop)
      alpha = fltarr(nloop)

      alphatry = findgen(21)*0.1 

      q = float(galfft[inside]/starfft[inside])
      var = (galvar0/ float(galfft[inside]*conj(galfft[inside])) + $
              starvar0 / float(starfft[inside]*conj(starfft[inside])) * q^2)

      for i=0,nloop-1 do begin
          IF testsigma[i] EQ 0 THEN broad = 1. ELSE BEGIN 
              fsig = 1.d/(2.*!dpi)/testsigma[i]
              broad = gauss_periodic(knums[inside], [1., 0., fsig], shft=1.)
          ENDELSE 

;        broad = exp(-(knums[inside]*testsigma[i] * 2.0 * PI)^2/2.0)

          alpha[i] = total(q * broad / var)/total(broad^2/var)
          qres = (q-alpha[i]*broad)
;          qresbar = total(qres/sqrt(var))/total(1./sqrt(var))
          qresbar = 0
          chi2[i] = total((qres-qresbar)^2/var)

      endfor

      findchi2min, testsigma, chi2, minchi2, minsigma, errsigma, $
	  deltachisq = deltachisq, doplot=doplot, npts= ninside

   ;   oplot, testsigma, alpha, ps=2

;      minc = min(chi2, alphaplace)
;      bestalpha =alpha[alphaplace]
      bestalpha = interpol(alpha, testsigma, minsigma)
;fsig= 1.d/(2.*!dpi)/minsigma
;print,fsig
    broad = gauss_periodic(knums[inside], [1., 0., fsig], shft=1.)
;stop
      return, [minchi2, minsigma, errsigma, bestalpha]
end 
