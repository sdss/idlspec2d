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
 testsigma2=testsigma2, lowlimit = lowlimit, highlimit=highlimit, $
 deltachisq=deltachisq, doplot=doplot, broadarr=broadarr

   if (NOT keyword_set(lowlimit)) then lowlimit = 1.0/80.0
   if (NOT keyword_set(highlimit)) then highlimit = 1.0/5.

   if (size(galfft, /tname) EQ 'DOUBLE') then PI = !dpi $
    else PI = !pi

   knums = fft_wavenums(n_elements(galfft))
   inside = where(abs(knums) GT lowlimit AND $
                  abs(knums) LT highlimit, ninside)

   if (inside[0] EQ -1) then begin
      print, 'No pixels in correct frequency range'
      return, -1
   endif

   if (n_elements(testsigma2) EQ 0) then testsigma2 = findgen(30)*0.2

;       Tonry and Davis show a simple expression to maximize
;       This is a slow minimizer, stepping through 30 sigmas to
;       find best one.  This method is the fft difference method
;       We need a routine for each method.

      nloop = n_elements(testsigma2)
      chi2 = fltarr(nloop)
      sigma = fltarr(nloop)
      alpha = fltarr(nloop)

      alphatry = findgen(21)*0.1 

      q = galfft[inside]/starfft[inside]

; Reject outliers on q (they can be quite large and drive the fit)
      qs = exp(smooth(alog(q), 25, /edge))
      dif = float(alog((q/qs)))
      djs_iterstat, dif, sigma=qsig
      wbad = where(abs(dif) GT qsig*5, ct)
      IF ct GT 0 THEN q[wbad] = qs[wbad]

; Now reevaluate smoothed q
      qs = exp(smooth(alog(q), 9, /edge))
qs = median(q, 75)

      var = (galvar0/ float(galfft[inside]*conj(galfft[inside])) + $
             starvar0 / float(starfft[inside]*conj(starfft[inside])) * $
             float(q)^2)

      ones = 1.+fltarr(n_elements(q))
      IF NOT keyword_set(broadarr) THEN BEGIN 
          broadarr = dblarr(n_elements(inside), nloop)
          for i=0,nloop-1 do begin
             IF testsigma2[i] EQ 0 THEN broad = ones ELSE BEGIN 
                fsig = 1.d/(2.*!dpi)/sqrt(abs(testsigma2[i]))
                broad = gauss_periodic(knums[inside], [1., 0., fsig], shft=1.)
                IF testsigma2[i] LT 0 THEN broad = 1./broad
             ENDELSE 
             broadarr[*, i] = broad
          ENDFOR     


      ENDIF 


      for i=0,nloop-1 do begin
          broad = broadarr[*, i]

;          alpha[i] = total(float(q) * broad / var)/total(broad^2/var)
          alpha[i] = total(float(qs) * broad)/total(broad^2)

;          qres = alog(float(qs/(alpha[i]*broad)))
          qres = float(qs-alpha[i]*broad)

;          qresbar = total(qres/sqrt(var))/total(1./sqrt(var))
;          chi2[i] = total((qres-qresbar)^2/var)

          chi2[i] = total(qres^2)

;          plot, knums[inside], qs, ps=3, yr=[-1, 2]
;          oplot, knums[inside], broad*alpha[i], ps=3
;          print, testsigma[i], chi2[i]
      endfor
plot, [1], [1], xr=[0, 1.5], yr=[-1, 5]
      findchi2min, testsigma2, chi2, minchi2, minsigma, errsigma, $
	  deltachisq = deltachisq, doplot=doplot, npts= ninside

      bestalpha = (interpol(alpha, testsigma2, minsigma))[0]
      fsig= 1.d/(2.*!dpi)/sqrt(minsigma)
      broad = gauss_periodic(knums[inside], [1., 0., fsig], shft=1.)

      plot, knums[inside], qs, ps=3, yr=[-1, 2]
      oplot, knums[inside], broad*bestalpha, ps=3


      return, [minchi2, minsigma, errsigma, bestalpha]
end 
