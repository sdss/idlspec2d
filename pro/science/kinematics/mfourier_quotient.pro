;+
; NAME:
;   mfourier_quotient
;
; PURPOSE:
;   Perform a chi2 fit to the fourier quotient of a single
;    galaxy and a broadened stellar template to calculate velocity dispersion
;    and uncertainty on velocity dispersion
;
; CALLING SEQUENCE:
;   answers = mfourier_quotient(galfft, starfft, galvar0, starvar0, $
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
;	Same inputs and outputs as mfourier_difference
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
;   27-Jun-2000  Completely rewritten - Finkbeiner & SWAT team
;-
;------------------------------------------------------------------------------
function mfourier_quotient, galfft, starfft, galvar0, starvar0, $
 testsigma2=testsigma2, lowlimit = lowlimit, highlimit=highlimit, $
 deltachisq=deltachisq, doplot=doplot, broadarr=broadarr

; testsigma2 is trial velocity dispersion^2

   if (NOT keyword_set(lowlimit)) then lowlimit = 1.0/80.0
   if (NOT keyword_set(highlimit)) then highlimit = 1.0/5.


   knums = fft_wavenums(n_elements(galfft))
   inside = where(abs(knums) GT lowlimit AND $
                  abs(knums) LT highlimit, ninside)

   if (inside[0] EQ -1) then begin
      print, 'No pixels in correct frequency range'
      return, -1
   endif

   if (n_elements(testsigma2) EQ 0) then testsigma2 = findgen(30)*0.2
   
   nloop = n_elements(testsigma2)
   chi2 = fltarr(nloop)
   sigma = fltarr(nloop)
   alpha = fltarr(nloop)

   alphatry = findgen(21)*0.1 

   galnorm = djsig(float(galfft[inside]), sigrej=5)
   starnorm =  djsig(float(starfft[inside]), sigrej=5)

   q = galfft[inside]/starfft[inside]*starnorm/galnorm

; Reject outliers on q (they can be quite large and drive the fit)
      qs = exp(smooth(alog(q), 25, /edge))
      dif = float(alog((q/qs)))
      djs_iterstat, dif, sigma=qsig, sigrej=5
      wbad = where(abs(dif) GT qsig*5, ct)
      IF ct GT 0 THEN q[wbad] = qs[wbad]

;	var = galvar0/ float(galfft*conj(galfft)) + $
;             starvar0 / float(starfft*conj(starfft))

;	 var1=  var[inside] * float(qs)^2
;      var= (median(var,150))[inside] * float(qs)^2

;      var = (galvar0/ float(galfft[inside]*conj(galfft[inside])) + $
;             starvar0 / float(starfft[inside]*conj(starfft[inside])) * $
;              float(q)^2)      
;             float(qs)^2)

;      var= exp(smooth(alog(var), 25, /edge))
;       var= (median(var,150))

;      dif = (float(q) - float(qs))^2/var

      

; Now reevaluate smoothed q

        qs = median(q, 75)

      var = (galvar0/ float(galfft[inside]*conj(galfft[inside])) + $
             starvar0/ float(starfft[inside]*conj(starfft[inside])) * $
;              float(q)^2)      
             float(qs)^2)


; Define broadarr (array of broadening functions) if not already defined
      ones = 1.+fltarr(n_elements(q))
      IF NOT keyword_set(broadarr) THEN BEGIN 
          broadarr = dblarr(n_elements(inside), nloop)
          for i=0,nloop-1 do begin
             IF testsigma2[i] EQ 0 THEN broad = ones ELSE BEGIN 
;                fsig = 1.d/(2.*!dpi)/sqrt(abs(testsigma2[i]))
	         fsig = 1.d/(2.*!dpi)/testsigma2[i]
                broad = gauss_periodic(knums[inside], [1., 0., fsig], shft=1.)
                IF testsigma2[i] LT 0 THEN broad = 1./broad
             ENDELSE 
             broadarr[*, i] = broad
          ENDFOR     


      ENDIF 

      for i=0,nloop-1 do begin
          broad = broadarr[*, i]
;          broad = broad/max(broad)
;           broad = broad*max(smooth(alog(qs), 25, /edge))
              
          alpha[i] = total(float(qs) * broad / var)/total(broad^2/var)

;          alpha[i] = total(float(qs) * broad)/total(broad^2)
; nul residual (for just noise)
          nul = sqrt(alpha[i]^2*broad^2+1)
;          res = broad*alpha[i]*starfft[inside]/starnorm-galfft[inside]/galnorm
 	  res= float(qs)-alpha[i]*broad

;IF i EQ 50 THEN 
;          chi2[i] = total(abs(res)^2)/n_elements(inside)/mean(nul)
           chi2[i] = total((res)^2/var)  


;         window, 1,retain=2   
;         broad= broad*alpha[i]
;         plot,knums[inside],q,ps=0,xr=[0,0.5],yr=[-1,2]
;         djs_oplot,knums[inside],qs,ps=0,color='green'
;         djs_oplot,knums[inside],broad,color='red'
;stop

      endfor

; Now do loop of negatives

; now plot

       varscale = (galvar0/ float(galfft[inside]*conj(galfft[inside])) + $
             starvar0/ float(starfft[inside]*conj(starfft[inside])) * $
              float(q)^2)

	 resscale= float(q-qs)
         chi2scale = total((resscale)^2/varscale)  

      deltachisq = 1./n_elements(inside)

      mfindchi2min, testsigma2, chi2, minchi2, minsigma, errsigma, $
	  deltachisq = deltachisq, doplot=doplot, npts= ninside, $
         chi2scale = chi2scale

;	print,'errsigma',errsigma

      bestalpha = (interpol(alpha, testsigma2, minsigma))[0]
       fsig= 1.d/(2.*!dpi)/minsigma
      broad = gauss_periodic(knums[inside], [1., 0., fsig], shft=1.)
;      res = broad*starfft[inside]-galfft[inside]
       res= float(q-bestalpha*broad)
	
       chi2 = total((res)^2/varscale) 

;	if (chi2scale NE 0) then $
;        errsigma = sqrt(chi2scale * deltachisq) * rms $
;         else $
;	errsigma = sqrt(minchi2 * deltachisq) * rms

;        uppersigma =  minsigma + errsigma
;        lowersigma = minsigma - errsigma 
;        errplot, lowersigma,-10, +10
;	errplot, uppersigma, -10, +10
;stop
;         window, 1,retain=2   
;         broad= broad*max(float(smooth(qs,40,/edge)))
; 	broad= broad*bestalpha
;         plot,knums[inside],qs,ps=0,xr=[0,0.5],yr=[-1,2]
;         djs_oplot,knums[inside],broad,color='red'


;      errsigma = errsigma*sqrt(minchi2)
      
;      IF keyword_set(doplot) THEN BEGIN 
;          plot, knums[inside], qs, ps=3, yr=[-1, 2]
;          oplot, knums[inside], broad*bestalpha, ps=3
;          plot, testsigma2, chi2, ps=-7
;          oplot,[0,100],[1,1]*minchi2+deltachisq
;     ENDIF 

      return, [minchi2, minsigma, errsigma, bestalpha]
end 
