;+
; NAME:
;   mfourier_difference
;
; PURPOSE:
;   Perform a chi2 fit to the fourier difference between a single
;   galaxy and a broadened stellar template to calculate velocity dispersion
;   and uncertainty on velocity dispersion
;
; CALLING SEQUENCE:
;   answers = mfourier_difference(galfft, starfft, galvar0, starvar0, $
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
;   2000-Sep-11  Considerable changes by the SWAT team
;-
;------------------------------------------------------------------------------
function alpha_chisq, alpha, gal, star, galvar, starvar, br, minchi2

   ntry = n_elements(alpha)
   if (ntry EQ 0) then return, -1

   chisq = fltarr(ntry)
   for i=0,ntry -1 do begin
      diff = float(gal) - float(star)*br*alpha[i]
      rediff = float(diff*conj(diff))
;      rediff = abs(diff)^2
      chisq[i] = total(rediff/(galvar + starvar*alpha[i]^2*br^2))
   endfor 

   mfindchi2min, alpha, chisq, minchi2, bestalpha

   return, bestalpha
end


;------------------------------------------------------------------------------
function get_alpha, galvar0norm, starvar0norm, broad, galfftnorm, $
        starfftnorm, maxiter=maxiter, chi2=chi2
      alpha_old = 1
      num1 = float(galfftnorm*conj(starfftnorm)*broad)
      num2 = abs(starfftnorm*broad)^2   
      for j=0,maxiter-1 do begin 
        denom  = galvar0norm + (starvar0norm*alpha_old^2)*broad^2
        alpha  = total(num1/denom)/total(num2/denom)
        alpha_old = alpha
;	print,j,alpha
      endfor

return, alpha
end


;------------------------------------------------------------------------------
; Suggestion on errors: 
;  The error on fourier modes can be estimated from the high k modes. 
;  Correlations between pixels leads to rolloff in the errors at high k,
;  compensated by increased errors at low k.  Suggest that 2-d output an
;  approximate sub-diagonal correlation term.  (Eisenstein)
;  Chi^2 can be calibrated by taking a random wavelength-space vector of 
;  unit normals and passing through cosbell (apodized) FFT, estimating
;  Fourier mode variance in same way as galaxy, computing chi^2 and
;  comparing result with analytic result. 

; Could also use realization of quoted errors from 2-D for galaxy spectra. 
; 
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
function newdiff, galfft, starfft, galvar0, starvar0, $
 testsigma=testsigma, lowlimit = lowlimit, highlimit=highlimit, $
 deltachisq=deltachisq, doplot=doplot, broadarr=broadarr

   if (NOT keyword_set(lowlimit)) then lowlimit = 1.0/80.0
   if (NOT keyword_set(highlimit)) then highlimit = 1.0/2.2

   if (size(galfft, /tname) EQ 'DOUBLE') then PI = !dpi $
    else PI = !pi

   knums = fft_wavenums(N_elements(galfft))

   inside = where(knums GT lowlimit AND $
                  knums LT highlimit, ninside)

   if (inside[0] EQ -1) then begin
      print, 'No pixels in correct frequency range'
      return, -1
   endif

   if (n_elements(testsigma) EQ 0) then testsigma = findgen(30)*0.2

   nloop = n_elements(testsigma)
   chi2diff = fltarr(nloop)
   sigma = fltarr(nloop)
   alpha = fltarr(nloop)

   galnorm = djsig(float(galfft[inside]), sigrej=5)
   starnorm =  djsig(float(starfft[inside]), sigrej=5)

   galfftnorm= float(galfft[inside]/galnorm)
   starfftnorm= float(starfft[inside]/starnorm)
   galvar0norm= galvar0/galnorm^2
   starvar0norm = starvar0/starnorm^2

   ones = 1.+fltarr(ninside)
   IF NOT keyword_set(broadarr) THEN BEGIN 
          broadarr = dblarr(ninside, nloop)
          for i=0,nloop-1 do begin
             IF testsigma[i] EQ 0 THEN broad = ones ELSE BEGIN 
                fsig = 1.d/(2.*!dpi)/testsigma[i]
                broad = gauss_periodic(knums[inside], [1., 0., fsig], shft=1.)
             ENDELSE 
             broadarr[*, i] = broad
          ENDFOR     
   ENDIF 


   for i=0,nloop-1 do begin
      broad = broadarr[*, i]
      alpha[i] = get_alpha(galvar0norm, starvar0norm, broad, $
        galfftnorm, starfftnorm, maxiter=4)

      denom = galvar0norm + starvar0norm*(alpha[i]*broad)^2
      chi2diff[i] = total(abs(galfftnorm-starfftnorm*(broad*alpha[i]))^2 / $
	                    denom)
   endfor

   deltachisq = 1./n_elements(inside)

   mfindchi2min, testsigma, chi2diff, minchi2, minsigma, errsigma, $
    deltachisq = deltachisq, doplot=doplot, npts=ninside

   broad = exp(-(knums[inside]*minsigma * 2.0 * PI)^2/2.0)
   bestalpha = get_alpha(galvar0norm, starvar0norm, broad, $
     galfftnorm, starfftnorm, maxiter=4)

;   oplot, testsigma, alpha, ps=2
   minc = min(chi2diff, alphaplace)
   bestalpha =alpha[alphaplace]


   return, [minchi2, minsigma, errsigma, bestalpha]
end
;------------------------------------------------------------------------------
