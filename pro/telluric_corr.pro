;+
; NAME:
;   telluric_corr
;
; PURPOSE:
;   Use SPECTROPHOTO and REDDEN_STD's to fit telluric features
;   between wavelengths 10^minw and 10^maxw
;
; CALLING SEQUENCE:
;   telluric_factor = telluric_corr(flux, fluxivar, wset, plugsort, $
;       contwave=contwave, contflux=contflux, contivar=contivar,    $
;       telluricbkpt=telluricbkpt, telluriccoeff=telluriccoeff, $
;       minw=minw, maxw=maxw, lower=lower, upper=upper, $
;       ncontbkpts = ncontbkpts, fibermask=fibermask)
;
; INPUTS:
;   flux         - sky-subtracted extracted spectra [nx,ntrace]
;   fluxivar     - corresponding inverse variance [nx,ntrace]
;   wset         - wavelength coefficients as a trace set
;   plugsort     - plugmap entries
;
; OPTIONAL KEYWORDS:
;   minw         - minimum wavelength to fit (Default 3.82)
;   maxw         - maximum wavelength to fit (Default 3.92)
;   ncontbkpts   - Number of bkpts to fit continuum in telluric region (5)
;   lower        - lower rejection threshold for telluric fitting
;   upper        - upper rejection threshold for telluric fitting
;   fibermask    - use to reject possible standards which have spectral troubles
;   
; OUTPUTS:
;   telluric_factor - Telluric correction for each pixel in flux array
;
; OPTIONAL OUTPUTS:
;   contflux     - Normalized flux of stars used in telluric fitting
;   contivar     - Corresponding inverse variance
;   contwave     - Corresponding wavelengths  (log lambda)
;   telluricbkpt - bkpts used in telluric absorption fit
;   telluriccoeff- best fit b-spline coefficients for telluric absorption
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   There is still some low order residual from flux correction
;     (I think due to telluric absorption), which is in turned
;   removed by the telluric_correction.  Although there is not
;   a clean break between the two steps, used together they seem
;   to correct the spectra properly
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   19-Oct-1999  Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
function telluric_corr,flux, fluxivar, wset, plugsort, $
        contwave=contwave, contflux=contflux, contivar=contivar,    $
        telluricbkpt=telluricbkpt, telluriccoeff=telluriccoeff, $
        minw=minw, maxw=maxw, lower=lower, upper=upper, $
        ncontbkpts = ncontbkpts, fibermask=fibermask

	if (NOT keyword_set(minw)) then minw = 3.82
        if (NOT keyword_set(maxw)) then maxw = 3.92
        if (NOT keyword_set(ncontbkpts)) then ncontbkpts = 5

        nfiber = (size(plugsort))[1]
        if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber)
;
;	Use SPECTROPHOTO_STD and REDDEN_STD to correct
;       telluric absorption
;

	ndim = (size(flux))[0]
	if (ndim NE 2) then message, 'expecting 2d flux array'

	npix = (size(flux))[1]
	ntrace = (size(flux))[2]

	tellcorr = fltarr(npix,ntrace) + 1.0

	itell = where((strtrim(plugsort.objtype,2) EQ 'SPECTROPHOTO_STD' OR $
                       strtrim(plugsort.objtype,2) EQ 'REDDEN_STD') AND $
                      (fibermask EQ 0), tellct)

	if (tellct EQ 0) then begin
            print, 'WARNING: No telluric correction stars'
	    return, tellcorr
        endif

;
;	Fill in wavelengths
;
        traceset2xy, wset, pixnorm, wave

	tellwave = wave[*,itell]
	tellflux = flux[*,itell]
	tellfluxivar = fluxivar[*,itell]
;
;	Fit continuum to each telluric standard
;	

        first = 1

        for i=0,tellct-1 do begin	
	   possible = where(tellwave[*,i] GT minw AND tellwave[*,i] LT maxw, $
                          npossible)

	   inside = where(tellwave[*,i] GT minw AND tellwave[*,i] LT maxw $
                          AND tellfluxivar[*,i] GT 0.0, ninside)

	   if (ninside GT 0.7 * npossible) then begin

	     tempwave = tellwave[inside,i]
	     tempflux = tellflux[inside,i]
	     tempivar = tellfluxivar[inside,i]

	     tellfeatures = 1 - (tempwave GT 3.855 AND tempwave LT 3.865 OR $
	                        tempwave GT 3.836 AND tempwave LT 3.842 OR $
	                        tempwave GT 3.875 AND tempwave LT 3.89)

             fullbkpt = slatec_splinefit(tempwave, median(tempflux,21), coeff, $
               maxiter=10, lower=1.0, upper=5.0, mask=mask, $
               invvar=tempivar*tellfeatures, nbkpts=ncontbkpts, rejper = 0.5)
           
	     continuum = slatec_bvalu(tempwave,fullbkpt,coeff)    
	
	     if (first NE 1) then begin 
	       contwave = [contwave, tempwave]
               contflux = [contflux, tempflux / continuum]
               contivar = [contivar, tempivar * continuum^2]
	       medivar = [medivar, fltarr(ninside) + $
                   median(tempivar * continuum^2)]
	     endif else begin
               first = 0
	       contwave = tempwave
               contflux = tempflux / continuum
               contivar = tempivar * continuum^2
	       medivar = fltarr(ninside) + median(contivar)
             endelse
          endif
	endfor

        if (first EQ 1) then begin
          splog,'WARNING: No telluric stars for this region'
          return, tellcorr
        endif
;
;	Scale the variances to give equal weight to each spectrum
;
	maxivar = max(medivar)
	contivar = contivar/(medivar/maxivar)


	ss = sort(contwave)
	contwave = contwave[ss]
	contflux = contflux[ss]
	contivar = contivar[ss]
;
;	Now bspline the features in contflux
;

         telluricbkpt = slatec_splinefit(contwave, contflux, telluriccoeff, $
            maxiter=10, lower=lower, upper=upper, $
            invvar=contivar, everyn=2*tellct)

	  
	 inside = where(wave GT minw AND wave LT maxw)
	 tellcorr[inside] = slatec_bvalu(wave[inside], $
                                telluricbkpt, telluriccoeff)

	return, tellcorr
end
	 


