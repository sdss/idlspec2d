function telluric_corr,flux, fluxivar, wset, plugsort, $
        minw=minw, maxw=maxw, bkspace=bkpace, lower=lower, upper=upper

	if (NOT keyword_set(minw)) then minw = 3.84
        if (NOT keyword_set(maxw)) then maxw = 3.92
        if (NOT keyword_set(bkspace)) then bkspace = 0.0001
;
;	Use SPECTROPHOTO_STD and REDDEN_STD to correct
;       telluric absorption
;

	ndim = (size(flux))[0]
	if (ndim NE 2) then message, 'expecting 2d flux array'

	npix = (size(flux))[1]
	ntrace = (size(flux))[2]

	tellcorr = fltarr(npix,ntrace) + 1.0

	tell = where(plugsort.objtype EQ 'SPECTROPHOTO_STD' OR $
                     plugsort.objtype EQ 'REDDEN_STD', tellct)

	if (tellct EQ 0) then begin
            print, 'No telluric correction stars'
	    return, tellcorr
        endif

;
;	Fill in wavelengths
;
        traceset2xy, wset, pixnorm, wave

	tellwave = wave[*,tell]
	tellflux = flux[*,tell]
	tellfluxivar = fluxivar[*,tell]
;
;	Fit continuum to each telluric standard
;	


        for i=0,tellct-1 do begin	
	   inside = where(tellwave[*,i] GT minw AND tellwave[*,i] LT maxw)

	   ss = sort(tellwave[inside,i])
	   tempwave = tellwave[inside[ss],i]
	   tempflux = tellflux[inside[ss],i]
	   tempivar = tellfluxivar[inside[ss],i]

	   tellfeatures = 1 - (tempwave GT 3.855 AND tempwave LT 3.865 OR $
	                        tempwave GT 3.875 AND tempwave LT 3.89)

           fullbkpt = slatec_splinefit(tempwave, median(tempflux,21), coeff, $
              maxiter=10, lower=0.8, upper=5.0, mask=mask, $
            invvar=tempivar*tellfeatures, nbkpts=5, rejper = 0.5)
           
	   continuum = slatec_bvalu(tempwave,fullbkpt,coeff)    
	
	   if (i GT 0) then begin 
	     contwave = [contwave, tempwave]
             contflux = [contflux, tempflux / continuum]
             contivar = [contivar, tempivar * continuum^2]
	   endif else begin
	     contwave = tempwave
             contflux = tempflux / continuum
             contivar = tempivar * continuum^2
           endelse

	endfor

	ss = sort(contwave)
	contwave = contwave[ss]
	contflux = contflux[ss]
	contivar = contivar[ss]
;
;	Now bspline the features in contflux
;

         fullbkpt = slatec_splinefit(contwave, contflux, coeff, $
            maxiter=10, lower=lower, upper=upper, $
            invvar=contivar, bkspace=bkspace)

	  
	 inside = where(wave GT 3.853 AND wave LT maxw)
	 tellcorr[inside] = (slatec_bvalu(wave[inside], fullbkpt, coeff) < 1)

	return, tellcorr
end
	 


