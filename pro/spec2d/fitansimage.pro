;
;	fitansimage takes the output from extract_image, and smooths
;          the corresponding parameters of nfibers and npoly with
;	   functions of order nord and nordscat, respectively
;
;
function fitansimage, ansimage, nparams, nfibers, npoly, nrows, yrow, $
        fluxm=fluxm, nord=nord, nordscat=nordscat, $
        ymin=ymin, ymax=ymax, fullrows=fullrows, crossfit=crossfit,  $
        scatimage = scatimage, scatfit=scatfit, mingood = mingood

  if (N_params() LT 6) then begin
      print, 'Syntax - fitansimage(ansimage, nparams, nfibers, npoly, '
      print, '  nrows, yrow, fluxm=fluxm, nord=nord, nordscat=nordscat, '
      print, '  ymin=ymin, ymax=ymax, fullrows=fullrows, crossfit=crossfit)'
      return, -1
   endif

	anssize = size(ansimage)

	if(NOT keyword_set(fluxm)) then $
 			fluxm = make_array(nparams,/long,value=1)
	if(NOT keyword_set(nord)) then nord=2	 ;quadratic
	if(NOT keyword_set(nordscat)) then nordscat=9	 ;scattered light
	if(NOT keyword_set(ymin)) then ymin=0.0	 
	if(NOT keyword_set(ymax)) then ymax=2047.0	 
	if(NOT keyword_set(fullrows)) then fullrows=2048	 
	if(keyword_set(crossfit)) then niter = 2 $
	else niter = 1

	flux = fltarr(nfibers,nrows)

        if (NOT keyword_set(mingood)) then mingood = nrows/3

	ynorm = (2.0*yrow-(ymax+ymin))/(ymax-ymin)
	yfnorm = (2.0*findgen(fullrows)-(ymax+ymin))/(ymax-ymin) 
	smallans = fltarr(nparams*nfibers+npoly,nrows)
	fitans = fltarr(nparams*nfibers+npoly,fullrows)
	smallflux = fltarr(nfibers,nrows)
	newflux = fltarr(nfibers,fullrows)
        iTrace = lindgen(nfibers)*nparams
        iParams = lindgen(nparams) 
	iFiber = findgen(nfibers)
	fitthis = fltarr(nrows,nparams*nfibers + npoly)

	for j=0,nparams-1 do flux = flux + fluxm[j] * ansimage[j+iTrace,*]

	for i=0,nfibers-1 do begin
	    mask = (flux[i,*] GT 0.0)
	    fitans[i*nparams,*] = 1.0
	    smallans[i*nparams,*] = 1.0
	  for j=1,nparams-1 do begin
	    good = where(mask)
            if (good[0] NE -1) then $
	      fitthis[good,j+i*nparams] = ansimage[j+i*nparams,good]/flux[i,good]
          endfor
       endfor

       for iter=1,niter do begin 
         for i=0,nfibers-1 do begin
	   for j=1,nparams-1 do begin
;
;		Iterate a few times to reject outliers
;
	      done = 0
	      while (done EQ 0) do begin
	        good = where(mask,ngood)
                if (ngood LT mingood) then done = 1 $
                else begin
	          tt = polyfitw(ynorm[good], fitthis[good,j+i*nparams], $
                    (flux[i,good] > 0), nord, yfit)
                  diff = fitthis[good,j+i*nparams] - yfit
	          worst = max(abs(diff),place)
;	          print, i, worst, place, 4*stddev(diff)
	          if (worst LE 4*stddev(diff)) then done = 1 $
                  else mask[good[place]] = 0
                endelse
	      endwhile

              if (ngood GE mingood) then begin
	        smallans[j+i*nparams,*] = poly(ynorm,tt)
	        fitans[j+i*nparams,*] = poly(yfnorm,tt)
              endif
	    endfor
	  endfor


      ;----------------------------------------------------------
      ;  Iterate one more time (crossfit) with a spline vs. fiber number to
      ;  obtain a much smoother parameter map

        if (keyword_set(crossfit) AND iter EQ 1) then begin

         for i=0,nfibers - 1 do $
	    smallflux[i,*] = fluxm # smallans[iParams+i*nparams,*]

	 for i=0,nrows-1 do begin
           print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])

           for j=1,nparams-1 do begin
;
;		Iterate a few times to reject outliers
;
	      mask = (smallflux[*,i] GT 0.0)
	      done = 0
	      while (done EQ 0) do begin
	        good = where(mask,ngood)
                if (ngood LT nfibers/3) then begin
                      done = 1 
                      splog, 'only a unmasked fibers in row ',i
                endif else begin
	          these = iTrace[good] + j
	          tt = polyfitw(iFiber[good], smallans[these,i], $
                    (smallflux[good,i] > 0), nord, yfit)
                  diff = (smallans[these,i] - yfit)*(smallflux[good,i] > 0)
	          worst = max(abs(diff),place)
;	          print, i, worst, place, 4*stddev(diff)
	          if (worst LE 4*stddev(diff)) then done = 1 $
                  else mask[good[place]] = 0
                endelse
	      endwhile

              if (ngood GE nfibers/3) then $
	         fitthis[i,iTrace+j] = poly(iFiber,tt)

;             fullbkpt = slatec_splinefit(iFiber, smallans[iTrace+j,i], $
;                coeff, nbkpt = 5, invvar=(smallflux[good,i] > 0))
;	     fitthis[i,iTrace+j] = slatec_bvalu(iFiber, fullbkpt, coeff)
	   endfor
         endfor
       endif


      endfor  ;niter loop


      ;--------------------------------------------------------
      ;  Normalize fitans with total flux
      ; 

      for i=0,nfibers - 1 do $
	    newflux[i,*] = fluxm # fitans[iParams+i*nparams,*]

      for j=0,nparams-1 do $
            fitans[j+iTrace,*] = fitans[j+iTrace,*] / newflux 


	;---------------------------------------------------
	;	Now do background terms
	;	First expand terms into nrows x nrows image


	scatimage = fltarr(fullrows,nrows)
	scatfit = fltarr(fullrows,fullrows)

	for i=0,nrows-1 do $
	  scatimage[*,i] = fchebyshev(yfnorm, nPoly, halfintwo=1) $
              # ansimage[nfibers*nparams:nfibers*nparams+npoly-1,i]  

	for i=0,fullrows-1 do begin
            fullbkpt = slatec_splinefit(ynorm, scatimage[i,*], coeff, $
                     nbkpts = 20)
;	    tt = poly_fit(ynorm, scatimage[i,*], nordscat, /double)
;	    diff = scatimage[i,*] - poly(ynorm,tt)
;	    good = where(abs(diff) LT 3*stddev(diff))
;            print, n_elements(good)
;	    tt2 = poly_fit(ynorm[good], scatimage[i,good], nordscat, /double)
;	    scatfit[i,*] = poly(yfnorm,tt2)
	    scatfit[i,*] = slatec_bvalu(yfnorm,fullbkpt,coeff)
	endfor

;	return, [[[scatfit]],[[scatimage]]]

	for i=0,fullrows-1 do begin
	    fitans[nfibers*nparams:*,i] =  $
               func_fit(yfnorm, scatfit[*,i], nPoly, $
                 function_name='fchebyshev', halfintwo=1)
	endfor
	return, fitans

	end   
	  
