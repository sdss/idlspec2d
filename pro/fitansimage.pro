;
;	fitansimage takes the output from extract_image, and smooths
;          the corresponding parameters of nfibers and npoly with
;	   functions of order nord and nordscat, respectively
;
;
function fitansimage, ansimage, nparams, nfibers, npoly, nrows, yrow, $
        fluxm=fluxm, nord=nord, nordscat=nordscat, $
        ymin=ymin, ymax=ymax, fullrows=fullrows

  if (N_params() LT 6) then begin
      print, 'Syntax - fitansimage(ansimage, nparams, nfibers, npoly, '
      print, '  nrows, yrow, fluxm=fluxm, nord=nord, nordscat=nordscat, '
      print, '  ymin=ymin, ymax=ymax, fullrows=fullrows)'
      return, -1
   endif

	anssize = size(ansimage)

;	if(anssize[0] NE 3) then $
;	  message,'Need a 3d ansimage [params,fibers,rows]'

;	nparams = anssize[1]
;	nfibers = anssize[2]
;	nrows = anssize[3]

	if(NOT keyword_set(fluxm)) then $
 			fluxm = make_array(nparams,/long,value=1)
	if(NOT keyword_set(nord)) then nord=2	 ;quadratic
	if(NOT keyword_set(nordscat)) then nordscat=9	 ;scattered light
	if(NOT keyword_set(ymin)) then ymin=0.0	 
	if(NOT keyword_set(ymax)) then ymax=2047.0	 
	if(NOT keyword_set(fullrows)) then fullrows=2048	 

	flux = fltarr(nfibers,nrows)

	ynorm = (2.0*yrow-(ymax+ymin))/(ymax-ymin)
	yfnorm = (2.0*findgen(fullrows)-(ymax+ymin))/(ymax-ymin) 
	fitans = fltarr(nparams*nfibers+npoly,fullrows)
	newflux = fltarr(nfibers,fullrows)
        iTrace = lindgen(nfibers)*nparams
	fitthis = fltarr(nrows,nparams*nfibers + npoly)

	for j=0,nparams-1 do flux = flux + fluxm[j] * ansimage[j+iTrace,*]

	for i=0,nfibers-1 do begin
	    mask = (flux[i,*] GT 0.0)
	    fitans[i*nparams,*] = 1.0
	    newflux[i,*] = 1.0
	  for j=1,nparams-1 do begin
	    good = where(mask)
	    fitthis[good,j+i*nparams] = ansimage[j+i*nparams,good]/flux[i,good]
;
;		Iterate a few times to reject outliers
;
	    done = 0
	    while (done EQ 0) do begin
	      good = where(mask)
	      tt = polyfitw(ynorm[good], fitthis[good,j+i*nparams], $
                    (flux[i,good] > 0), nord, yfit)
              diff = fitthis[good,j+i*nparams] - yfit
	      worst = max(abs(diff),place)
;	      print, i, worst, place, 4*stddev(diff)
	      if (worst LT 4*stddev(diff)) then done = 1 $
              else mask[good[place]] = 0
	    endwhile

	    fitans[j+i*nparams,*] = poly(yfnorm,tt)
	    newflux[i,*] = newflux[i,*] + fluxm[j] * fitans[j+i*nparams,*] 
	  endfor
	endfor

	for j=0,nparams-1 do $
            fitans[j+iTrace,*] = fitans[j+iTrace,*] / newflux 

	;
	;	Now do background terms
	;	First expand terms into nrows x nrows image


	scatimage = fltarr(nrows,nrows)
	scatfit = fltarr(nrows,fullrows)

	for i=0,nrows-1 do $
	  scatimage[*,i] = fchebyshev(ynorm, nPoly, halfintwo=1) $
              # ansimage[nfibers*nparams:*,i]  

	for i=0,nrows-1 do begin
	    tt = poly_fit(ynorm, scatimage[i,*], nordscat, /double)
	    diff = scatimage[i,*] - poly(ynorm,tt)
	    good = where(abs(diff) LT 3*stddev(diff))
	    tt2 = poly_fit(ynorm[good], scatimage[i,good], nordscat, /double)
	    scatfit[i,*] = poly(yfnorm,tt2)
	endfor

;	return, [[[scatfit]],[[scatimage]]]

	for i=0,fullrows-1 do begin
	    fitans[nfibers*nparams:*,i] =  $
               func_fit(ynorm, scatfit[*,i], nPoly, $
                 function_name='fchebyshev', halfintwo=1)
	endfor
	return, fitans

	end   
	  
