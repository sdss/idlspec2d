
function fitansimage, ansimage, yrow, fluxm=fluxm, npoly=npoly, ymin=ymin, $
	ymax=ymax

	anssize = size(ansimage)

	if(anssize[0] NE 3) then $
	  message,'Need a 3d ansimage [params,fibers,rows]'

	nparams = anssize[1]
	nfibers = anssize[2]
	nrows = anssize[3]

	if(NOT keyword_set(fluxm)) then $
 			fluxm = make_array(nparams,/long,value=1)
	if(NOT keyword_set(npoly)) then npoly=2	 ;quadratic
	if(NOT keyword_set(ymin)) then ymin=0.0	 
	if(NOT keyword_set(ymax)) then ymax=2047.0	 

	flux = fltarr(nfibers,nrows)
	for j=0,nparams-1 do flux = flux + fluxm[j] * ansimage[j,*,*]

	ynorm = (2.0*yrow-(ymax+ymin))/(ymax-ymin)
	fullrows = 2048		;; !!!! Hard wired 2048 !!!!	
	yfnorm = (2.0*findgen(2048)-(ymax+ymin))/(ymax-ymin) 
	fitans = fltarr(nparams,nfibers,fullrows)

	for i=0,nfibers-1 do begin
	    good = where(flux[i,*] NE 0.0)
	    fitans[0,i,*] = 1.0
	  for j=1,nparams-1 do begin
	    fitthis = fltarr(nrows)
	    fitthis(good) = ansimage[j,i,good]/flux[i,good]
	    tt = polyfitw(ynorm, fitthis, flux[i,*], npoly)
	    fitans[j,i,*] = poly(yfnorm,tt)
	  endfor
	endfor

	newflux = fltarr(nfibers,fullrows)
	for j=0,nparams-1 do newflux = newflux + fluxm[j] * fitans[j,*,*]
	for j=0,nparams-1 do fitans[j,*,*] = fitans[j,*,*] / newflux 

	return, fitans
	end   
	  
