; Inputs: xnew, ycen, lambda
; Outputs: goodlines, wset, invset
; optional input keywords ymin, ymax, func, linesearch

; linesearch  - set for quick linesearch.  Once a line has been rejected in 
;                  one fiber, never use it again. 


pro fit_tset, xnew, ycen, lambda, goodlines, wset, invset, $
	ymin=ymin, ymax=ymax, func=func, linesearch=linesearch


   if keyword_set(func) eq 0 then func = 'legendre'
   if keyword_set(ymin) eq 0 then ymin = wset.xmin
   if keyword_set(ymax) eq 0 then ymax = wset.xmax

    ntrace=(size(xnew))[1]
    ncoeff=(size(wset.coeff))[1]
    icoeff=(size(invset.coeff))[1]
    npix=2048	

    ymid = 0.5*(invset.xmax + invset.xmin)
    yrange = invset.xmax - invset.xmin
    xx = dindgen(2048)
    pixarray = 2.0d0*dindgen(npix)/(npix-1) - 1.0d0

    if (func EQ 'legendre')  then f_pixarray = flegendre(pixarray,ncoeff)
    if (func EQ 'chebyshev') then f_pixarray = fchebyshev(pixarray,ncoeff)

    if (func EQ 'legendre')  then function_name = 'flegendre'
    if (func EQ 'chebyshev') then function_name = 'fchebyshev'

;wavnorm = 2.0d*xnew/(npix-1) - 1.0d
	wavnorm = 2.0d*xnew/(ymax-ymin) - 1.0d

        nline = n_elements(lambda)

	goodlines = lonarr(nline,nTrace) + 1
	mastergood = lonarr(nline) + 1

	for i=0,nTrace-1 do begin

	  done = 0

	  while (done EQ 0) do begin
	    
            if keyword_set(linesearch) then begin 
               use = where((goodlines[*,i] NE 0) and (mastergood NE 0))
            endif else begin 
	       use = where(goodlines[*,i] NE 0)
            endelse 

	    res = svdfit((wavnorm[i,use])[*], lambda[use], ncoeff, $
               function_name=function_name, singular=singular,yfit=yfit)
	    diff = yfit - lambda[use]

;
;	Take lines within 20 km/s - throw out 1 at a time
;
	    bad = where(abs(diff) GT 3.0d-5, nbad)
	    if (nbad EQ 0) then done = 1 $ 
	    else begin
	      maxdiff = max(abs(diff),badplace)
	      goodlines[use[badplace],i] = 0
              mastergood[use[badplace]] = 0
	    endelse 
	  endwhile

	  wset.coeff[*,i] = res

;
;	Now fit inverse set
;

	yy = f_pixarray # res
;
;	Fit to the 2048 pixels with the wavelengths as the dependent variable.
;
          yynorm = 2.0*(yy-ymid)/yrange
          invset.coeff[*,i] = func_fit(yynorm, xx, icoeff, $
              function_name=function_name)

          print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])
	endfor	


  return
end





