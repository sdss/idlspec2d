pro locateskylines, skylinefile, fimage, invvar, invset, $
       xsky, ysky, skywaves

	readcol,skylinefile,skywaves

	nlines = n_elements(skywaves)

	ncoeff = (size(invset.coeff))[1]
;
;	Make xsky, ysky, pairs from skywaves
;
	nrows = (size(fimage))[2]
	xnorm = (2.0*alog10(skywaves) - (invset.xmin + invset.xmax))/ $
	          (invset.xmax - invset.xmin)

	ysky = fltarr(nrows,nlines)
	xstart = fltarr(nrows,nlines)
	for i=0,nrows-1 do begin
	  ysky[i,*] = float(i)
	  if (invset.func EQ 'legendre') then $
              xstart[i,*] = flegendre(xnorm,ncoeff) # invset.coeff[*,i]
	  if (invset.func EQ 'chebyshev') then $
              xstart[i,*] = fchebyshev(xnorm,ncoeff) # invset.coeff[*,i]
 
	endfor	  

	xsky =  trace_fweight(fimage, xstart, ysky, invvar=invvar) 

	return
end 
