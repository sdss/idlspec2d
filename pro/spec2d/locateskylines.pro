pro locateskylines, skylinefile, second, secondinvvar, invset, $
       xsky, ysky, skywaves

	readcol,skylinefile,skywaves

	nlines = n_elements(skywaves)

	ncoeff = (size(invset.coeff))[1]
;
;	Make xsky, ysky, pairs from skywaves
;
	nrows = (size(second))[2]
	xnorm = (2.0*alog10(skywaves) - (invset.xmin + invset.xmax))/
	          (invset.xmax - invset.xmin)

	ysky = fltarr(nrows,nlines)
	xsky = fltarr(nrows,nlines)
	for i=0,nrows-1 do begin
	  ysky[i,*] = ysky[i,*] + i
	  if (invset.func EQ 'legendre') then $
              xsky[i,*] = flegendre(xnorm,ncoeff) # invset.coeff[*,i]
	  if (invset.func EQ 'chebyshev') then $
              xsky[i,*] = fchebyshev(xnorm,ncoeff) # invset.coeff[*,i]
 
	endfor	  

	return
end 
