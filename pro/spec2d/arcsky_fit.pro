
function arcsky_fit, x, y, numarcs, arccoeff, skycoeff, $
          function_name=function_name, yfit=yfit
	if (N_params() LT 5) then begin
          print, 'function arcsky_fit(xarc, yarc, arccoeff, xsky, ysky, '
          print, ' ycoeff, function_name=function_name)'
	  print,'function_name can be flegendre or fchebyshev'
	  return, -1
	endif

	
	if(NOT keyword_set(function_name)) then $
 	      function_name = 'flegendre'

	if(function_name EQ 'flegendre') then begin
            arclegarr = flegendre(x, arccoeff)
            skylegarr = flegendre(x, skycoeff)
;
;	Zero out these components so the arc lines don't affect the
;	second set of coefficients
;
	    skylegarr[0:numarcs-1,*] = 0.0
	    legarr = [[arclegarr],[skylegarr]]
	endif
	if(function_name EQ 'fchebyshev') then begin
            arclegarr = fchebyshev(x, arccoeff)
            skylegarr = fchebyshev(x, skycoeff)
	    skylegarr[0:numarcs-1,*] = 0.0
	    legarr = [[arclegarr],[skylegarr]]
	endif
	
	beta = transpose(y # legarr)

	alpha = transpose(legarr)#legarr

	svdc, alpha, w, u, v, /double
	
	res = svsol(u, w, v, beta, /double)

	yfit = legarr # res
	return, res
end
	
