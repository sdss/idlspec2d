
function arcsky_fit, xarc, yarc, arccoeff, xsky, ysky, ycoeff, $
          wskyset, invskyset, function_name=function_name
	if (N_params() LT 6) then begin
          print, 'function arcsky_fit, xarc, yarc, arccoeff, xsky, ysky, '
          print, ' ycoeff, function_name=function_name'
	  print,'function_name can be legendre or chebyshev'
	  return, 0
	endif

	numarcs = n_elements(xarc)
	
	if(NOT keyword_set(function_name)) then $
 	      function_name = 'flegendre'

	if(function_name EQ 'flegendre') then begin
            arclegarr = flegendre([xarc, xsky], arccoeff)
            skylegarr = flegendre([xarc, xsky], skycoeff)
;
;	Zero out these components so the arc lines don't affect the
;	second set of coefficients
;
	    skylegarr[0:numarcs-1,*] = 0.0
	    legarr = [[arclegarr],[skylegarr]]
	endif
	if(function_name EQ 'fchebyshev') then begin
            arclegarr = fchebyshev([xarc, xsky], arccoeff)
            skylegarr = fchebyshev([xarc, xsky], skycoeff)
	    skylegarr[0:numarcs-1,*] = 0.0
	    legarr = [[arclegarr],[skylegarr]]
	endif
	
	beta = transpose(y # legarr)

	alpha = transpose(legarr)#legarr

	svdc, alpha, w, u, v, /double
	
	res = svsol(u, w, v, beta, /double)

	return,res
end
	
