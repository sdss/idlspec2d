
pro arcsky_fit, xarc, yarc, arccoeff, xsky, ysky, skycoeff, $
          wskyset, invskyset, function_name=function_name
	if (N_params() LT 6) then begin
          print, 'function arcsky_fit, xarc, yarc, arccoeff, xsky, ysky, '
          print, ' ycoeff, function_name=function_name'
	  print,'function_name can be legendre or chebyshev'
	  return
	endif

	xarc=(xarc)[*]
	xsky=(xsky)[*]
	x=(2.0d*[xarc, xsky]-2047.0d)/2047.0d

	numarcs = n_elements(xarc)
	
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
        y = [yarc, ysky]


	
	beta = transpose(y # legarr)

	alpha = transpose(legarr)#legarr

	svdc, alpha, w, u, v, /double
	
	res = svsol(u, w, v, beta, /double)

        arcfit = flegendre(x, arccoeff) # res[0:arccoeff-1]
        skyres = res[0:arccoeff-1]
        skyres[0:skycoeff-1] = skyres[0:skycoeff-1] +res[arccoeff:*]
        skyfit = flegendre(x, arccoeff) # skyres

wskyset=res

	return
end
	
