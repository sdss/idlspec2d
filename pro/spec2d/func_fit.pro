
function func_fit, x, y, ncoeff, function_name=function_name

	if (N_params() LT 3) then begin
	  print,'function func_fit, x, y, ncoeff, function_name=function_name' 
	  print,'function_name can be legendre or chebyshev'
	  return, 0
	endif

	
	if(NOT keyword_set(function_name)) then $
 	      function_name = 'flegendre'

	if(function_name EQ 'flegendre') then legarr = flegendre(x, ncoeff)
	if(function_name EQ 'fchebyshev') then legarr = fchebyshev(x, ncoeff)

	
	beta = transpose(y # legarr)

	alpha = transpose(legarr)#legarr

	svdc, alpha, w, u, v, /double
	
	res = svsol(u, w, v, beta, /double)

	return,res
end
	
