;
;	Use fitans to tweak trace and sigma
;
pro tweaktrace, xin, sigmain, centershift, sigmashift, $
                xout, sigmaout

	if (n_params() LT 4) then begin
          print, 'Syntax - tweaktrace, xin, sigmain, centershift, sigmashift,'
          print, '        xout, sigmaout'
	  return
	endif

;
;	For proftypes 1 and 2, the following represents first order shifts
;	
	xout = centershift * sigmain + xin
	sigmaout = (sigmashift + 1.0) * sigmain 
