;
;	Use fitans to tweak trace and sigma
;
pro tweaktrace, x, sigma, centershift, sigmashift

     if (n_params() LT 4) then begin
       print, 'Syntax - tweaktrace, x, sigma, centershift, sigmashift,'
       return
     endif

;     if (NOT keyword_set(maxshift)) then maxshift=1.0

;
;	For proftypes 1 and 2, the following represents first order shifts
;	
     x = x - centershift * sigma 
     sigma = (sigmashift + 1.0) * sigma 

    return
end
