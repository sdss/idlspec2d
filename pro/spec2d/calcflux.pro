pro calcflux, ansrow, prow, fluxrow, finvrow, wfixed, proftype, lTrace,nCoeff, $
            sigmacur, xcencur, corcalc 

     if(proftype EQ 3 OR proftype EQ 4) then begin
       fluxrow = ansrow[0,*]
       if(nCoeff GE 2) then fluxrow = fluxrow + ansrow[1,*]
       if(nCoeff GE 3) then fluxrow = fluxrow + ansrow[2,*]
       
; best estimate we can do
       finvrow = prow[lTrace*nCoeff]*prow[lTrace*nCoeff] 

;
;	Higher order symmetric terms are not yet included
;
       return
     endif

     if(proftype EQ 1 OR proftype EQ 2) then begin
;       print, 'Analyzing row', cur, '     With Gaussian Profile', wfixed

       fluxrow = ansrow[0,*]
       finvrow = prow[lTrace*nCoeff]*prow[lTrace*nCoeff] ; best estimate we can do
					      ; without covariance matrix

       if(nCoeff GE 2) then begin 	      ; add in symmetric term if present
	  widthrow = ansrow[1,*]
          widthinvvar = prow[lTrace*nCoeff + 1]*prow[lTrace*nCoeff + 1]
          fluxrow = fluxrow + widthrow

;
;	Estimate new widths if specified
;
          if(wfixed[1] GT 0 AND keyword_set(sigmacor)) then begin 
;             print, 'Calculating SIGMACOR...'

;
;	 Make a guess at an underestimated error
;
	     safe = where(fluxrow GT 0.0, safecount)
             if (safecount GT 0) then begin
                r = widthrow[safe]/fluxrow[safe]
                rinvvar = fluxrow[safe] * fluxrow[safe] * $
                        finvrow[safe]  * widthinvvar[safe]
	        nz = where(rinvvar NE 0.0, nonzerocount)
	 
                if (nonzerocount GT 0) then $
                  rinvvar[nz] = rinvvar[nz] / (r[nz]*r[nz] * $
                     finvrow[safe[nz]] + widthinvvar[safe[nz]])

;
;		Only take corrections significant at 2 sigma
;
	        check = where(rinvvar GT 4.00 AND abs(r) LT 0.4, count)
                if(count GT 0) then $
	          sigmacur[safe[check]] = sigmacur[safe[check]] * $ 
                   (r[check]+ 1.0)
             endif
             sigmacor[iy,*] = sigmacur
          endif

       endif

;
;	Estimate new centroids if specified
;
       if(nCoeff GE 3) then  $ ; calculate asymmetric term if present
       if(wfixed[2] GT 0 AND keyword_set(xcencor)) then begin 
;          print, 'Calculating XCENCOR...'

	 centerrow = ansrow[2,*]
         centinvvar = prow[lTrace*nCoeff + 2]*prow[lTrace*nCoeff + 2]

;
;	 Make a guess at an underestimated error
;
	 safe = where(fluxrow GT 0.0, safecount)
         if (safecount GT 0) then begin
            r = centerrow(safe)/fluxrow(safe)
            rinvvar = fluxrow[safe] * fluxrow[safe] * $
                        finvrow[safe]  * centinvvar[safe]
            nz = where(rinvvar NE 0.0, nonzerocount)

            if (nonzerocount GT 0) then $
               rinvvar[nz] = rinvvar[nz] / (r[nz]*r[nz] * $
                     finvrow[safe[nz]] + centinvvar[safe[nz]])
;
;		Only take corrections significant at 2 sigma
;
	    check = where(rinvvar GT 4.00 AND abs(r) LT 0.4, count)
            xcencurrent = fltarr(nTrace)
            if(count GT 0) then $
	       xcencurrent(safe(check)) = r(check) * sigmacur(safe(check))
          endif
          xcencor[iy,*] = xcencurrent
       endif
      return
   endif

   return
end
;------------------------------------------------------------------------------
