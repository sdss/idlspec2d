function shift_trace, image, xcen, ycen, lagrange=lagrange, lagstep=lagstep, $
        maxiter=maxiter

        if (NOT keyword_set(lagrange)) then lagrange = 1.0
        if (NOT keyword_set(lagstep)) then lagstep = 0.1
        if (NOT keyword_set(maxiter)) then maxiter = 3

	nlag = 2L*long(lagrange/lagstep) + 1
        lagx = findgen(nlag)*lagstep - lagrange

        nfiber = (size(xcen))[2]
        
        totalflux = fltarr(nfiber, nlag)

        xstart = xcen
        xshift = 0.0

        for iter=0,maxiter - 1 do begin
	  for i=0, nlag - 1 do begin
            xcur = xstart + lagx[i]

	    fex = extract_boxcar(image, xcur, ycen, radius=1.0)

	    totalflux[*,i] = djs_median(fex,1)
          endfor

          fulltotal = total(totalflux,1)
	  finenlag = 20L*long(lagrange/lagstep) + 1
          finelagx = findgen(finenlag)*lagstep*0.1 - lagrange

          ; --------------------------------------------
          ;  This is just for fun
          ;
          maxt = fltarr(nfiber)
          for i=0,nfiber -1 do begin
            a = poly_fit(lagx, totalflux[i,*], 5)
            peak = max(poly(finelagx, a), place)
            maxt[i] = finelagx[place]
          endfor

          a = poly_fit(lagx, fulltotal, 5)
          totalfit = poly(finelagx, a)
	  peak = max(totalfit, place)

          if (place EQ 0) then begin
            splog, 'max is at start of lagrange...shifting'
            xshift = xshift - lagrange
            xstart = xcen - lagrange
            
          endif else if (place EQ finenlag - 1) then begin
            splog, 'max is at end of lagrange...shifting'
            xshift = xshift + lagrange
            xstart = xstart + lagrange

;
;	things look good return best shift
;
	  endif else return, finelagx(place) + xshift
       endfor

       splog, 'No good shift found after '+string(maxiter)+ ' iterations'
       return, 0.0
end

	
