function fitwithmx, invset, lambda, xpos, nord=nord
;
;	lambda is log 10 wavelength
;	
        if (NOT keyword_set(nord)) then nord=4
	ndim = (size(xpos))[0]
	if (ndim NE 2) then $
            message, 'xpos is not 2d' 
	nfiber = (size(xpos))[1]
	nline = (size(xpos))[2]

; evaluate invset at every lambda
        pix1 = traceset2pix(invset,lambda)

        x=findgen(nfiber)/float(nfiber)
	xnew = xpos

;---------------------------------------------------------------------------
; Poly fits for each arcline
;---------------------------------------------------------------------------

        for i=0,nline-1 do begin
           mx=pix1[*,i]
           dif=xpos[*,i]-mx
           if (!version.release LT '5.4') then begin
              dum=poly_fit(x,dif,nord,yfit)
           endif else begin
              dum=poly_fit(x,dif,nord,yfit=yfit,/double)
           endelse

           res1=yfit-dif
           good=abs(res1) lt 4*stddev(res1,/double)
           good=abs(res1) lt 4*stddev(res1*good,/double)
           if (!version.release LT '5.4') then begin
              kent=polyfitw(x,dif,good,nord,yfit)
           endif else begin
              kent=polyfitw(x,dif,good,nord,yfit,/double)
           endelse

           xnew[*,i] = mx+yfit
        endfor

	return, xnew
end


	

