function fitwithmx, invset, lambda, xpos, nord=nord
;
;	lambda is log 10 wavelength
;	
        if (NOT keyword_set(nord)) then nord=4
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
           dum=poly_fit(x,dif,nord,yfit)
           res1=yfit-dif
           good=abs(res1) lt 4*stddev(res1)
           good=abs(res1) lt 4*stddev(res1*good)
           kent=polyfitw(x,dif,good,nord,yfit)   
           xnew[*,i] = mx+yfit
        endfor

	return, xnew
end


	

