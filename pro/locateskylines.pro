pro locateskylines, skylinefile, fimage, invvar, invset, $
       xsky, ysky, skywaves, lambda=lambda

	if keyword_set(lambda) then begin 
	  skywaves=10.d^lambda
        endif else begin 
          print,'Reading ',skylinefile
          readcol,skylinefile,skywaves,/silent,form='d'
 	endelse

	nlines = n_elements(skywaves)

	ncoeff = (size(invset.coeff))[1]
;
;	Make xsky, ysky, pairs from skywaves
;
	nrows = (size(fimage))[2]
	xnorm = (2.0d*alog10(skywaves) - (invset.xmin + invset.xmax))/ $
	          (invset.xmax - invset.xmin)

	ysky = dblarr(nrows,nlines)
	xstart = dblarr(nrows,nlines)
	for i=0,nrows-1 do begin
	  ysky[i,*] = float(i)
	  if (invset.func EQ 'legendre') then $
              xstart[i,*] = flegendre(xnorm,ncoeff) # invset.coeff[*,i]
	  if (invset.func EQ 'chebyshev') then $
              xstart[i,*] = fchebyshev(xnorm,ncoeff) # invset.coeff[*,i]
 
	endfor	  


	good = total((xstart LE 0),1) EQ 0 
        gind = where(good)
        if gind[0] eq -1 then message, 'No good fibers!!!'
        xstart=xstart[*,gind]
	ysky=ysky[*,gind]
        skywaves=skywaves[gind]

	ntrace=(size(fimage))[2]
	fmed=fimage-fimage
	for i=0,nTrace-1 do fmed[*,i]=median(fimage[*,i],17)

	xskytmp = trace_fweight(fimage, xstart, ysky, invvar=invvar) 
; iterate trace_fweight ?!?
	xsky =  trace_fweight(fimage, xskytmp, ysky, invvar=invvar) 

;	xsky =  trace_gauss(fimage, xstart, ysky, invvar=invvar) 



	nskyline=(size(xsky))[2]
	mean=fltarr(nskyline)
	sigma=fltarr(nskyline)
	xres=xstart-xsky
	for i=0,nskyline-1 do begin 
	  djs_iterstat,xres[3:150,i],mean=mn,sigma=sig
	  mean[i]=mn
	  sigma[i]=sig
        endfor 
        mean1=mean & sigma1=sigma
	for i=0,nskyline-1 do begin 
	  djs_iterstat,xres[170:316,i],mean=mn,sigma=sig
	  mean[i]=mn
	  sigma[i]=sig
        endfor 
        mean2=mean & sigma2=sigma
	plot,skywaves,mean,xr=[5500,9500],yr=[-.2,.2]
	errplot,skywaves,mean-sigma,mean+sigma
	
	dum=poly_fit(skywaves,mean,2,yfit)
	plot,skywaves,mean-yfit,xr=[5500,9500],yr=[-.3,.3]
	errplot,skywaves,mean-yfit-sigma,mean-yfit+sigma
	print,stdev(mean-yfit)
	
pix=traceset2pix(invset,alog10(skywaves))
pix=pix[160,*]
	dum=poly_fit(skywaves,mean,2,yfit)
	plot,pix,mean-yfit,yr=[-.3,.3]
	errplot,pix,mean-yfit-sigma,mean-yfit+sigma
	print,stdev(mean-yfit)
	pix_space,del


stop
return



stop
trace_gauss, fimage, xstart, ysky, width=width, invvar=invvar

	dif0=xstart[*,0]-xsky[*,0]
	dif1=xstart[*,1]-xsky[*,1]
	dif2=xstart[*,2]-xsky[*,2]
	dif3=xstart[*,3]-xsky[*,3]
	dif4=xstart[*,4]-xsky[*,4]
	dif5=xstart[*,5]-xsky[*,5]
	djs_iterstat, (dif1+dif2)/2, sigma=sigma

kent=gaussfit(x,y,a)

logwav=readfits('redwave.fits')
wav=10.^logwav[*,160]
xs=xstart[160,*]
plot,wav,fimage[*,160]
oplot,skywaves,xs-xs+300,ps=4

; gaussfit
xest=xstart[160,*]
x=findgen(2048)
y=fimage(*,160)
plot,x,y,xr=[1900,2000]
kent=gaussfit(x,y,a,nterms=4,estimate=[800,1904.,1,100])




	return
end 

;gaussfit stuff
;a0 amplitude
;a1 center
;a2 sigma
;a3 offset
;a4 linear terim