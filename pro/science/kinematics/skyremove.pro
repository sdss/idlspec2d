Pro skyremove, objflux, objwave, lmin, lmax, ppp, kcor

        pppold=lindgen(n_elements(objflux))
        sky=where(10^objwave GT lmin AND 10^objwave LT lmax)
	outsky=where((10^objwave GT lmin-1000 AND 10^objwave LT lmin) OR $
	             (10^objwave GT lmax AND 10^objwave LT lmax+1000)) 
	fluxsmooth=smooth(objflux,10)
	dif=objflux[outsky]-fluxsmooth[outsky]
	djs_iterstat,dif, sigma=rms, sigrej=5
	maxsky=max(objflux[sky])
        minsky=min(objflux[sky]) 

        if (maxsky GT 5*rms OR abs(minsky) GT 5*rms) then begin

	ppp=where(10^objwave LT lmin OR 10^objwave GT lmax) 
        kcor=1
	endif else begin
        ppp=pppold
        kcor=0
	endelse

return
end 