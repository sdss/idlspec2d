pro finalarcfit, x, loglam, wset, ncoeff, ic, nsetcoeff=nsetcoeff, $
         maxsig=maxsig, plot=plot, _EXTRA = extra

    nfiber = (size(x))[1]
    nlines = (size(x))[2]

    lmatrix = loglam # (dblarr(nfiber)+1)
    xy2traceset, transpose(x), lmatrix, wset, ncoeff=ncoeff, _EXTRA = extra, $
       yfit=yfit, xmin = wset.xmin, xmax=wset.xmax

    wsave = wset

     ;------
    ; Fit and replace all the coefficients numbered 1 through NCOEFF-1
    ; with their fit value.  Fit a Chebyshev with 8 coefficients, and
    ; a split in the baseline at the central fibers.

    fitcoeff = ncoeff - ic 
    xy2traceset, dindgen(nfiber) # (dblarr(fitcoeff) + 1.0), $
       transpose(wset.coeff[ic:ncoeff-1,*]), tmpset, func='chebyshev', $
       ncoeff=nsetcoeff, maxsig=maxsig, yfit=yfit, /halfintwo
    wset.coeff[ic:ncoeff-1,*] = transpose(yfit)


    ; Fit the first ic coefficients, keep the others fixed

     ia = bytarr(ncoeff)
     ia[0:ic-1] = 1

    xy2traceset, transpose(x), lmatrix, wset, xmin=wset.xmin, xmax=wset.xmax, $
       ncoeff=ncoeff, yfit=yfit, inputans = wset.coeff, ia=ia, $
       maxsig=2.0, _EXTRA = extra

    fibermed = fltarr(nfiber)
    fibererr = fltarr(nfiber)
    for i=0,nfiber - 1 do fibermed[i] = median(lmatrix[*,i] - yfit[*,i])
    for i=0,nfiber - 1 do fibererr[i] = $
          sqrt(total((lmatrix[*,i] - yfit[*,i] - fibermed[i])^2)/nfiber)

    if (keyword_set(plot)) then begin
      plot,findgen(nfiber) ## (dblarr(nlines) + 1),lmatrix-yfit,ps=3, $
       yr=[-max(fibererr),max(fibererr)],xr=[-1,nfiber],/xst,/yst
      oplot,fibermed,ps=1
      errplot,findgen(nfiber),fibermed-fibererr, fibermed+fibererr
  
      plot, 10^lmatrix,(lmatrix-yfit)*(3.0e5/70.0),ps=3,yr=[-0.1,0.1] 
    endif

    return
end




