;------------------------------------------------------------------------------
pro finalarcfit, x, loglam, wset, ncoeff, ic, nsetcoeff=nsetcoeff, $
 fibermask=fibermask, xweight=xweight, yfit=yfit, $
 maxsig=maxsig, plot=plot, _EXTRA=extra

   dims = size(x, /dimens)
   nfiber = dims[0]
   nlines = dims[1]
   if (n_elements(fibermask) NE nfiber) then $
    fibermask = bytarr(nfiber) 
   if (n_elements(xweight) NE n_elements(x)) then $
    xweight = bytarr(nfiber,nlines) + 1

   lmatrix = loglam # (dblarr(nfiber)+1)
   xy2traceset, transpose(x), lmatrix, wset, ncoeff=ncoeff, $
    invvar=transpose(xweight), _EXTRA=extra, $
    yfit=yfitfirst, xmin=wset.xmin, xmax=wset.xmax

   wsave = wset

   ;------
   ; Fit and replace all the coefficients numbered 1 through NCOEFF-1
   ; with their fit value.  Fit a Chebyshev with NSETCOEFF coefficients, and
   ; a split in the baseline at the central fibers.
  
   fitcoeff = ncoeff - ic 

   coeffmask = (fibermask EQ 0) # (fltarr(fitcoeff) + 1)
   xy2traceset, dindgen(nfiber) # (dblarr(fitcoeff) + 1.0), $
    transpose(wset.coeff[ic:ncoeff-1,*]), tmpset, $
    invvar=coeffmask, func='chebyshev_split', $
    ncoeff=nsetcoeff, maxsig=maxsig, yfit=yfitcoeff
   wset.coeff[ic:ncoeff-1,*] = transpose(yfitcoeff)

   ; Fit the first ic coefficients, keep the others fixed

   ia = bytarr(ncoeff)
   ia[0:ic-1] = 1

   xy2traceset, transpose(x), lmatrix, wset, xmin=wset.xmin, xmax=wset.xmax, $
    ncoeff=ncoeff, invvar=transpose(xweight), $
    yfit=yfit, inputans=wset.coeff, ia=ia, $
    maxsig=2.0, _EXTRA=extra

   fibermed = fltarr(nfiber)
   fibererr = fltarr(nfiber)
   for i=0, nfiber-1 do fibermed[i] = median(lmatrix[*,i] - yfit[*,i])
   for i=0, nfiber-1 do fibererr[i] = $
    sqrt(total((lmatrix[*,i] - yfit[*,i] - fibermed[i])^2)/nfiber)

   if (keyword_set(plot)) then begin
      plot, findgen(nfiber) ## (dblarr(nlines) + 1), lmatrix-yfit, psym=3, $
       yr=[-max(fibererr),max(fibererr)],xr=[-1,nfiber], /xstyle, /ystyle
      oplot, fibermed, psym=1
      errplot, findgen(nfiber),fibermed-fibererr, fibermed+fibererr
  
      plot, 10^lmatrix, (lmatrix-yfit)*(3.0e5/70.0), psym=3, yrange=[-0.1,0.1] 
   endif

   return
end
;------------------------------------------------------------------------------
