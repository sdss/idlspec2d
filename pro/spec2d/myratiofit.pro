
;------------------------------------------------------------------------------
; Compute the polynomial function to multiply by YDATA2 to get YDATA1.
; ??? Need to iterate with the values of COEFF, and applying those to
;     the errors in YDATA2.
function myratiofit, xdata, ydata1, yivar1, ydata2, yivar2, ncoeff=ncoeff

   ndata = n_elements(xdata)
   acoeff = dblarr(ncoeff)
   acoeff[0] = 1.0

   ; Renormalize X to [-1,1]
   xnorm = 2.0 * (xdata - min(xdata)) / (max(xdata) - min(xdata)) -1.0

   ; Construct the matrix A, such that the first row is YDATA1,
   ; the second row is YDATA1 * XNORM, the third is YDATA1 * XNORM^2
   amatrix = fltarr(ndata,ncoeff)
   for icoeff=0, ncoeff-1 do $
    amatrix[*,icoeff] = ydata2 * xnorm^icoeff

   mask1 = yivar1 GT 0 ; =1 for good
   mask2 = yivar2 GT 0 ; =1 for good
   invsig = 1.0 / sqrt( 1.0 / (yivar1 + (1-mask1)) $
                      + 1.0 / (yivar2 + (1-mask2)) ) * mask1 * mask2

   ;----------
   ; Rejection loop

   iiter = 0
   maxiter = 10
   ymodel = 0
   outmask = bytarr(ndata) + 1b ; Begin with all points good
   while (NOT keyword_set(qdone) AND iiter LE maxiter) do begin
      ; At most, reject 5% of all points, distributed in 5 regions
      qdone = djs_reject(ydata1, ymodel, invvar=invsig^2, $
       groupsize=ceil(ndata/5), maxrej=ceil(0.01*ndata), $
       outmask=outmask, upper=5, lower=5)

      mmatrix = amatrix
      for icoeff=0, ncoeff-1 do $
       mmatrix[*,icoeff] = mmatrix[*,icoeff] * invsig * outmask
      bvec = ydata1 * invsig * outmask

      mmatrixt = transpose( mmatrix )
      mm = mmatrixt # mmatrix

      ; Use SVD to invert the matrix
;      mmi = invert(mm, /double)
      if (ncoeff EQ 1) then begin
         mmi = 1.0 / mm
      endif else begin
         svdc, mm, ww, uu, vv, /double
         mmi = 0 * vv
         for i=0L, ncoeff-1 do mmi[i,*] = vv[i,*] / ww[i]
         mmi = mmi ## transpose(uu)
      endelse

      acoeff = mmi # (mmatrixt # bvec)
      yratio = fltarr(ndata) + acoeff[0]
      for icoeff=1, ncoeff-1 do $
       yratio = yratio + acoeff[icoeff] * xnorm^icoeff

      ymodel = yratio * ydata2

      iiter = iiter + 1
   endwhile

   return, yratio
end

;------------------------------------------------------------------------------
