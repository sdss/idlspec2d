function match_trace, image, invvar, xcen, xpoly=xpoly, ypoly=ypoly, $
   first=first, maxiter=maxiter

  if NOT keyword_set(xpoly) then xpoly=3
  if NOT keyword_set(ypoly) then ypoly=3
  if NOT keyword_set(maxiter) then maxiter=10

  nparam = xpoly*ypoly

  ny = (size(xcen))[1]
  ntrace = (size(xcen))[2]

  ycen = findgen(ny) # replicate(1,ntrace)

  tmp_xpos = trace_fweight(image, xcen, ycen, invvar=invvar)
  tmp_xpos = trace_fweight(image, tmp_xpos, ycen, invvar=invvar)
  tmp_xpos = trace_fweight(image, tmp_xpos, ycen, invvar=invvar)
  first = trace_fweight(image, tmp_xpos, ycen, xerr=errfirst, $
         invvar=invvar)

  invvarfirst = ycen * 0.0
  good = where(errfirst NE 999, ngood)
  if ngood LT nparam*10 then begin
      splog, 'Can not recenter on new image'
      return, xcen
  endif

  invvarfirst[good] = 1.0/errfirst[good]^2  


  diff = first - xcen 



;
;	let's do x first
; 

  xmid = 0.5 * (min(xcen) + max(xcen))
  xrange = max(xcen) - min(xcen)
  xnorm = 2.0 * (xcen - xmid) / xrange ; X positions renormalized
  xbasis = fchebyshev(xnorm[*], xpoly)


;
;	let's do y first
; 

  ymid = 0.5 * (min(ycen) + max(ycen))
  yrange = max(ycen) - min(ycen)
  ynorm = 2.0 * (ycen - ymid) / yrange ; Y positions renormalized
  ybasis = fchebyshev(ynorm[*], ypoly)

  full1 = fltarr(nparam,ny*ntrace)
  full2 = fltarr(ny*ntrace, nparam)
  ivar = invvarfirst 
  shift = diff

  for iiter=0, maxiter - 1 do begin
 
 
    for i=0,xpoly - 1 do $
      for j=0, ypoly - 1 do $
        full1[i*ypoly+j,*] = xbasis[*,i] * ybasis[*,j] * ivar
       
    for i=0,xpoly - 1 do $
      for j=0, ypoly - 1 do $
        full2[*,i*ypoly+j] = xbasis[*,i] * ybasis[*,j] 

  
    alpha = full1 # full2
    beta =  full1 # diff[*]

    choldc, alpha, p
    ans = cholsol(alpha,p,beta)
    shift[*] = full2 # ans

    outmask = 0
    qdone = djs_reject(diff, shift, outmask=outmask, $
                 invvar=ivar,upper=8,lower=8)
    ivar = ivar * outmask

    bbad = where(outmask EQ 0, nbad)
    print, iiter, nbad, qdone
    if qdone EQ 1 then iiter = maxiter
  endfor


  return,  xcen + shift
end
 
