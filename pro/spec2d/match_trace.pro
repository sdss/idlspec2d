;+
; NAME:
;   match_trace
;
; PURPOSE:
;   Tweak flat field trace to match object trace
;    with a 2d continuous surface fit to low order
;
; CALLING SEQUENCE:
;   xnow = match_trace(image, invvar, xtrace, [xpoly=, ypoly=, $
;       first=, maxiter= ])
;
; INPUTS:
;   image      - two-dimensional object image
;   invvar     - associated inverse variance image
;   xtrace     - array of best fit flat-field traces [nrow, nfiber]
;
; OPTIONAL KEYWORDS:
;   xpoly      - Order of chebyshev polynomial in first dimension (default 3)
;   ypoly      - Order of chebyshev polynomial in second dimension (default 3)
;   first      - Final fweight centroids of object image
;   maxiter    - Maximum number of rejection iterations for 2d surface fit (default 10)
;
; OUTPUTS:
;   xnow       - Best fit trace positions when tweaked to match image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xnow = match_trace(image, invvar, xtrace)
;   bestlag = median(xnow-xtrace)
;   splog, 'Shifting traces by match_trace ', bestlag
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_reject()
;   fchebyshev()
;   trace_fweight()
;
; REVISION HISTORY:
;   16-Oct-2000  Written by S. Burles, FNAL
;   25-Feb-2002  Added radius keyword, radius=2 is much better when wings of
;                     bright neighboring fibers affects centroiding
;-
;------------------------------------------------------------------------------
function match_trace, image, invvar, xcen, xpoly=xpoly, ypoly=ypoly, $
   first=first, maxiter=maxiter, radius=radius

  if NOT keyword_set(xpoly) then xpoly=3
  if NOT keyword_set(ypoly) then ypoly=3
  if NOT keyword_set(maxiter) then maxiter=10
  if NOT keyword_set(radius) then radius=2

  nparam = xpoly*ypoly

  ny = (size(xcen))[1]
  ntrace = (size(xcen))[2]

  ycen = findgen(ny) # replicate(1,ntrace)

  tmp_xpos = trace_fweight(image, xcen, ycen, invvar=invvar,radius=radius)
  tmp_xpos = trace_fweight(image, tmp_xpos, ycen, invvar=invvar,radius=radius)
  tmp_xpos = trace_fweight(image, tmp_xpos, ycen, invvar=invvar,radius=radius)
  first = trace_fweight(image, tmp_xpos, ycen, xerr=errfirst, $
         invvar=invvar,radius=radius)
  diff = first - xcen 

  invvarfirst = ycen * 0.0
  good = where(errfirst NE 999 AND abs(diff) LT radius, ngood)
  if ngood LT nparam*10 then begin
      splog, 'Can not recenter on new image'
      return, xcen
  endif

  invvarfirst[good] = 1.0/errfirst[good]^2  


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

  sumbad = 0L

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
    sumbad = sumbad + nbad
    splog, "Iteration: ",iiter, " Pixels Rejected: ", nbad, $
           " Finished? ", qdone, format='(a,i3,a,i6,a,i2)'
    if qdone EQ 1 then iiter = maxiter
  endfor

  if sumbad GT 20000L then $
    splog, 'Warning: Large number of pixels rejected! ', sumbad

  return,  xcen + shift
end
 
