; D. Finkbeiner
; 15 Oct 1999

; Generate QA plots

pro qaplot_arcline, xdif, lambda, arcname

  pmulti = !p.multi
  !p.multi = [0,1,2]
  nline = (size(xdif))[2]
  arcnum = strmid(arcname, 0, strpos(arcname, '.fit'))

  plot,xdif[*,0]*1e3+lambda[0],ps=3,yr=[5000,9000],/yst, $
	xr=[-10,330],/xst, xtit='Fiber Number', ytit='Lambda (A)', $
	title='Arcline Fit for  '+arcnum
  for k=1,nline-1 do oplot,xdif[*,k]*1e3+lambda[k],ps=3


  sig=fltarr(nline)
  for k=0,nline-1 do sig[k]=djsig(xdif[*,k])*1000.

  meandif=fltarr(nline)
  for k=0,nline-1 do begin 
     djs_iterstat, xdif[*,k], mean=mn
     meandif[k]=mn
  endfor

  fiber=140
  xres=xdif[fiber,*]*1e3

  plot, lambda, xres, yr=[-100,100], xtit='lambda', $
	ytit='Delta X (mPix)'
  errplot, lambda, xres-sig, xres+sig

  xyouts, 0.95, 0., systime(), /normal, align=1, chars=0.5

  !p.multi = pmulti
  return
end

