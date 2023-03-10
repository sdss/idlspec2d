pro plot_extraction_profiles, xcen, sigmacur, fluxvec, fmodel, plottitle, nx, iy, nplotrow=nplotrow
	if not keyword_set(nplotrow) then nplotrow = 9
	device, /color, /portrait, xsize=8.0, ysize=10.0,xoff=0.25, yoff=0.5, /inch
	!p.multi = [0,1,nplotrow]
	xvec = findgen(n_elements(fmodel))
	fmax = max(fmodel)
	yrange = [0,1.1*fmax]
	newline = '!C'
	for j=0, nplotrow-1 do begin
		xrange=nx*[j,j+1]/float(nplotrow)
		xrange[1] = xrange[1]+5
		if (j EQ 0) then title=plottitle + ' at pixel='+strtrim(iy,2)+newline $ 
		else title = ''
		indx = where(xvec GE xrange[0] AND xvec LE xrange[1])
		yrange=[-1.1*abs(min(fmodel[indx])), 1.1*max(fmodel[indx])]
		djs_plot, xvec[indx], fluxvec[indx], $
			xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
			ymargin=[1.5,1.5], xmargin=[1, 1], $
			title=title, charsize=1.2, ytickname=replicate(' ',10), $
			xticklen=0, ticklen=0, yticklen=.005,  XTICK_GET = v
			;ticklen=0, xtickname=replicate(' ',10)
		djs_oplot, xvec[indx], fmodel[indx], color='green'
		if keyword_set(v) then djs_oplot,  v,fltarr(n_elements(v))+.001, psym=1, symsize=.3
		if j eq 0 then begin
			AL_Legend, ['Data', 'Model', 'Sigma', 'Sigma=0'], PSym=[0,0,0,7],charsize=.5, $
				Color=['black','green','red','red'], Position=[5,.95*yrange[1]],SYMSIZE=1
		endif
		findx = where(xcen GE xrange[0] AND xcen LE xrange[1],ct)
		if ct ne 0 then begin
			foreach fi, findx do begin
				x = xcen[fi]
				y = yrange[1]
				yz = .05*yrange[1]
				sig = sigmacur[fi]
				djs_oplot, [x], [y],   color='black', psym=1
				djs_xyouts,[x], [y], strtrim(fi+1,2), color='black', charsize=0.4, ORIENTATION=90
				if sig ne 0 then djs_oplot, [x-sig,x+sig], [yz,yz], color='red' $
				else djs_oplot, [x], [yz], color='red', pysm=7
			endforeach
		endif
	endfor
	!p.multi = [0,1,1]
	device, /color, /portrait;, xsize=7.0, ysize=5.0, xoff=0, yoff=0.5, /inch
end
