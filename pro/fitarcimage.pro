function gaussprofile, width, halfsize

        if (NOT keyword_set(width)) then return, 0
        if (NOT keyword_set(halfsize)) then halfsize = long(4 * width)

	fullsize = 2*halfsize + 1
	profile = fltarr(fullsize)
	ii = lindgen(halfsize+1) 
	submm = fltarr(halfsize+1,5)
	for i=0,4 do submm[*,i] = ii + (i-2)*0.2
	denom = 1.06448 * width
	ee = (submm*submm)/(width*width/2.77254) 
	
	subprofile = exp(-ee)/denom
	for i=0,halfsize do begin
           profile[-i + halfsize] = total(subprofile[i,*])/5.0
           profile[i + halfsize] = total(subprofile[i,*])/5.0
	endfor

	return, profile
end

; lampfit,blue10, [[waves],[inten]], guess 
; IDL> lampfit, blue10, [[waves],[inten]], 8.3, [[-0.0003, -0.0002, -0.00001],[1.0e-8,2.0e-8,1.0e-9]]



function lampfit, spec, linelist, guess0, guesshi, width=width, lagwidth = lagwidth, lag=lag, ftol=ftol
	common lampstuff, start, loglamlist, intensity, speccorr
	common lagstuff, llag, bestlag, bestcorr, zeroterm
	common pixstuff, pixarray, profile, middle

	zeroterm = guess0

	if (NOT keyword_set(width)) then width = 2.5
	if (NOT keyword_set(lagwidth)) then lagwidth = 300
	if (NOT keyword_set(ftol)) then ftol = 1.0e-4

	profilesize = 2 * long(width*3) + 1

	profile = gaussprofile(width)
	middle = n_elements(profile)/2

	
	;
	;	Guess is an array of polynomial coefficients
	;

	npix = n_elements(spec)
	pixarray = 2.0*findgen(npix)/(npix-1) - 1.0


	nlines = (size(linelist))[1]
	loglamlist = (alog10(linelist[*,0]))[*]
	intensity = (linelist[*,1])[*]
	pix = fltarr(nlines)

	nsmooth = npix+2*lagwidth
	specsmooth = convol(spec,profile)
	speccorr = fltarr(nsmooth)
	start = lagwidth
	speccorr[start:start+npix-1] = specsmooth
	lag = lindgen(2*lagwidth) - lagwidth
	llag = lag

	p0 = (guesshi[0,*])[*]
	scale = (guesshi[1,*])[*]
	res = amoeba(ftol, p0 = p0, scale=scale, function_name='corrlamps', nmax=200)
	
	ans = [guess0,res]
	nparams = n_elements(ans)
	goffset = fltarr(nparams)

	pixoff = 2.0*bestlag/(npix-1) 
	goffset[0] = ans[1]*pixoff
	if(nparams LT 3) then return, ans-goffset

 	goffset[0] = goffset[0] - ans[2]*pixoff*pixoff	
 	goffset[1] = 2.0*ans[2]*pixoff	

	return, ans-goffset
end	

function corrlamps, a

	common lampstuff, start, loglamlist, intensity, speccorr
	common lagstuff, lag, bestlag, bestcorr, zeroterm
	common pixstuff, pixarray, profile, middle

	guess = [zeroterm, a]
	loglam = poly(pixarray,guess)
	model = speccorr*0.0
	nlines = n_elements(loglamlist)
	for i=0,nlines-1 do begin
           diff = loglam - loglamlist[i]
	   val = min(abs(diff), place)

	   if(place GT 0 AND place LT 2047) then begin	
	     x1 = place+start-middle
	     x2 = place+start+middle
	     model[x1:x2] = model[x1:x2] + profile*intensity[i]
	     endif
	  endfor

	bestcorr = 0.0
	bestlag = 0
	model = sqrt(model)

;	stop
        if (total(model) GT 0) then begin
	  res = c_correlate(model,speccorr, lag)
	  bestcorr = max(res,val)
	  bestlag = lag[val]
	endif
	print, guess, bestcorr, bestlag

	return, 1.0 - bestcorr
end

function fullfit, spec, linelist, guess
	common lagstuff, lag, bestlag, bestcorr, zeroterm

	guess0 = guess[0]
	p0 = guess[1:*]
	scale = abs(p0)*0.3
	scale[0] = scale[0]*0.05

	first = lampfit(spec, linelist, guess0, transpose([[p0],[scale]]), $
	   width = 25.0, lagwidth=250, ftol=1.0e-4)
	
	final = first

	while (abs(bestlag) GT 5)  do begin
	  guess0 = final[0]
	  p0 = final[1:*]
	  scale = abs(p0)*0.1
	  scale[0] = scale[0]*0.2
	  final = lampfit(spec, linelist, guess0, transpose([[p0],[scale]]), $
	     width = 10.0, lagwidth=100, ftol=1.0e-4)
	endwhile

	return,final
end

pro fitarcimage, arc, side, linelist, xnew, ycen, tset, invset, $
                  func=func, ncoeff=ncoeff, ans=ans, lambda=lambda, $
                  thresh=thresh, row=row, goodlines=goodlines

   common lagstuff, lag, bestlag, bestcorr, zeroterm

   if (NOT keyword_set(func)) then func = 'legendre'
   if (NOT keyword_set(ncoeff)) then ncoeff = 5
   if (NOT keyword_set(ans)) then ans = 0

   icoeff = ncoeff + 2

   ndim = size(arc, /n_dim)
   if (ndim NE 2) then $
	message, 'expecting 2-d arc image'
   dims = size(arc, /dim)
   npix = dims[0]
   nTrace = dims[1]

   if (NOT keyword_set(row)) then row = (nTrace-30)/2

   if (NOT keyword_set(thresh)) then begin
     if (side EQ 'blue') then thresh = 200
     if (side EQ 'red') then thresh = 500
   endif
;
;	First trace
;

	xcen = trace_crude(arc, yset=ycen, nave=1, nmed=1, thresh = thresh, $
               maxshifte=1.0)
	xnew = trace_fweight(arc, xcen, ycen)

	bad = where(abs(xnew-xcen) GT 3.0)
	if(bad[0] NE -1) then xnew[bad] = xcen[bad]


	
;
;	Now find wavelength solution
;

	spec = arc[*,row]
;	spec = djs_median(arc[*,row-1:row+1],2)

	if (ans[0] EQ 0) then begin

	  guess = 0
	  if (side EQ 'blue') then guess = [3.68, -0.106, -0.005, 0.005]	
	  if (side EQ 'red') then guess = [3.87, 0.10, -0.003]	
	  if (guess[0] EQ 0) then begin
	    print,'please choose side (red or blue)'
	    return
	  endif

          ans = fullfit(spec, linelist, guess)

	  if (side EQ 'blue' AND bestcorr LT 0.6) then $
            print, 'Initial wavelength solution looks suspicious'
	  if (side EQ 'red' AND bestcorr LT 0.5) then $
            print, 'Initial wavelength solution looks suspicious'

	endif

	wnorm = 2.0*xnew/(npix-1) - 1.0
	wpeak = poly(wnorm, ans)		
	
	thispeak = (wpeak[row,*])[*]
	nlines = n_elements(thispeak)
	lambda = fltarr(nlines)
	val = fltarr(nlines)
	loglamlist = (alog10(linelist[*,0]))[*]
	for i=0,nlines-1 do begin
	  diff = loglamlist - thispeak[i]
	  val[i] = min(abs(diff),place)
;	  print, i, val[i], 10^thispeak[i], 10^loglamlist[place]
	  if(val[i] LT 0.0003 AND linelist[place,2] NE 0.0) then $
              lambda[i] = loglamlist[place]
	endfor

	nonzero = where(lambda GT 0.0, oldcount)
	highones = where(lambda GT 3.9, highct)
	if(oldcount LT 6) then $
	  message, 'only '+string(oldcount)+ ' good arclines found'

	print, 'Found ', oldcount, ' good arc lines'
	if (highct GT 0) then print, '----', highct, ' are above 8000 A'

	xnew = xnew[*,nonzero]
	ycen = ycen[*,nonzero]
        lambda = lambda[nonzero]
	wnorm = 2.0*xnew/(npix-1) - 1.0
	
;
;	Now store best log lambda solutions in tset
;	

     tset = $      
      { func    :    func              , $
        xmin    :    0.0               , $
        xmax    :    0.0               , $
        coeff   :    fltarr(ncoeff, nTrace) $
      }
     invset = $      
      { func    :    func              , $
        xmin    :    0.0               , $
        xmax    :    0.0               , $
        coeff   :    fltarr(icoeff, nTrace) $
      }

	tset.xmin = 0.0
	tset.xmax = 1.0*(npix - 1)

	  if (side EQ 'blue') then begin
	    invset.xmin = 3.57
	    invset.xmax = 3.80
	  endif
	  if (side EQ 'red') then begin
	    invset.xmin = 3.75
	    invset.xmax = 3.97
	  endif
	ymid = 0.5*(invset.xmax + invset.xmin)
	yrange = invset.xmax - invset.xmin
	xx = findgen(2048)
	pixarray = 2.0*findgen(npix)/(npix-1) - 1.0

      if (func EQ 'legendre') then function_name = 'flegendre'
      if (func EQ 'chebyshev') then function_name = 'fchebyshev'

	goodlines = lonarr(oldcount,nTrace) + 1
        

	for i=0,nTrace-1 do begin

	  done = 0

	  while (done EQ 0) do begin
	    
	    use = where(goodlines[*,i] NE 0)
	    res = svdfit((wnorm[i,use])[*], lambda[use], ncoeff, $
               function_name=function_name, singular=singular,yfit=yfit)
	    diff = yfit - lambda[use]

;
;	Take lines within 20 km/s
;
	    bad = where(abs(diff) GT 3.0e-5, nbad)
	    if (nbad EQ 0) then done = 1 $ 
	    else begin
	      maxdiff = max(abs(diff),badplace)
	      goodlines[use[badplace],i] = 0
	    endelse 
	  endwhile

	  tset.coeff[*,i] = res

;
;	Now fit inverse set
;
          if (func EQ 'legendre') then yy = flegendre(pixarray,ncoeff) # res
          if (func EQ 'chebyshev') then yy = fchebyshev(pixarray,ncoeff) # res

          yynorm = 2.0*(yy-ymid)/yrange
          invset.coeff[*,i] = func_fit(yynorm, xx, icoeff, $
              function_name=function_name)

	
          print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])
	endfor	

	return
  end

