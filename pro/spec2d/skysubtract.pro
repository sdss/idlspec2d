
function skysubtract, tt, skyfinal

;
;	small set is [xmin, xmax, coeffs] for wavelengths in log10
;

	npix = (size(tt.data))[1]
	nTrace = (size(tt))[1]

	skies = where (tt.plugmap.objtype EQ 'SKY')
	if (skies[0] EQ -1) then $
	  message, 'No objType SKY in plugmap'

	
	skyfinal = fltarr(npix, nTrace)
	pp = findgen(2048) 

	skystruct = tt[skies]
	nsky = (size(skystruct))[1]

	skyspline = fltarr(npix,nsky)
	for i=0,nsky - 1 do $
	    skyspline[*,i] = spl_init(pp,skystruct[i].data[*,0])


;
;	Make a copy of tt
;	

	ttsub = tt

	for j=0,nTrace -1 do begin
          print, format='($, ".",i4.4,a5)',j,string([8b,8b,8b,8b,8b])

	  smallset = tt[j].coeff
	  smallsize = n_elements(smallset)
          nparams = smallsize - 2
	  xrange = smallset[1] - smallset[0]
          pixarray = (2.0*findgen(npix)-total(smallset[0:1]))/xrange
          waves = flegendre(pixarray, nparams) # smallset[2:*]

;	rebin each sky with proper wavelength solution
	  skyrebin = fltarr(npix,nsky)
	  skymask = make_array(npix,nsky,/long,value=1)


	  for i=0,nsky - 1 do begin
	    wavemin = skystruct[i].icoeff[0]
	    wavemax = skystruct[i].icoeff[1]
	    wavecoeff = skystruct[i].icoeff[2:*]
	    iparams = n_elements(wavecoeff)
	    wavenorm = (2.0*waves - (wavemin+wavemax))/(wavemax-wavemin)
	    pixFromWave = flegendre(wavenorm, iparams) # wavecoeff
	    skyrebin[*,i] = $
              spl_interp(pp,skystruct[i].data[*,0],skyspline[*,i],pixFromWave)

	    offtheedge = where(pixFromWave LT smallset[0] $
                        OR pixFromWave GT smallset[1])
	    if (offtheedge[0] NE -1) then skymask[offtheedge,i] = 0
	  
	  endfor

	  for i=0,npix-1 do begin
	    goodies = where(skymask[i,*])
	    if (goodies[0] NE -1) then  $
	      skyfinal[i,j] = median(skyrebin[i,goodies])
	  endfor

	  ttsub[j].data[*,0] = ttsub[j].data[*,0] - skyfinal[*,j]
	endfor
	
	return, ttsub
	end	
	  
		

