pro skysubtract, obj, objivar, plugmap, wset, skysub, skysubivar, $
            nbkpt=nbkpt, nord=nord, fibermask=fibermask, allwave=allwave, $
	    allsky=allsky, allfit=allfit

	objsize = size(obj)
	ndim = objsize[0]

	if ndim NE 2 then message, 'obj is not 2d'
	if (size(objivar))[0] NE 2 then message, 'objivar is not 2d'

	ncol = objsize[1]
	nrow = objsize[2]

        if (n_elements(fibermask) NE nrow) then $
           fibermask = bytarr(nrow) + 1

	if (size(plugmap))[1] NE nrow then $
	   message, 'plugmap does not have same size as nrow'

	if (size(wset.coeff))[2] NE nrow then $
	   message, 'wset does not have same size as nrow'

	traceset2xy,wset,pixnorm,wave
	
	;
	; Find sky fibers
	;
  	
	skies = where(plugmap.objtype EQ 'SKY' AND plugmap.fiberid GT 0 AND $
                      fibermask, nskies)
	if skies[0] EQ -1 then message, 'no sky fibers in plugmap'

        allwave    =  (wave[*,skies])[*]
        allsky     =  (obj[*,skies])[*]
        allskyivar =  (objivar[*,skies])[*]

;
;	Sort sky points by wavelengths
;
        allwavesort = sort(allwave)
	allwave    = allwave(allwavesort)
	allsky     = allsky(allwavesort)
	allskyivar = allskyivar(allwavesort)
;
;	find nice spline
;
        fullbkpt   = slatec_splinefit(allwave, allsky, coeff, $
                     invvar=allskyivar, maxIter=maxIter, upper=upper, $
                     lower=lower, everyn=nskies)
        allfit  = slatec_bvalu(allwave, fullbkpt, coeff)

	fullfit = slatec_bvalu(wave, fullbkpt, coeff) 
	skysub = obj - fullfit * (objivar GT 0.0)

        ; blue plot
	if (min(wave) LT alog10(5600.0)) then begin
          plot, 10^allwave, allfit, ps=3, xr=[5570,5590], $
           title = 'Sky fibers'
          oplot, 10^allwave, allsky, ps=3
          djs_oplot, 10^allwave, allfit, color='red'
          plot, 10^allwave, allsky-allfit, ps=3, xr=[5570,5590], $
           title = 'Sky subtracted sky fibers ', yr=[-1000,1000]

          plot, 10^allwave, allfit, ps=3, xr=[5570,5590], $
           title = 'All fibers'
          oplot, 10^wave, obj, ps=3 
          djs_oplot, 10^allwave, allfit, color='red

;
;	Let's check flux in 5577
;

	  flux5577 = fltarr(nrow)
	  for i=0,nrow-1 do begin
	    inside = where(wave[*,i] GT alog10(5571) AND $
	                    wave[*,i] LT alog10(5588) AND $
                            objivar[*,i] GT 0.0) 
	    if (inside[0] NE -1 AND fibermask[i]) then begin

	      aa = [16000.0, 3.74655, 1.0e-4, 500.0]
;
;		Might need curvefit for bad pixels
;	      bb = curvefit(wave[inside,i], obj[inside,i], $
;                           objivar[inside,i],aa,function_name='gaussfunct')
                           
	      bb = gaussfit(wave[inside,i], obj[inside,i], aa, nterms=4, $
                           estimates=[16000.0, 3.74655, 1.0e-4, 500.0])
              flux5577[i] = total(bb-aa[3])
            endif

          endfor

	  djs_iterstat,flux5577, median=fluxmed, sigma=fluxsigma
          plot,flux5577,ps=1,yr=[fluxmed-3.0*fluxsigma,fluxmed+3.0*fluxsigma]

	  r2 = plugmap.xFocal^2 + plugmap.yFocal^2
	  plot, r2, flux5577, ps=1,$
            yr=[fluxmed-3.0*fluxsigma,fluxmed+3.0*fluxsigma]

        endif

        ; red plot
	if (max(wave) GT alog10(8000.0)) then begin
          plot, 10^allwave, allsky, ps=3, xr=[8000,8100], $
           title = 'Sky fibers'
          oplot, 10^allwave, allfit
          plot, 10^allwave, allsky-allfit, ps=3, xr=[8000,8100], $
           title = 'Sky subtracted sky fibers '

          plot, 10^wave, obj, ps=3, xr=[8000,8100], $
           title = 'All fibers'
          oplot, 10^allwave, allfit
        endif

;	if (max(wave) GT alog10(8500.0)) then begin

;
;	Now attempt to model variance with residuals on sky
;	This is difficult since variance has noise!
;

	diff = abs(allfit-allsky)*sqrt(allskyivar)

	binsize = nskies
	nn = ncol
	diffr = reform(diff,binsize,nn)
	rivar = reform(allskyivar,binsize,nn)
	rwave = allwave(lindgen(nn)*nskies+nskies/2)

;
;	Try plotting sky residuals as a function of wavelength
;
	plot, 10^allwave, diff, ps=3, yr=[0,10], $
               title='Skysubtraction chi squared'

	pos67 = 2*binsize/3
	diff67 = fltarr(nn)
	alphav = fltarr(nn)
	for i=0,nn-1 do diff67[i] = (diffr[sort(diffr[*,i]),i])[pos67]
	for i=0,nn-1 do alphav[i] = median(rivar[*,i])
        djs_oplot, 10^rwave, diff67, color = 'red'

	alpha = diff67 - median(diff67,2*(nn/40) + 1)
	
	skysubivar = objivar	
	good = where(alphav NE 0.0)
	if good[0] NE -1 then begin
	  deltav = fltarr(nn)
	  deltav[good] = alpha[good]/alphav[good]


;
;	Spline unless we come up with something better
;	
        fullbkpt   = slatec_splinefit(rwave, deltav, coeff, $
                     maxIter=maxIter, upper=30, $
                     lower=30, nbkpts=nn/2)
	within = where(wave GE rwave[0] AND wave LE rwave[nn-1] $
                         AND skysubivar NE 0.0)

        deltaans = slatec_bvalu(wave[within], fullbkpt, coeff)
	skysubivar[within] = 1.0/(1.0/skysubivar[within] + abs(deltaans))

	endif

	return
end
	
             
