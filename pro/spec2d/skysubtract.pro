pro skysubtract, obj, objivar, plugmap, wset, skysub, skysubivar, $
            nbkpt=nbkpt, nord=nord, everyn=everyn, allwave=allwave, $
	    allsky=allsky, allfit=allfit

	objsize = size(obj)
	ndim = objsize[0]

	if ndim NE 2 then message, 'obj is not 2d'
	if (size(objivar))[0] NE 2 then message, 'objivar is not 2d'

	ncol = objsize[1]
	nrow = objsize[2]

	if (size(plugmap))[1] NE nrow then $
	   message, 'plugmap does not have same size as nrow'

	if (size(wset.coeff))[2] NE nrow then $
	   message, 'wset does not have same size as nrow'

	traceset2xy,wset,pixnorm,wave
	
	;
	; Find sky fibers
	;
  	
	skies = where(plugmap.objtype EQ 'SKY' AND plugmap.fiberid GT 0, nskies)
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
                     lower=lower, nbkpts=3*ncol/2)
        allfit  = slatec_bvalu(allwave, fullbkpt, coeff)

	skysub = obj - slatec_bvalu(wave, fullbkpt, coeff)


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
	pos67 = 2*binsize/3
	diff67 = fltarr(nn)
	alphav = fltarr(nn)
	for i=0,nn-1 do diff67[i] = (diffr[sort(diffr[*,i]),i])[pos67]
	for i=0,nn-1 do alphav[i] = median(rivar[*,i])

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
	
             
