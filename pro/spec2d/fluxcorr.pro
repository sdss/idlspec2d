
function fluxcorr, flux, fluxivar, wset, plugsort, $
        spectrostd=spectrostd, bkptfile=bkptfile, lower=lower, upper=upper

	if (NOT keyword_set(lower)) then lower = 2.5
	if (NOT keyword_set(upper)) then upper = 10
;
;	We'll need a smoothly splined version of the intrinsic
;	f8 subdwarf.
;
	if (keyword_set(spectrostd)) then $
	  filename = findfile(spectrostd) $
	else $
	  filename = djs_locate_file('f8vspline.dat')

	if (filename EQ '') then $
             message, 'cannot fluxcorrect, no intrinsic spectrum'

	readcol, filename[0], f8wave, f8flux

;
;	Now try to find the best spectrophoto_std on this CCD
;	
	
	spectrophoto = where(plugsort.objtype EQ 'SPECTROPHOTO_STD')
	if (spectrophoto[0] EQ -1) then begin
	  print,'FLUXCORR: no spectrophoto stds on this side, trying reddening'
	  spectrophoto = where(plugsort.objtype EQ 'REDDEN_STD')
	  if (spectrophoto[0] EQ -1) then $
	    message, 'FLUXCORR: can not find reddening standards either'
	endif

	nstds = n_elements(spectrophoto)	

;	Now with list of spectrophoto, return best choice by color
;	difference with BD+17 4708:

	bd17mag = [10.56, 9.64, 9.35, 9.25, 9.23]	
	bd17color = bd17mag[0:3]-bd17mag[1:4]
	
	spectrocolor = (plugsort[spectrophoto].mag)[0:3,*] - $
	                   (plugsort[spectrophoto].mag)[1:4,*] 

	bd17big = (fltarr(nstds) + 1.0) ## bd17color	
	colordiff = sqrt(total((spectrocolor - bd17big)^2,1))

	score = min((plugsort[spectrophoto].mag)[2,*]+colordiff*10.0, bestcolor)

	spectroflux = flux[*,spectrophoto[bestcolor]]
	spectrofluxivar = fluxivar[*,spectrophoto[bestcolor]]

	traceset2xy, wset, pixnorm, wave
	spectrowave = wave[*,spectrophoto[bestcolor]]

	mask = spectrowave*0.0 + 1.0

;
;	Mask out absorption line regions
;
       
        absreg = [[3.582, 3.5845],[3.5885, 3.591],[3.594, 3.596], $
                  [3.597, 3.60],[3.605, 3.62],[3.63, 3.645],[3.684, 3.688]]
        nregions = (size(absreg))[2] 
	for i=0,nregions - 1 do begin
          absline = where(spectrowave GT absreg[0,i] AND  $
                           spectrowave LT absreg[1,i])
          if (absline[0] NE -1) then  mask[absline] = 0
        endfor 
;
;	Now find limits of wmin, wmax
;

	wmin = min(wave)
	wmax = max(wave)

;
;	Read in bkpt file, which should have best spacings
;	
	if (keyword_set(bkptfile)) then $
	  filename = findfile(bkptfile) $
	else $
	  filename = djs_locate_file('fluxcorr.bkpts')

	if (filename EQ '') then begin
             print, 'FLUXCORR cannot find bkpt file, filling in here'
	     allbkpts = [ 3.7578, 3.768, 3.772, 3.782, 3.785, $
                   3.79, 3.793, 3.797, 3.805, 3.81, 3.82, 3.83, 3.845, $
                   3.87, 3.90, 3.93, 3.95, 3.97]
	endif else begin
	  readcol, filename, allbkpts
	endelse

	thesebkpts = where(allbkpts GE wmin AND allbkpts LE wmax, numbkpts)	
	if (numbkpts LT 4) then $
	     message, 'FLUXCORR: bkpts are screwed up'

	bkpt = [wmin, allbkpts[thesebkpts[1:numbkpts-2]], wmax]

;
;	Ready to spline spectro std
;	

	fullbkpt = slatec_splinefit(spectrowave, spectroflux, coeff, $
            maxiter=10, lower=lower, upper=upper, $
            invvar=spectrofluxivar*mask, bkpt=bkpt, rejper=0.4)

	intrinspl = spl_init(alog10(f8wave), f8flux)
	fullf8v = spl_interp(alog10(f8wave), f8flux, intrinspl, wave) 

	
;	cheap scaling
;
;	at v=0, flux at 5556A is 956/3.42 photons/cm/A/s
;	        = 9.72e-10 ergs/cm/s/A
;	scale to r', with 10^(-r'/2.5)
;	and return in units to 1e-17 ergs/cm/s/A
;	so factor in exponent is 10^((20.0 - r')/2.5)

	scaling = 10^((20.0 - plugsort[spectrophoto[bestcolor]].mag[2])/2.5)

	fluxfactor = fullf8v / slatec_bvalu(wave, fullbkpt, coeff) * scaling

	return, fluxfactor
end

