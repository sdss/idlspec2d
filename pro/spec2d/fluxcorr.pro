
function fluxcorr, flux, fluxivar, wset, plugsort, color=color, $
        spectrostd=spectrostd, bkptfile=bkptfile, lower=lower, upper=upper, $
        fibermask=fibermask

        nord=3
	if (NOT keyword_set(lower)) then lower = 5
	if (NOT keyword_set(upper)) then upper = 5

        ncol = (size(flux))[1] 
        nrow = (size(flux))[2] 
        if (n_elements(fibermask) NE nrow) then $
           fibermask = bytarr(nrow) + 1

;
;	We'll need a smoothly splined version of the intrinsic
;	f8 subdwarf.
;
	if (keyword_set(spectrostd)) then $
	  filename = findfile(spectrostd) $
	else $
	  filename = getenv('IDLSPEC2D_DIR') + '/etc/f8vspline.dat'

	if (filename EQ '') then $
             message, 'cannot fluxcorrect, no intrinsic spectrum'

	readcol, filename[0], f8wave, f8flux

;
;	Now try to find the best spectrophoto_std on this CCD
;	
	
	spectrophoto = where(plugsort.objtype EQ 'SPECTROPHOTO_STD' $
               AND fibermask)
	if (spectrophoto[0] EQ -1) then begin
	  splog, 'WARNING: No spectrophoto stds on this side, trying reddening'
	  spectrophoto = where(plugsort.objtype EQ 'REDDEN_STD' $
               AND fibermask)
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

	score = min((plugsort[spectrophoto].mag)[2,*]+colordiff*40.0, bestcolor)

	spectroflux = flux[*,spectrophoto[bestcolor]]
	spectrofluxivar = fluxivar[*,spectrophoto[bestcolor]]

	traceset2xy, wset, pixnorm, wave
	spectrowave = wave[*,spectrophoto[bestcolor]]

	mask = spectrowave*0.0 + 1.0

;
;	Mask out absorption line regions
;

	filename = ''
	if (keyword_set(absfile)) then $
	  filename = findfile(bkptfile) $
        else $
	  filename = getenv('IDLSPEC2D_DIR') + '/etc/f8v.abs'
      
        if (filename EQ '') then message, 'no F8V Abs file found'
 
        readcol, filename, absmin, absmax

	absreg = transpose([[absmin],[absmax]])
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
;	Read in bkpt file, which should have best spacing
;	
        filename = ''
	if (keyword_set(bkptfile)) then $
	  filename = findfile(bkptfile) $
	else if (keyword_set(color)) then begin 
          if (color EQ 'red') then $
	    filename = getenv('IDLSPEC2D_DIR') + '/etc/red.bkpts'
          if (color EQ 'blue') then $
	    filename = getenv('IDLSPEC2D_DIR') + '/etc/blue.bkpts'
        endif

	if (filename EQ '') then message, 'FLUXCORR cannot find bkpt file'

	readcol, filename, allbkpts

        ; Switch to log10 wavelengths
        allbkpts = alog10(allbkpts)

        ; Add buffer to min and max wavelengths to select
        ; breakpoints interior to extrema
        ; 3.0e-3 in log10 is approx 30 pixels

	thesebkpts = where(allbkpts GE wmin + 3.0e-3 $
                       AND allbkpts LE wmax - 3.0e-3, numbkpts)	

	if (numbkpts LT 4) then $
	     message, 'FLUXCORR: bkpts are screwed up'

	bkpt = [wmin, allbkpts[thesebkpts[1:numbkpts-2]], wmax]

;
;	Ready to spline spectro std
;	


	fullbkpt = slatec_splinefit(spectrowave, spectroflux, coeff, $
            maxiter=10, lower=lower, upper=upper, nord=nord, $
            invvar=spectrofluxivar*mask, bkpt=bkpt, rejper=0.4)

	intrinspl = spl_init(alog10(f8wave), f8flux)
	fullf8v = spl_interp(alog10(f8wave), f8flux, intrinspl, wave) 

        spectrofit = slatec_bvalu(spectrowave, fullbkpt, coeff)

;	cheap scaling
;
;	at v=0, flux at 5556A is 956/3.42 photons/cm/A/s
;	        = 9.72e-10 ergs/cm^2/s/A
;	scale to r', with 10^(-r'/2.5)
;	and return in units to 1e-17 ergs/cm/s/A
;	so factor in exponent is 10^((20.0 - r')/2.5)

	scaling = 10^((20.0 - plugsort[spectrophoto[bestcolor]].mag[2])/2.5)

        fullfit = slatec_bvalu(wave, fullbkpt, coeff)

	negs = where(fullfit LE 0.0, nnegs)
	if (nnegs GT 0) then message, 'Flux factor has negative elements'

	fluxfactor = fullfit / (fullf8v * scaling)

	return, fluxfactor
end

