;+
; NAME:
;   fluxcorr
;
; PURPOSE:
;   Use sky-subtracted SPECTROPHOTO_STD's to calculate flux correction
;   as a function of wavelength.  We assume the correction is uniform
;   for all fibers (This may change).
;
; CALLING SEQUENCE:
;   ff = fluxcorr(flux, fluxivar, wset, plugsort, color=color, $
;        spectrostd=spectrostd, bkptfile=bkptfile, lower=lower, upper=upper, $
;        fibermask=fibermask)
;
; INPUTS:
;   flux         - sky-subtracted extracted spectra [nx,ntrace]
;   fluxivar     - corresponding inverse variance [nx,ntrace]
;   wset         - wavelength coefficients as a trace set
;   plugsort     - plugmap entries
;
; REQUIRED KEYWORDS:
;   color        - 'red' or 'blue'
;
; OPTIONAL KEYWORDS:
;   spectrostd   - file with spectrophoto continuum, set to 1.0 at 5556 Ang.
;   bkptfile     - file with wavelengths of specific bkpts 
;                     default is (red.bkpts or blue.bkpts)
;   absfile      - file delineating regions of absorption to avoid in
;                     in continuum fitting
;   lower        - lower rejection threshold for continuum fitting
;   upper        - upper rejection threshold for continuum fitting
;   fibermask    - use to reject possible standards which have spectral troubles
;   
; OUTPUTS:
;   ff           - Flux calibration to correct each pixel in flux array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;	  What we really need is an accurate prediction of
;	  f-star continua as a function of u',g',r', and i'
;         Right now, all we can do is estimate.
;         Also need to rewrite to make use of all f-stars.
;       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;
; PROCEDURES CALLED:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/f8vspline.dat
;   $IDLSPEC2D_DIR/etc/f8v.abs
;   $IDLSPEC2D_DIR/etc/blue.bkpts
;   $IDLSPEC2D_DIR/etc/red.bkpts
;
; REVISION HISTORY:
;   19-Oct-1999  Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
function fluxcorr, flux, fluxivar, wset, plugsort, color=color, $
        spectrostd=spectrostd, bkptfile=bkptfile, absfile=absfile, $
        lower=lower, upper=upper, fibermask=fibermask

        nord=3
	if (NOT keyword_set(lower)) then lower = 5
	if (NOT keyword_set(upper)) then upper = 5

        ncol = (size(flux))[1] 
        nrow = (size(flux))[2] 
        if (n_elements(fibermask) NE nrow) then $
           fibermask = bytarr(nrow) 

;
;	We'll need a smoothly splined version of the intrinsic
;	f8 subdwarf.
;
        if (keyword_set(spectrostd)) then $
         filename = spectrostd $
        else $
         filename = filepath('f8vspline.dat', $
          root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')

        filename = (findfile(filename))[0]
        if (filename EQ '') then $
         message, 'Cannot fluxcorrect, no intrinsic spectrum'

        readcol, filename[0], f8wave, f8flux

;
;	Now try to find the best spectrophoto_std on this CCD
;	
	
	spectrophoto = where(strtrim(plugsort.objtype,2) EQ 'SPECTROPHOTO_STD' $
               AND (fibermask EQ 0))
	if (spectrophoto[0] EQ -1) then begin
	  splog, 'WARNING: No spectrophoto stds on this side, trying reddening'
	  spectrophoto = where(strtrim(plugsort.objtype,2) EQ 'REDDEN_STD' $
               AND (fibermask EQ 0))
	  if (spectrophoto[0] EQ -1) then $
	    message, 'FLUXCORR: can not find reddening standards either'
	endif

	nstds = n_elements(spectrophoto)
        splog, 'Number of possible spectrophoto standards = ', nstds

;	Now with list of spectrophoto, return best choice by color
;	difference with BD+17 4708:

	bd17mag = [10.56, 9.64, 9.35, 9.25, 9.23]
	bd17color = bd17mag[0:3]-bd17mag[1:4]

	spectrocolor = (plugsort[spectrophoto].mag)[0:3,*] - $
	                   (plugsort[spectrophoto].mag)[1:4,*] 

	bd17big = (fltarr(nstds) + 1.0) ## bd17color
	colordiff = sqrt(total((spectrocolor - bd17big)^2,1))

	score = min((plugsort[spectrophoto].mag)[2,*]+colordiff*40.0, bestcolor)
	splog, 'Spectrophoto star is fiberid = ', $
         plugsort[spectrophoto[bestcolor]].fiberid
	splog, 'Spectrophoto mag = ', $
         plugsort[spectrophoto[bestcolor]].mag
	splog, 'Spectrophoto score = ', score

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
         filename = absfile $
        else $
         filename = filepath('f8v.abs', $
          root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')

        filename = (findfile(filename))[0]
        if (filename EQ '') then $
         message, 'No F8V Abs file found'
 
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
        if (keyword_set(bkptfile)) then $
         filename = bkptfile $
        else if (keyword_set(color)) then begin 
           if (color EQ 'red') then $
            filename = filepath('red.bkpts', $
             root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
          if (color EQ 'blue') then $
           filename = filepath('blue.bkpts', $
             root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
        endif

        filename = (findfile(filename))[0]
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
;	at v=0, flux at 5556A is 956 photons/cm^2/A/s
;	        = 3.52e-9 ergs/cm^2/s/A
;	scale to r', with 10^(-r'/2.5)
;	and return in units to 1e-17 ergs/cm/s/A
;	so factor in exponent is 10^((21.37 - r')/2.5)
;
;	AB magnitude, the scaling this assumes
;       the AB magnitude of BD+17 is 9.4 at 5560
;
;	we then use f_nu = 10^-0.4*(AB+48.6)
;       and then f_lambda = f_nu * c / (lambda)^2
;  
;       c is 3.0e18 Ang/s
;       lambda is in Ang
;
;
;
;	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;	  What we really need is an accurate prediction of
;	  f-star continua as a function of u',g',r', and i'
;         Right now, all we can do is estimate.
;         Also need to rewrite to make use of all f-stars.
;       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;
;


	scaling = 10^((21.37 - plugsort[spectrophoto[bestcolor]].mag[2])/2.5)

        fullfit = slatec_bvalu(wave, fullbkpt, coeff)

	negs = where(fullfit LE 0.0, nnegs)
	if (nnegs GT 0) then $
	 splog, 'WARNING: Flux factor has negative elements:', nnegs

	fluxfactor = fullfit / (fullf8v * scaling)

	return, fluxfactor
end

