function makelabel, hdr

	plate = strtrim(string(sxpar(hdr, 'PLATEID')),2)
	camera = strtrim(sxpar(hdr, 'CAMERAS'),2)
	mjd = strtrim(string(sxpar(hdr, 'MJD')),2)
	seqid =  strtrim(string(sxpar(hdr, 'SEQID')),2)
	expos =  strtrim(string(sxpar(hdr, 'EXPOSURE')),2)

	return, plate+'-'+camera+'-'+mjd+'-'+seqid+'-'+expos
end

pro combine2dout, filenames, outputfile, bin, zeropoint, nord=nord, $
        ntrials=ntrials, fullspec=fullspec, fullerr=fullerr, $
        fullwave=fullwave, output=output

	if (NOT keyword_set(bin)) then bin = 0.0001d	
	if (NOT keyword_set(zeropoint)) then zeropoint = 3.5d
	if (NOT keyword_set(nord)) then nord = 4
	if (NOT keyword_set(ntrials)) then ntrials = 25

; filenames = findfile('s-2b-*050*')

	nfiles = (size(filenames,/dimen))[0]

	if (nfiles EQ 0) then return
;
;	read in first one
;
	spec = read2dout(filenames[0], wave, hdr=hdr)
	fullspec = spec[*,0]
	fullerr = spec[*,1]
	fullwave = wave
	bluered = bytarr(n_elements(wave)) + $
             (strpos(sxpar(hdr,'cameras'),'r') EQ 0)

	label =  makelabel(hdr)
	exptime = sxpar(hdr,'EXPTIME')

	for i=1,nfiles - 1 do begin
	  spec = read2dout(filenames[i], wave, hdr=hdr)
	  fullspec = [fullspec, spec[*,0]] 
	  fullerr = [fullerr, spec[*,1]] 
	  fullwave = [fullwave, wave] 
	  bluered = [bluered, bytarr(n_elements(wave)) + $
                (strpos(sxpar(hdr,'cameras'),'r') EQ 0)] 
	  label = [label, makelabel(hdr)]
	  exptime = exptime + sxpar(hdr,'EXPTIME')
	endfor

	totalpix = n_elements(fullspec)
	redpix = where(bluered, numred)
	bluepix = where(bluered EQ 0, numblue)

	if (numblue GT 0 AND numred GT 0) then exptime = exptime * 0.5

	

;
;	Use medians to merge red and blue here
;
	maxblue = max(fullwave(where(bluered EQ 0)))
	minred = min(fullwave(where(bluered EQ 1)))

	scale = 1.0
	if (minred LT maxblue) then begin
          bluecross = where(bluered EQ 0 and fullwave GT minred)
          redcross = where(bluered EQ 1 and fullwave LT maxblue)
	  djs_iterstat, fullspec[bluecross], median=bluemed
	  djs_iterstat, fullspec[redcross], median=redmed
	  scale = bluemed/redmed
	  print, 'COMBINE2DOUT: scaling red by', scale
	  fullspec[redpix] = fullspec[redpix]*scale
	  fullerr[redpix] = fullerr[redpix]*scale
	endif

	spotmin = fix((min(fullwave) - zeropoint)/bin) + 1
	spotmax = fix((max(fullwave) - zeropoint)/bin) 
	wavemin = spotmin * bin + zeropoint
	wavemax = spotmax * bin + zeropoint

	npix = spotmax - spotmin + 1

	newwave = float(dindgen(npix)*bin + wavemin)
	bkpt = float(dindgen(npix/2 + 1)*bin*2.0 + wavemin)
	

;
;	Need to construct ivar for bspline
;
	fullivar = fullerr * 0.0
	nonzero = where(fullerr GT 0.0)
	if (nonzero[0] EQ -1) then $
	   message, 'no good points, all have 0.0 or negative sigma'

	fullivar[nonzero] = 1.0/(fullerr[nonzero]^2)

;
;	Using newwave as breakpoints
;		

	ss = sort(fullwave)
	
	fullbkpt = slatec_splinefit(fullwave[ss], $
              fullspec[ss], coeff, bkpt=bkpt, invvar=fullivar[ss])

	bestguess = slatec_bvalu(newwave, fullbkpt, coeff)

;
;	Below is very dirty Monte Carlo to estimate errors in b-spline
;
	trials = fltarr(ntrials,npix)
	iseed = long(systime(1))
	for i=0,ntrials-1 do begin
	  tempspec = randomu(iseed,totalpix,/normal)*fullerr + fullspec
	  fullbkpt = slatec_splinefit(fullwave[ss], tempspec[ss], $
               coeff, invvar=fullivar[ss], fullbkpt=fullbkpt)
	  trials[i,*] = slatec_bvalu(newwave, fullbkpt, coeff)
	endfor

	besterr = bestguess*0.0
	for i=0,npix-1 do besterr[i] = stddev(trials[*,i])

	output = [[newwave],[bestguess],[besterr]]

;
;	Fix up new header, any one should do to start with
;
	newhdr = hdr

	ncoeff = sxpar(newhdr, 'NWORDER')

	sxaddpar,newhdr, 'NWORDER', 2, 'Linear-log10 coefficients'
	sxaddpar,newhdr, 'WFITTYPE', 'LINEAR-LOG', $
            'Linear-log10 dispersion'
	sxaddpar,newhdr, 'COEFF0', wavemin, $
            'center wavelength (log10) of first pixel'
	sxaddpar,newhdr, 'COEFF1', bin, $
            'log10 dispersion per pixel'
        for i=2,ncoeff-1 do $
          sxdelpar,newhdr, 'COEFF'+strtrim(string(i),2)

        sxdelpar,newhdr, 'PIXMIN'
        sxdelpar,newhdr, 'PIXMAX'

	sxaddpar,newhdr, 'CREATORS', 'Burles & Schlegel (1999) IDLspec', $
                        AFTER='SDSS'

;
;	Now get rid of exposure, and add list of exposures
;

	sxdelpar, newhdr, 'EXPOSURE'
	sxdelpar, newhdr, 'SEQID'
	sxaddpar, newhdr, 'NEXP', nfiles, $
               'Number of exposure in this file', AFTER='TELESCOP'
	for i=0,nfiles-1 do $
	   sxaddpar, newhdr, 'EXPID'+strtrim(string(i),2), label[i], $
               'ID String for exposure '+strtrim(string(i),2), $
                BEFORE='EXPTIME'

	sxaddpar, newhdr, 'EXPTIME', exptime, 'total exposure time (seconds)'
	sxaddpar, newhdr, 'REDSCAL',scale,'Red scaling to match blue overlap', $
                AFTER='EXPTIME'

	writefits, outputfile, [[bestguess],[besterr]], newhdr
	return
end
	  
	  


	
