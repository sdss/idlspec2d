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
        fullwave=fullwave, output=output, dosky=dosky, wavemin = wavemin, $
        bkptbin = bkptbin, montecarlo=montecarlo

;
;	Set to 50 km/s for now to match 1d
;	Better guess would be 69 km/s
;
	if (NOT keyword_set(bin)) then bin = (50.0/300000.0) / 2.30258

	if (NOT keyword_set(zeropoint)) then zeropoint = 3.5d
	if (NOT keyword_set(nord)) then nord = 4
	if (NOT keyword_set(ntrials)) then ntrials = 25
	if (NOT keyword_set(bkptbin)) then bkptbin = 0.0002

; filenames = findfile('s-2b-*050*')

	nfiles = (size(filenames,/dimen))[0]

	if (nfiles EQ 0) then return
;
;	read in first one
;
	spec = read2dout(filenames[0], wave, hdr=hdr)

	objtype = strtrim(sxpar(hdr, 'OBJTYPE'),2)
	if(objtype EQ 'SKY' AND NOT keyword_set(dosky)) then return

	fullspec = spec[*,0]
	fullerr = spec[*,1]
	fullwave = wave
	nwave = n_elements(wave)
        medwidth =  (wave[2:*] - wave[0:*])/2.0
        width = abs([wave[1] - wave[0],medwidth,wave[nwave-1] - wave[nwave-2]])
	specnum = bytarr(nwave)
	bluered = bytarr(nwave) + (strpos(sxpar(hdr,'cameras'),'r') EQ 0)

	label =  makelabel(hdr)
	exptime = sxpar(hdr,'EXPTIME')

	for i=1,nfiles - 1 do begin
	  spec = read2dout(filenames[i], wave, hdr=hdr)
	  fullspec = [fullspec, spec[*,0]] 
	  fullerr = [fullerr, spec[*,1]] 
	  fullwave = [fullwave, wave] 
	  nwave = n_elements(wave)
          medwidth =  (wave[2:*] - wave[0:*])/2.0
          width = [width, abs([wave[1]-wave[0], medwidth, $
               wave[nwave-1]-wave[nwave-2]])]
	  specnum = [specnum, bytarr(nwave) + i]
	  bluered = [bluered, bytarr(nwave) + $
               (strpos(sxpar(hdr,'cameras'),'r') EQ 0)] 
	
	  label = [label, makelabel(hdr)]
	  exptime = exptime + sxpar(hdr,'EXPTIME')
	endfor

	totalpix = n_elements(fullspec)
	redpix = where(bluered, numred)
	bluepix = where(bluered EQ 0, numblue)

	if (numblue GT 0 AND numred GT 0) then begin
           exptime = exptime * 0.5

;
;	Use medians to merge red and blue here
;
	  maxblue = max(fullwave(where(bluered EQ 0)))
	  minred = min(fullwave(where(bluered EQ 1)))

	  scale = 1.0
	  if (minred LT maxblue) then begin
            bluecross = where(bluered EQ 0 and fullwave GT minred $ 
                      AND fullerr GT 0.0)
            redcross = where(bluered EQ 1 and fullwave LT maxblue $
                      AND fullerr GT 0.0)
	    if (redcross[0] NE -1 AND bluecross[0] NE -1) then begin 
	      djs_iterstat, fullspec[bluecross], median=bluemed
	      djs_iterstat, fullspec[redcross], median=redmed
	      scale = bluemed/redmed

	      if (scale LT 0.1 OR scale GT 10.0) then scale = 1.0
	      print, 'COMBINE2DOUT ', outputfile, ': scaling red by', $
                scale, bluemed/redmed
	      fullspec[redpix] = fullspec[redpix]*scale
	      fullerr[redpix] = fullerr[redpix]*scale
	    endif
	  endif

	endif

	if (NOT keyword_set(wavemin)) then begin
	  spotmin = fix((min(fullwave) - zeropoint)/bin) + 1
	  spotmax = fix((max(fullwave) - zeropoint)/bin) 
	  wavemin = spotmin * bin + zeropoint
	  wavemax = spotmax * bin + zeropoint
	  bkptmin = wavemin
	  bkptmax = wavemax
	endif else begin
	  spotmin = 0
	  spotmax = fix((max(fullwave) - wavemin)/bin) 
	  wavemax = spotmax * bin + wavemin
	  bkptmin = min(fullwave)
	  bkptmax = wavemax
	endelse

	npix = spotmax - spotmin + 1
	nbkpt = fix((bkptmax - bkptmin)/bkptbin) + 1

	newwave = dindgen(npix)*bin + wavemin
	bkpt = dindgen(nbkpt)*bkptbin + bkptmin
	

;
;	Need to construct ivar for bspline
;
	fullivar = fullerr * 0.0
	nonzero = where(fullerr GT 0.0)
	if (nonzero[0] EQ -1) then begin
	   print, 'no good points, all have 0.0 or negative sigma'
	   return
        endif

	fullivar[nonzero] = 1.0/(fullerr[nonzero]^2)

;
;	Using newwave as breakpoints
;		

	ss = sort(fullwave)
	fullbkpt = slatec_splinefit(fullwave[ss], fullspec[ss], coeff, $
              bkpt=float(bkpt), invvar=fullivar[ss], mask=mask, /silent)

	mask[ss] = mask

	bestguess = fltarr(npix)
	inside = where(newwave GE bkptmin AND newwave LE bkptmax, numinside)
	if (inside[0] EQ -1) then $
           message, 'No wavelengths inside breakpoints'

	fwave = float(newwave[inside])
        bestguess[inside] = slatec_bvalu(fwave,fullbkpt,coeff)

      if (NOT keyword_set(montecarlo)) then begin
;
;	Instead of MonteCarlo, let's just guess error with invvar
;	Let's use interpol, hope it's fast enough

	bestivar = bestguess*0.0
	besterr = bestivar

	for i=0,nfiles-1 do begin
	  these = where(specnum EQ i)
	  if (these[0] NE -1) then begin
	    inbetween = where(newwave GE min(fullwave[these]) AND $
	                      newwave LE max(fullwave[these]))
	    if (inbetween[0] NE -1) then begin

	      result = interpol(fullivar[these] * mask[these] / width[these], $
                      fullwave[these], newwave[inbetween])
	
	      bestivar[inbetween] = bestivar[inbetween] + result * bin 
	    endif
	  endif
        endfor

	nonzero = where(bestivar GT 0.0)
	if (nonzero[0] NE -1) then $
          besterr[nonzero] = 1.0/sqrt(bestivar[nonzero])

      endif else begin

;
;	Below is very dirty Monte Carlo to estimate errors in b-spline
;
	trials = fltarr(ntrials,npix)
	iseed = long(systime(1))
	for i=0,ntrials-1 do begin
	  tempspec = randomu(iseed,totalpix,/normal)*fullerr + fullspec
	  fullbkpt = slatec_splinefit(fullwave[ss], tempspec[ss], $
               coeff, invvar=fullivar[ss], fullbkpt=fullbkpt)
	  trials[i,inside] = slatec_bvalu(fwave, fullbkpt, coeff)
	endfor

	besterr = bestguess*0.0
	for i=0,numinside-1 do besterr[inside[i]] = stddev(trials[*,inside[i]])

      endelse

	output = [[newwave],[bestguess],[besterr]]

;
;	Fix up new header, any one should do to start with
;
	newhdr = hdr

	ncoeff = sxpar(newhdr, 'NWORDER')

	sxaddpar,newhdr, 'NWORDER', 2, 'Linear-log10 coefficients'
	sxaddpar,newhdr, 'WFITTYPE', 'LOG-LINEAR', $
            'Linear-log10 dispersion'
	sxaddpar,newhdr, 'COEFF0', wavemin, $
            'center wavelength (log10) of first pixel'
	sxaddpar,newhdr, 'COEFF1', bin, $
            'log10 dispersion per pixel'
        for i=2,ncoeff-1 do $
          sxdelpar,newhdr, 'COEFF'+strtrim(string(i),2)

        sxaddpar,newhdr, 'PIXMIN', 0.000, 'Place holding'
        sxaddpar,newhdr, 'PIXMAX', float(npix - 1), 'Place holding'

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

	sxaddpar, newhdr, 'COMBINE2', systime(), $
                'COMBINE2DOUT finished', AFTER='EXPTIME'

	sxaddpar, newhdr, 'NAXIS1', n_elements(bestguess)
	sxaddpar, newhdr, 'WAT0_001', 'system=linear'
	sxaddpar, newhdr, 'WAT1_001', $
              'wtype=linear label=Wavelength units=Angstroms'
	sxaddpar, newhdr, 'CRVAL1', wavemin, 'Iraf zero point'
	sxaddpar, newhdr, 'CD1_1', bin, 'Iraf dispersion'
	sxaddpar, newhdr, 'CRPIX1', 1, 'Iraf starting pixel'
	sxaddpar, newhdr, 'CTYPE1', 'LINEAR   ' 
	sxaddpar, newhdr, 'WCSDIM', 2
	sxaddpar, newhdr, 'DC-FLAG', 1, 'Log-linear flag'

	writefits, outputfile, [[bestguess],[besterr]], newhdr
	return
end
	  
	  


	
