function makelabel, hdr

	plate = strtrim(string(sxpar(hdr, 'PLATEID')),2)
	camera = strtrim(sxpar(hdr, 'CAMERAS'),2)
	mjd = strtrim(string(sxpar(hdr, 'MJD')),2)
	seqid =  strtrim(string(sxpar(hdr, 'SEQID')),2)
	expos =  strtrim(string(sxpar(hdr, 'EXPOSURE')),2)

	return, plate+'-'+camera+'-'+mjd+'-'+seqid+'-'+expos
end

pro combine2dout, filenames, outputroot, bin, zeropoint, nord=nord, $
        ntrials=ntrials, fullspec=fullspec, fullerr=fullerr, $
        fullwave=fullwave, output=output, dosky=dosky, wavemin = wavemin, $
        bkptbin = bkptbin, everyn=everyn, display=display, window=window

;
;	Set to 50 km/s for now to match 1d
;	Better guess would be 69 km/s
;
	if (NOT keyword_set(bin)) then bin = (69.0/299792.5) / 2.30258

	if (NOT keyword_set(zeropoint)) then zeropoint = 3.5d
	if (NOT keyword_set(nord)) then nord = 3
	if (NOT keyword_set(ntrials)) then ntrials = 25
	if (NOT keyword_set(bkptbin)) then bkptbin = bin

; filenames = findfile('s-2b-*050*')

	nfiles = (size(filenames,/dimen))[0]

	if (nfiles EQ 0) then return
;
;	read in first one
;

	flux     = mrdfits(filenames[0], 0, hdr)
	fluxivar = mrdfits(filenames[0], 1)
	plugmap  = mrdfits(filenames[0], 2)
	wset     = mrdfits(filenames[0], 3)
        traceset2xy, wset, pixnorm, wave

        npix  = (size(flux))[1]
        nfiber = (size(flux))[2]
        specnum = bytarr(npix)
	bluered = bytarr(npix) + (strpos(sxpar(hdr,'cameras'),'r') EQ 0)
	label =  makelabel(hdr)
	exptime = sxpar(hdr,'EXPTIME')

	for i=1,nfiles - 1 do begin
          tempflux = mrdfits(filenames[i], 0, hdr)
	  tempivar = mrdfits(filenames[i], 1)
	  tempplug = mrdfits(filenames[i], 2)
	  tempwset = mrdfits(filenames[i], 3)
          traceset2xy, tempwset, pixnorm, tempwave

          npix     = (size(tempflux))[1]
	  flux     = [flux, tempflux]

          if (keyword_set(window)) then $
               tempivar[0:window] = tempivar[0:window]*findgen(window+1)/window

	  fluxivar = [fluxivar, tempivar]
          wave     = [wave,tempwave]
          specnum = [specnum, bytarr(npix) + i]
	  bluered  = [bluered, bytarr(npix) + $ 
                           (strpos(sxpar(hdr,'cameras'),'r') EQ 0)] 
	  
	  label = [label, makelabel(hdr)]
	  exptime = exptime + sxpar(hdr,'EXPTIME')
        endfor

	totalpix = (size(flux))[1]
        nfiber = (size(flux))[2]

	redpix = where(bluered, numred)
	bluepix = where(bluered EQ 0, numblue)
	if (numblue GT 0 AND numred GT 0) then $
               exptime = exptime * 0.5
;
;	Fix up new header, any one should do to start with
;

	ncoeff = sxpar(hdr, 'NWORDER')
        for i=2,ncoeff-1 do sxdelpar,hdr, 'COEFF'+strtrim(string(i),2)


        sxaddpar,hdr, 'PIXMIN', 0.000, 'Place holding'
        sxaddpar,hdr, 'PIXMAX', float(npix - 1), 'Place holding'

	sxaddpar,hdr, 'CREATORS', 'Burles & Schlegel (1999) IDLspec', $
                        AFTER='SDSS'

;
;	Now get rid of exposure, and add list of exposures
;

	sxdelpar, hdr, 'EXPOSURE'
	sxdelpar, hdr, 'SEQID'
	sxaddpar, hdr, 'NEXP', nfiles, $
               'Number of exposure in this file', AFTER='TELESCOP'
	for i=0,nfiles-1 do $
	   sxaddpar, hdr, 'EXPID'+strtrim(string(i),2), label[i], $
               'ID String for exposure '+strtrim(string(i),2), $
                BEFORE='EXPTIME'

	sxaddpar, hdr, 'EXPTIME', exptime, 'total exposure time (seconds)'
	sxaddpar, hdr, 'COMBINE2', systime(), $
                'COMBINE2DOUT finished', AFTER='EXPTIME'

        scale = fltarr(nfiber)
        blueflux = fltarr(nfiber)
        redflux = fltarr(nfiber)

        for i=0,nfiber - 1 do begin
 
	  scale[i] = 1.0
;	  if (strtrim(plugmap[i].objtype,2) EQ 'SKY' AND $
;             NOT keyword_set(dosky)) then $
;                splog, ' skipping sky on fiber ', i+1 $
;          else if (strtrim(plugmap[i].objtype,2) EQ 'NA') then $
;               splog, ' skipping bad fiber ', i+1 $
;          else begin

          splog, plugmap[i].fiberid, ' ', plugmap[i].objtype, plugmap[i].mag, $ 
		      format = '(i4.3, a, a, f6.2, f6.2, f6.2, f6.2, f6.2)'
          fullwave = wave[*,i] 
          fullspec = flux[*,i] 
          fullivar = fluxivar[*,i] 

          outputfile = outputroot+string(format='(i3.3,a)',i+1)+'.fit'

	  bad = 0
	  nonzero = where(fullivar GT 0.0)
	  if (nonzero[0] EQ -1) then begin
	    splog, 'no good points, all have 0.0 or negative sigma'
            bad = 1

          endif else begin 
               
            minfullwave = min(fullwave[nonzero])
            maxfullwave = max(fullwave[nonzero])

;
;	Use medians to merge red and blue here
;
;	    if (numblue GT 0 AND numred GT 0) then begin
;	       maxblue = max(fullwave[where(bluered EQ 0)])
;	       minred = min(fullwave[where(bluered EQ 1)])
;
;	       if (minred LT maxblue) then begin
;                  bluecross = where(bluered EQ 0 and fullwave GT minred $ 
;                      AND fullivar GT 0.0)
;                  redcross = where(bluered EQ 1 and fullwave LT maxblue $
;                      AND fullivar GT 0.0)
;	          if (redcross[0] NE -1 AND bluecross[0] NE -1) then begin 
;	             djs_iterstat, fullspec[bluecross], median=bluemed, $
;                          sigma=bluesigma
;	             djs_iterstat, fullspec[redcross], median=redmed, $
;                          sigma=redsigma
;
;;	             scale[i] = bluemed/redmed
;	             blueflux[i] = bluemed
;                     redflux[i] = redmed
;
;		     if (bluemed - 0.5*bluesigma LE 0) then scale[i] = 1.0
;		     if (redmed - 0.5*redsigma LE 0) then scale[i] = 1.0
;
;	             
;	             splog, i+1, ' Blue:', bluemed, bluesigma, $
;                      ' Red: ', redmed, redsigma,  ' scale: ', scale[i], $
;		      format = '(i4.3, a, f6.2, f6.2, a, f6.2, f6.2, a, f6.2)'
;	             fullspec[redpix] = fullspec[redpix]*scale[i]
;	             fullivar[redpix] = fullivar[redpix]/(scale[i]^2)
;	          endif
;	       endif
;
;	   endif

;
;	get max and min from good pixels
;

	   if (NOT keyword_set(wavemin)) then begin
	     spotmin = fix((minfullwave - zeropoint)/bin) + 1
	     spotmax = fix((maxfullwave - zeropoint)/bin) 
	     wavemin = spotmin * bin + zeropoint
	     wavemax = spotmax * bin + zeropoint
	     bkptmin = wavemin
	     bkptmax = wavemax
	   endif else begin
	     spotmin = 0
	     bkptmin = minfullwave
	     if (NOT keyword_set(wavemax)) then begin 
	       spotmax = fix((maxfullwave - wavemin)/bin) 
               wavemax = spotmax * bin + wavemin
	       bkptmax = wavemax
             endif else begin
	       spotmax = fix((wavemax - wavemin)/bin)
	       bkptmax = maxfullwave
             endelse
	   endelse

	   npix = spotmax - spotmin + 1
	   nbkpt = fix((bkptmax - bkptmin)/bkptbin) + 1

	   newwave = dindgen(npix)*bin + wavemin

;
;		Fit b-spline to each individual spectrum
;
;           indbkpt = slatec_splinefit(fullwave, fullspec, $
;              eachgroup=1, everyn=3, invvar=fullivar, /silent)


           if (keyword_set(everyn)) then begin 
             everyn = (nfiles + 1)/2 
             bkpt = 0
           endif else bkpt = dindgen(nbkpt)*bkptbin + bkptmin
	

;
;	Need to construct ivar for bspline
;
           

;
;	Using newwave as breakpoints
;		

	   ss = sort(fullwave)
	   fullbkpt = slatec_splinefit(fullwave[ss], fullspec[ss], coeff, $
              nord=nord, eachgroup=1, maxiter=20, upper=3.0, lower=3.0, $
              bkpt=bkpt, everyn=everyn, invvar=fullivar[ss], $
              mask=mask, /silent)

	   if (total(coeff) EQ 0.0) then begin
              splog, 'All b-splines coeffs have been set to ZERO!!!!'
           endif

	   splog, 'Masked ', fix(total(1-mask)), ' of', $
                n_elements(mask), ' pixels'

	   mask[ss] = mask

	   bestguess = fltarr(npix)
	   inside = where(newwave GE bkptmin AND newwave LE bkptmax, numinside)
	   if (inside[0] EQ -1) then $
              message, 'No wavelengths inside breakpoints'

   	   fwave = float(newwave[inside])
           bestguess[inside] = slatec_bvalu(fwave,fullbkpt,coeff)

	   bestivar = bestguess*0.0
	   besterr = bestivar

	   for j=0,nfiles-1 do begin
	     these = where(specnum EQ j)
	     if (these[0] NE -1) then begin
	       inbetween = where(newwave GE min(fullwave[these]) AND $
	                         newwave LE max(fullwave[these]))
	       if (inbetween[0] NE -1) then begin
;
;		let's conserve inverse variance
;
                 totalbefore = total(fullivar[these] * mask[these])
	         result = interpol(fullivar[these] * mask[these], $
                      fullwave[these], newwave[inbetween])

                 conservevariance = totalbefore / total(result) 
	         bestivar[inbetween] = bestivar[inbetween] + $
                   result * conservevariance

	       endif
	     endif
           endfor

	   nonzero = where(bestivar GT 0.0)
	   if (nonzero[0] NE -1) then $
             besterr[nonzero] = 1.0/sqrt(bestivar[nonzero])

          if (keyword_set(display)) then begin
            djs_plot, 10^fullwave, fullspec, ps=3, xr=[3700,4800], yr=[-2,10], $
                  title=string(format='(i4,x,a,5(f7.3))', $
                     plugmap[i].fiberid, plugmap[i].objtype, plugmap[i].mag), $
                   xtitle='\lambda [A]', ytitle='Flux (10^-17 cgs)'
  
             djs_oplot, 10^newwave, bestguess,color='red',ps=10
	    endif

          endelse
	  newhdr = hdr

          sxaddpar, newhdr, 'OBJID', string(format='(5(i))', $
                  plugmap[i].objid)
          sxaddpar, newhdr, 'MAG', string(format='(5(f8.3))', $
                plugmap[i].mag)
          sxaddpar, newhdr, 'RAOBJ', plugmap[i].ra, 'RA (deg) of object'
          sxaddpar, newhdr, 'DECOBJ', plugmap[i].dec, 'DEC (deg) of object'
          sxaddpar, newhdr, 'OBJTYPE', plugmap[i].objtype
;
;	Need to be more clever when multiple plugmap's are used
; 	   This is true throughout this entire routine!
;
          sxaddpar, newhdr, 'XFOCAL', plugmap[i].xfocal
          sxaddpar, newhdr, 'YFOCAL', plugmap[i].yfocal
          sxaddpar, newhdr, 'SPECID', plugmap[i].spectrographId

          sxaddpar, newhdr, 'PRIMTARG', plugmap[i].primtarget
          sxaddpar, newhdr, 'SECTARGE', plugmap[i].sectarget
          sxaddpar, newhdr, 'FIBERID', plugmap[i].fiberId

	  sxaddpar,newhdr, 'NWORDER', 2, 'Linear-log10 coefficients'
	  sxaddpar,newhdr, 'WFITTYPE', 'LOG-LINEAR', $
               'Linear-log10 dispersion'
	  sxaddpar,newhdr, 'COEFF0', wavemin, $
               'center wavelength (log10) of first pixel'
	  sxaddpar,newhdr, 'COEFF1', bin, $
               'log10 dispersion per pixel'
	  sxaddpar, newhdr, 'REDSCAL',scale[i],$
                'Red scaling to match blue overlap', AFTER='EXPTIME'

	  sxaddpar, newhdr, 'NAXIS1', n_elements(bestguess)
	  sxaddpar, newhdr, 'NAXIS2', 2
	  sxaddpar, newhdr, 'WAT0_001', 'system=linear'
	  sxaddpar, newhdr, 'WAT1_001', $
              'wtype=linear label=Wavelength units=Angstroms'
	  sxaddpar, newhdr, 'CRVAL1', wavemin, 'Iraf zero point'
	  sxaddpar, newhdr, 'CD1_1', bin, 'Iraf dispersion'
	  sxaddpar, newhdr, 'CRPIX1', 1, 'Iraf starting pixel'
	  sxaddpar, newhdr, 'CTYPE1', 'LINEAR   ' 
	  sxaddpar, newhdr, 'WCSDIM', 2
	  sxaddpar, newhdr, 'DC-FLAG', 1, 'Log-linear flag'


;
;	This is place holding for a fiber with no counts at all!
;
	  if (bad EQ 1) then begin
            bestguess = fltarr(npix)
            besterr = fltarr(npix)
	  endif

          inff = where(finite(bestguess) EQ 0 OR finite(besterr) EQ 0)
          if (inff[0] NE -1) then begin
            stop
            bestguess[inff] = 0.0
            besterr[inff] = 0.0
          endif

	  writefits, outputfile, [[bestguess],[besterr]], newhdr

          if (keyword_set(display)) then begin
              plot, 10^newwave, bestguess, /xstyle, yr=[-3,10]
              djs_oplot, 10^newwave, besterr, color='red'
          endif

;        endelse 
      endfor

;	plots for scaling, if we've done scaling
;      plot,blueflux,scale, ps=1, ytitle='Scale', xtitle='Blue end flux'
;      plot,scale,ps=10, ytitle='Scale', xtitle='Fiber #'
       

      return
end
	  
	  


	
