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
        bkptbin = bkptbin, everyn=everyn, display=display

;
;	Set to 50 km/s for now to match 1d
;	Better guess would be 69 km/s
;
	if (NOT keyword_set(bin)) then bin = (50.0/300000.0) / 2.30258

	if (NOT keyword_set(zeropoint)) then zeropoint = 3.5d
	if (NOT keyword_set(nord)) then nord = 2
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
          tempflux = mrdfits(filenames[0], 0, hdr)
	  tempivar = mrdfits(filenames[0], 1)
	  tempplug = mrdfits(filenames[0], 2)
	  tempwset = mrdfits(filenames[0], 3)
          traceset2xy, tempwset, pixnorm, tempwave

          npix     = (size(tempflux))[1]
	  flux     = [flux, tempflux]
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

        for i=0,nfiber - 1 do begin
  
	  if (plugmap[i].objtype EQ 'SKY' AND NOT keyword_set(dosky)) then $
               splog, ' skipping sky on fiber ', i $
          else if (plugmap[i].objtype EQ 'NA') then $
               splog, ' skipping bad fiber ', i $
          else begin
           
          fullwave = wave[*,i] 
          fullspec = flux[*,i] 
          fullivar = fluxivar[*,i] 

          outputfile = outputroot+string(format='(i3.3,a)',i)+'.fit'
;
;	Use medians to merge red and blue here
;
	    scale = 1.0
	    if (numblue GT 0 AND numred GT 0) then begin
               exptime = exptime * 0.5

	       maxblue = max(fullwave[where(bluered EQ 0)])
	       minred = min(fullwave[where(bluered EQ 1)])

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
	             splog, outputfile, ': scaling red by', scale, bluemed/redmed
	             fullspec[redpix] = fullspec[redpix]*scale
	             fullivar[redpix] = fullivar[redpix]/(scale^2)
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
	     bkptmin = min(fullwave)
	     if (NOT keyword_set(wavemax)) then begin 
	       spotmax = fix((max(fullwave) - wavemin)/bin) 
               wavemax = spotmax * bin + wavemin
	       bkptmax = wavemax
             endif else begin
	       spotmax = fix((wavemax - wavemin)/bin)
	       bkptmax = max(fullwave)
             endelse
	   endelse

	   npix = spotmax - spotmin + 1
	   nbkpt = fix((bkptmax - bkptmin)/bkptbin) + 1

	   newwave = dindgen(npix)*bin + wavemin

         if (keyword_set(everyn)) then $
           everyn = (nfiles + 1)/2 $
         else bkpt = dindgen(nbkpt)*bkptbin + bkptmin
	

;
;	Need to construct ivar for bspline
;
	   nonzero = where(fullivar GT 0.0)
	   if (nonzero[0] EQ -1) then begin
	      print, 'no good points, all have 0.0 or negative sigma'
	      return
           endif

;
;	Using newwave as breakpoints
;		

	   stop
	   ss = sort(fullwave)
	   fullbkpt = slatec_splinefit(fullwave[ss], fullspec[ss], coeff, $
              bkpt=bkpt, everyn=everyn, invvar=fullivar[ss], mask=mask, /silent)

	   stop
	   mask[ss] = mask

	   bestguess = fltarr(npix)
	   inside = where(newwave GE bkptmin AND newwave LE bkptmax, numinside)
	   if (inside[0] EQ -1) then $
              message, 'No wavelengths inside breakpoints'

   	   fwave = float(newwave[inside])
           bestguess[inside] = slatec_bvalu(fwave,fullbkpt,coeff)

	   bestivar = bestguess*0.0
	   besterr = bestivar

	   for i=0,nfiles-1 do begin
	     these = where(specnum EQ i)
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

                 conservevariance = total(result) / totalbefore	
	         bestivar[inbetween] = bestivar[inbetween] + result * conservevariance

	       endif
	     endif
           endfor

	   nonzero = where(bestivar GT 0.0)
	   if (nonzero[0] NE -1) then $
             besterr[nonzero] = 1.0/sqrt(bestivar[nonzero])

	   newhdr = hdr
	   sxaddpar,newhdr, 'NWORDER', 2, 'Linear-log10 coefficients'
	   sxaddpar,newhdr, 'WFITTYPE', 'LOG-LINEAR', $
               'Linear-log10 dispersion'
	   sxaddpar,newhdr, 'COEFF0', wavemin, $
               'center wavelength (log10) of first pixel'
	   sxaddpar,newhdr, 'COEFF1', bin, $
               'log10 dispersion per pixel'
	    sxaddpar, newhdr, 'REDSCAL',scale,'Red scaling to match blue overlap', $
                AFTER='EXPTIME'

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


 	    output = [[newwave],[bestguess],[besterr]]
 

	  writefits, outputfile, [[bestguess],[besterr]], newhdr

          if (keyword_set(display)) then begin
            djs_plot, 10^newwave, bestguess, /xstyle, /ystyle
            djs_oplot, 10^newwave, besterr, color='red'
          endif
        endelse 
      endfor

      return
end
	  
	  


	
