
function read2dout, filename, wave, hdr=hdr
	spectra = readfits(filename, hdr)	
	if (NOT ARG_PRESENT(wave)) then return, 0

	npix = (size(spectra))[1]
	
	wfit = sxpar(hdr, 'WFITTYPE')
	if (wfit EQ '') then return, 0

	ncoeff = sxpar(hdr, 'NWORDER')
        if (ncoeff EQ 0) then return, 0

	coeff = fltarr(ncoeff)
	for i=0,ncoeff-1 do $
	  coeff[i] = sxpar(hdr, 'COEFF'+strtrim(string(i),2))

	xmin = sxpar(hdr, 'PIXMIN')
	xmax = sxpar(hdr, 'PIXMAX')
	
	pixnorm = (2.d0*findgen(npix) - (xmax + xmin))/(xmax - xmin)

;
;	We should try to store the flegendre(pixnorm,ncoeff) since we'll
;	calculate it so many times (millions and millions)
;
	if (strpos(strupcase(wfit),'LEGENDRE') GT -1) then $
	  wave = flegendre(pixnorm,ncoeff) # coeff
	if (strpos(strupcase(wfit),'CHEBYSHEV') GT -1) then $
	  wave = fchebyshev(pixnorm,ncoeff) # coeff

	return, spectra
end
	  
