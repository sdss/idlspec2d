;+
; NAME:
;   fitansimage
;
; PURPOSE:
;   Convert output from extract_image (ansimage) into a smooth
;    function of fiber number and row
;
; CALLING SEQUENCE:
;    fitans = fitansimage(ansimage, nterms, ntrace, npoly, nfirst, yrow, $
;            tempflux, fluxm = [1,1,0], scatfit=scatfit, $
;            smallflux, fluxm=fluxm, nord=nord, nscatbkpts=nscatbkpts, $
;            ymin=ymin, ymax=ymax, fullrows=fullrows)
;
; INPUTS:
;     ansimage  -  Keyword Output from extract_image
;     nterms    -  number of profile terms
;     ntrace    -  Number of fibers (currently hardwired to 320)
;     npoly     -  Order of chebyshev polynomial in scattered light
;     yrow      -  Array of rows extracted in first pass 
;     smallflux -  Flux extracted in first pass
;
; OPTIONAL KEYWORDS:
;     fluxm    -   Factor profile terms contribute to total flux
;                   (asymmetric terms are set to zero)
;
; OUTPUTS:
;    fitans    -  Smooth version of ansimage (expanded to 2048 rows)
;
; OPTIONAL OUTPUTS:
;    scatfit   - Image of scattered light from smoothing polynomials
;
; COMMENTS:
;	fitansimage takes the output from extract_image, and smooths
;          the corresponding parameters of nfibers and npoly with
;	   functions of order nord 
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   traceset2xy
;   xy2traceset
;
; REVISION HISTORY:
;   ??-Oct-1999  Written by S. Burles, Chicago
;   25-Feb-2000  Modified to be more robust (yeah!)
; 		   but requires exactly 320 fibers (boo) ???
;-
;------------------------------------------------------------------------------
function fitansimage, ansimage, nparams, nfibers, npoly, yrow, $
        smallflux, fluxm=fluxm, nord=nord, nscatbkpts=nscatbkpts, $
        ymin=ymin, ymax=ymax, fullrows=fullrows,  $
        scatfit=scatfit

  if (N_params() LT 6) then begin
      print, 'Syntax - fitansimage(ansimage, nparams, nfibers, npoly, '
      print, '  yrow, tempflux, fluxm=fluxm, nord=nord, nscatbkpts=nscatbkpts, '
      print, '  ymin=ymin, ymax=ymax, fullrows=fullrows)'
      return, -1
   endif

	if(NOT keyword_set(fluxm)) then $
 			fluxm = make_array(nparams,/long,value=1)
	if(NOT keyword_set(nord)) then nord=3	 ;quadratic
	if(NOT keyword_set(nscatbkpts)) then nscatbkpts=20 ;scattered light
	if(NOT keyword_set(ymin)) then ymin=0.0	 
	if(NOT keyword_set(ymax)) then ymax=2047.0	 
	if(NOT keyword_set(fullrows)) then fullrows=2048	 

        nrows = n_elements(yrow)


	ynorm = (2.0*yrow-(ymax+ymin))/(ymax-ymin)
	yfnorm = (2.0*findgen(fullrows)-(ymax+ymin))/(ymax-ymin) 
	fitans = fltarr(nparams*nfibers+npoly,fullrows)
	newflux = fltarr(nfibers,fullrows)
        iTrace = lindgen(nfibers)*nparams
        iParams = lindgen(nparams) 
	iFiber = findgen(nfibers)

	for i=0,nfibers-1 do fitans[i*nparams,*] = 1.0

        fillans = where(smallflux GT 0)
        if (nfibers NE 320 OR fillans[0] EQ -1) then $
	  splog, 'ABORT: ansimage not well defined for traces' $
        else begin

	  answeight = fltarr(nrows, 20, 16)
	  answeight[fillans] = 1.0

          ;  smallans will be used to fit a smooth function for fitans
	  for i=1,nparams-1 do begin

	    smallans = fltarr(nrows, 20, 16)
	    smallans[fillans] = $
              (transpose(ansimage[itrace+i, *]))[fillans] $
                                 / smallflux[fillans]

	    squashweight = total(answeight,2)
            xy2traceset,yrow # replicate(1,16),djs_median(smallans,2),set1, $
               ncoeff=nord, xmin=ymin, xmax=ymax, invvar= squashweight

            x = findgen(fullrows) # replicate(1,16)
            traceset2xy, set1, x, tempans

            x = ((findgen(16) + 0.5)*20) # replicate(1,fullrows)
            xy2traceset,x, transpose(tempans), finalsigmaset, $
             ncoeff=nord+2, xmin=0, xmax=nfibers

            x = findgen(nfibers) # replicate(1,fullrows)
            traceset2xy, finalsigmaset, x, finalshift

            fitans[iTrace+i,*] = finalshift
          endfor
        endelse

      ;--------------------------------------------------------
      ;  Normalize fitans with total flux
      ; 

      for i=0,nfibers - 1 do $
	    newflux[i,*] = fluxm # fitans[iParams+i*nparams,*]

      for j=0,nparams-1 do $
            fitans[j+iTrace,*] = fitans[j+iTrace,*] / newflux 


      ;---------------------------------------------------
      ;	Now do background terms
      ;	First expand terms into nrows x nrows image
      ;    Without the step function at halfway

      scatfit = fltarr(fullrows,fullrows)

      for i=0, nPoly-1 do begin
        fullbkpt = slatec_splinefit(ynorm, ansimage[nfibers*nparams+i,*], $
                     coeff, nbkpts = nscatbkpts)
        fitans[nfibers*nparams+i,*] = slatec_bvalu(yfnorm, fullbkpt, coeff)
      endfor

      fullcheb = fchebyshev(yfnorm, nPoly)
      for i=0,fullrows-1 do $
	  scatfit[*,i] = fullcheb # fitans[nfibers*nparams:*,i]  

      return, fitans

end   
	  
