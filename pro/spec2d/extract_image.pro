;+
; NAME:
;   extract_image
;
; PURPOSE:
;   Extract the fiber profile flux for an entire image
;
; CALLING SEQUENCE:
;   extract_image(fimage, invvar, xcen, sigma, flux, [error, yrow=yrow,
;              ymodel=ymodel, fscat=fscat, proftype = proftype, 
;              wfixed=wfixed, sigmacor=sigmacor, xcencor=xcencor, mask=mask,
;              nPoly=nPoly, maxIter=maxIter, highrej=highrej, lowrej=lowrej])
;
; INPUTS:
;   fimage     - Image[nCol, nRow]
;   invvar     - Inverse Variance[nCol, nRow]
;   xcen       - Initial guesses for X centers[nRow, nFibers]
;   sigma      - sigma of gaussian profile; default to 1.0 
;                  (scalar or [nFibers] or [nRow, nFibers])
;
; OPTIONAL KEYWORDS:
;   yrow       - long array specifying which rows to extract, default is all
;   proftype   - currently, one can only use 1: Gaussian (scalar)
;   wfixed     - array of 1's and zero's which set which parameters are fixed.
;                e.g. [1] just gaussian's with fixed width sigma
;                     [1, 1] fit gaussian + sigma correction
;                     [1, 0, 1] fit gaussian + center correction
;                     [1, 1, 1] fit gaussian + sigma and center corrections.   
;   sigmacor   - new estimates of sigma, must have second element of wfixed set
;   xcencor    - new estimates of xcen, must have third element of wfixed set
;   mask       - byte mask: 1 is good and 0 is bad [nCol,nRow] 
;   nPoly      - order of chebyshev scattered light background; default to 5
;   maxIter    - maximum number of profile fitting iterations; default to 10
;   highrej    - positive sigma deviation to be rejected (default 5.0)
;   lowrej     - negative sigma deviation to be rejected (default 5.0)
;
; OUTPUTS:
;   flux       - Total extracted flux in each profile [nRowExtract, nFibers]
;
; OPTIONAL OUTPUTS:
;   mask       - modified by setting the values of bad pixels to 0
;   error      - Estimated total error in each profile [nRowExtract, nFibers]
;   ymodel     - model best fit of row[nCol, nRow]
;   fscat      - scattered light contribution in each fiber[nRow, nFibers]
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;    extract_row.pro
;
; REVISION HISTORY:
;    8-Aug-1999  Version 0.0 Scott Burles, Chicago 
;-
;------------------------------------------------------------------------------
pro extract_image, fimage, invvar, xcen, sigma, flux, error, yrow=yrow, $
               ymodel=ymodel, fscat=fscat, proftype = proftype,  $
               wfixed=wfixed, sigmacor=sigmacor, xcencor=xcencor, mask=mask, $
               nPoly=nPoly, maxIter=maxIter, highrej=highrej, lowrej=lowrej 

   ; Need 5 parameters
   if (N_params() LT 5) then begin
      print, 'Syntax - extract_image(fimage, invvar, xcen, sigma, flux, [error,'
      print, ' yrow=yrow, ymodel=ymodel, fscat=fscat, proftype = proftype, '
      print, ' wfixed=wfixed, sigmacor=sigmacor, xcencor=xcencor, mask=mask, '
      print, ' nPoly=nPoly, maxIter=maxIter, highrej=highrej, lowrej=lowrej])'
      return
   endif

;
;	fimage should have [nCol,nRow]
;
   fimagesize = size(fimage)
   if (fimagesize[0] NE 2) then message, 'FIMAGE must be 2 dimensional'

   invvarsize = size(invvar)
   if (invvarsize[0] NE 2) then message, 'INVVAR must be 2 dimensional'

   xcensize = size(xcen)
   if (xcensize[0] NE 2) then message,'XCEN must be 2 dimensional [nRow,nTrace]'

;
;	Check dimensions
;

   nx = fimagesize[1]
   if (invvarsize[1] NE nx) then $
    message, 'Number of cols in FIMAGE and INVVAR must be equal'

   ny = fimagesize[2]
   if (invvarsize[2] NE ny) then $
    message, 'Number of rows in FIMAGE and INVVAR must be equal'
   if (xcensize[1] NE ny) then $
    message, 'Number of cols in xcen must equal number of rows in FIMAGE'

;
;	Xcen should have dimensions [nRows, nTrace]
;
   nTrace = xcensize[2]

;
;	For this procedure, we want to work with transposes:)
;	That is [nTrace, nRow] since we work row by row.
;	But all answers will be returned as [nRow, nTrace]
;

   xcenuse = transpose(xcen)
  
   sigmasize = size(sigma)

   if (sigmasize[0] EQ 0) then begin
      sigma1 = xcenuse*0.0 + sigma
   endif else if (sigmasize[0] EQ 1) then begin
      if (sigmasize[1] EQ nTrace) then $ 
         sigma1 = rebin(sigma,nTrace,ny) $
      else if (sigmasize[1] EQ ny) then $
         sigma1 = transpose(rebin(sigma,ny,nTrace)) $
      else message, 'Number of elements in sigma does not equal nTrace nor nRow'
   endif else if (sigmasize[0] EQ 2) then begin
      if (sigmasize[1] NE ny OR sigmasize[2] NE nTrace) then $
         message, '2d sigma array must have same dimensions as XCEN'
      sigma1 = transpose(sigma)
   endif else message, 'Sigma must be scalar, 1d, or 2d array'

   nRowExtract = ny          ; default first to total number of rows
   if (NOT keyword_set(yrow)) then yrow = lindgen(ny) $
   else begin
       nRowExtract = n_elements(yrow)
       check = where(yrow LT 0 OR yrow GE ny, count)
       if (count GT 0) then $
          message, 'YROW has elements which do not correspond to rows in FIMAGE'
       if (nRowExtract GT ny) then $
          message, 'YROW has more elements than FIMAGE'
   endelse

   if (N_elements(nPoly) EQ 0) then nPoly = 5      ; order of background
   if (NOT keyword_set(maxIter)) then maxIter = 10
   if (NOT keyword_set(highrej)) then highrej = 5.0
   if (NOT keyword_set(lowrej)) then lowrej = 5.0 
   if (NOT keyword_set(wfixed)) then wfixed = [1]  ; Zeroth order term
   if (NOT keyword_set(proftype)) then proftype = 1  ; Gaussian
   if (NOT keyword_set(ymodel)) then ymodel = fltarr(nx,ny) 

   masksize = size(mask)
   if (NOT keyword_set(mask)) then mask = make_array(nx,ny, /byte, value=1) $
      else if (masksize[0] NE 2) then $
         message, 'MASK is not 2 dimensional' $
      else if (masksize[1] NE nx) then $
         message, 'Number of cols in FIMAGE and MASK must be equal' $
      else if (masksize[2] NE ny) then $
         message, 'Number of rows in FIMAGE and MASK must be equal'

   nCoeff = n_elements(wfixed)       ;Number of parameters per fibers

   if (keyword_set(sigmacor)) then $
      if((size(sigmacor))[0] NE 2) then $
         message, 'MASK is not 2 dimensional' $
      else if ((size(sigmacor))[1] NE nRowExtract) then $
         message, 'Number of cols in SIGAMCOR and YROW must be equal' $
      else if ((size(sigmacor))[2] NE nTrace) then $
         message, 'Number of rows in SIGAMCOR and XCEN must be equal'

   nPoly = LONG(nPoly)
   ma = nPoly + nTrace*nCoeff
   maxIter = LONG(maxIter)
   proftype = LONG(proftype)

   ; Allocate memory for C routines
   p = fltarr(ma)         ; diagonal errors
   ans = fltarr(ma)       ; parameter values
   ymodelrow = fltarr(nx)
   fscatrow = fltarr(nTrace)
   lTrace = lindgen(nTrace)

;
;	Prepare Output arrays
;

   flux = fltarr(nRowExtract, nTrace)
   error = fltarr(nRowExtract, nTrace)

;
;	Now loop over each row specified in YROW 
;	and extract with rejection with a call to extract_row
;       Check to see if keywords are set to fill optional arrays
;

   for iy=0, nRowExtract-1 do begin
     cur = yrow[iy]
     xcencurrent = xcenuse[*,cur]
     sigmacur = sigma1[*, cur]
;
;	Check that xcen is sorted in increasing order
;	with separations of at 3 pixels.
;
     check = where(xcencurrent[0:nTrace-1] GE xcencurrent[1:nTrace-2] - 3,count)
     if (count GT 0) then $
        message, 'XCEN is not sorted or not separated by greater than 3 pixels.'

     ansrow = extract_row(fimage[*,cur], invvar[*,cur], xcencurrent, $
      sigmacur, ymodel=ymodelrow, fscat=fscatrow, proftype=proftype, $
      wfixed=wfixed, mask=mask[*,cur], diagonal=prow, nPoly=nPoly, $
      maxIter=maxIter, highrej=highrej, lowrej=lowrej, calcCovar=0)

     if(proftype EQ 1) then begin
        print, format='($, ".",i4.4)',cur
;       print, 'Analyzing row', cur, '     With Gaussian Profile', wfixed

       fluxrow = ansrow[0,*]
       errorrow = 1.0 / prow[lTrace*nCoeff] ; best estimate we can do
					      ; without covariance matrix

       if(nCoeff GE 2) then begin 	      ; add in symmetric term if present
	  widthrow = ansrow[1,*]
          errorwidth = 1.0 / prow[lTrace*nCoeff + 1]
          fluxrow = fluxrow + widthrow

;
;	Estimate new widths if specified
;
          if(wfixed[1] GT 0 AND keyword_set(sigmacor)) then begin 
;             print, 'Calculating SIGMACOR...'

;
;	 Make a guess at an underestimated error
;
	     safe = where(fluxrow GT 0.0, safecount)
             if (safecount GT 0) then begin
                r = widthrow[safe]/fluxrow[safe]
                rerror = sqrt((errorrow[safe]*r)^2 + errorwidth[safe]^2) / $
                               fluxrow[safe]

;
;		Only take corrections significant at 2 sigma
;
	        check = where(rerror LT 0.50 AND abs(r) LT 0.4, count)
                if(count GT 0) then $
	          sigmacur[safe[check]] = sigmacur[safe[check]] * $ 
                   (r[check]+ 1.0)
             endif
             sigmacor[iy,*] = sigmacur
          endif

       endif

;
;	Estimate new centroids if specified
;
       if(nCoeff GE 3) then  $ ; calculate asymmetric term if present
       if(wfixed[2] GT 0 AND keyword_set(xcencor)) then begin 
;          print, 'Calculating XCENCOR...'

	 centerrow = ansrow[2,*]
         errorcent = 1.0/prow(lTrace*nCoeff + 2)
;
;	 Make a guess at an underestimated error
;
	 safe = where(fluxrow GT 0.0, safecount)
         if (safecount GT 0) then begin
            r = centerrow(safe)/fluxrow(safe)
            rerror = sqrt((errorrow(safe)*r)^2 + errorcent(safe)^2) / $
                      fluxrow(safe) 
;
;		Only take corrections significant at 2 sigma
;
	    check = where(rerror LT 0.50 AND abs(r) LT 0.4, count)
            if(count GT 0) then $
	       xcencurrent(safe(check)) = xcencurrent(safe(check)) + $ 
                 r(check) * sigmacur(safe(check))
          endif
          xcencor[iy,*] = xcencurrent
       endif
     endif

     if(keyword_set(ymodel)) then ymodel[*,cur] = ymodelrow
     if(keyword_set(fscat)) then fscat[iy,*] = fscatrow

     flux[iy,*] = fluxrow 
     error[iy,*] = errorrow 
   endfor	  

   return
end
;------------------------------------------------------------------------------
