;+
; NAME:
;   extract_image
;
; PURPOSE:
;   Extract the fiber profile flux for an entire image
;
; CALLING SEQUENCE:
;   extract_image(fimage, invvar, xcen, sigma, flux, [finv, yrow=yrow,
;              ymodel=ymodel, fscat=fscat,proftype = proftype,ansimage=ansimage,
;              wfixed=wfixed, mask=mask,
;              nPoly=nPoly, maxIter=maxIter, highrej=highrej, lowrej=lowrej,
;              calcCovar=calcCovar, fitans=fitans, whopping=whopping,relative=relative])
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
;              - or                          2: Exp Cubic
;              - or                          3: Double Gaussian
;              - or              4: Exp Cubic with doublewide Gaussian
;   ansimage   - return the coefficients of fit for each row [nCoeff,nRow]
;   wfixed     - array of 1's and zero's which set which parameters are fixed.
;                e.g. [1] just gaussian's with fixed width sigma
;                     [1, 1] fit gaussian + sigma correction
;                     [1, 0, 1] fit gaussian + center correction
;                     [1, 1, 1] fit gaussian + sigma and center corrections.   
;   mask       - byte mask: 1 is good and 0 is bad [nCol,nRow] 
;   nPoly      - order of chebyshev scattered light background; default to 4
;   maxIter    - maximum number of profile fitting iterations; default to 10
;   highrej    - positive sigma deviation to be rejected (default 10.0)
;   lowrej     - negative sigma deviation to be rejected (default 10.0)
;   calcCovar  - calculate Full covariance matrix
;   fitans     - ratio of profiles to do in single profile fitting
;   relative   - Scale rejection thresholds by reduced chi-squared (default 0)
;   whopping   - traces which have WHOPPINGingly high counts, and need extra
;                background terms
;
; OUTPUTS:
;   flux       - Total extracted flux in each profile [nRowExtract, nFibers]
;
; OPTIONAL OUTPUTS:
;   mask       - modified by setting the values of bad pixels to 0
;   finv       - Estimated inverse variance each profile [nRowExtract, nFibers]
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
pro extract_image, fimage, invvar, xcen, sigma, flux, finv, yrow=yrow, $
               ymodel=ymodel, fscat=fscat,proftype=proftype,ansimage=ansimage, $
               wfixed=wfixed, mask=mask, $
               nPoly=nPoly, maxIter=maxIter, highrej=highrej, lowrej=lowrej, $
	       calcCovar=calcCovar, fitans=fitans, whopping=whopping, $
               relative=relative

   ; Need 5 parameters
   if (N_params() LT 5) then begin
      print, 'Syntax - extract_image(fimage, invvar, xcen, sigma, flux, [finv,'
      print, ' yrow=yrow, ymodel=ymodel, fscat=fscat, proftype = proftype, '
      print, ' ansimage = ansimage, calcCovar=calcCovar, fitans=fitans,relative=relative'
      print, ' wfixed=wfixed, mask=mask, '
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

;   xcenuse = transpose(xcen)
  
   sigmasize = size(sigma)

   if (sigmasize[0] EQ 0) then sigma1 = fltarr(nTrace) + sigma $
   else if (sigmasize[0] EQ 1) then sigma1 = sigma $
   else if (sigmasize[0] EQ 2) then begin
      if (sigmasize[1] NE ny OR sigmasize[2] NE nTrace) then $
         message, '2d sigma array must have same dimensions as XCEN'
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
   if (NOT keyword_set(highrej)) then highrej = 15.0
   if (NOT keyword_set(lowrej)) then lowrej = 20.0 
   if (NOT keyword_set(wfixed)) then wfixed = [1]  ; Zeroth order term
   if (NOT keyword_set(proftype)) then proftype = 1  ; Gaussian
   if (NOT keyword_set(calcCovar)) then calcCovar=0
   if (NOT keyword_set(whopping)) then whopping = -1
   relative = keyword_set(relative)


   if (ARG_PRESENT(ymodel)) then ymodel = fltarr(nx,ny) 

   masksize = size(mask)
   if (NOT keyword_set(mask)) then mask = make_array(nx,ny, /byte, value=1) $
      else if (masksize[0] NE 2) then $
         message, 'MASK is not 2 dimensional' $
      else if (masksize[1] NE nx) then $
         message, 'Number of cols in FIMAGE and MASK must be equal' $
      else if (masksize[2] NE ny) then $
         message, 'Number of rows in FIMAGE and MASK must be equal'

   nCoeff = n_elements(wfixed)       ;Number of parameters per fibers

   nPoly = LONG(nPoly)
   ma = nPoly + nTrace*nCoeff
   maxIter = LONG(maxIter)
   proftype = LONG(proftype)

   ; Allocate memory for C routines
   if (ARG_PRESENT(ansimage)) then $
            ansimage = fltarr(ma,nRowExtract)       ; parameter values
   ymodelrow = fltarr(nx)
   fscatrow = fltarr(nTrace)
   lTrace = lindgen(nTrace)

;
;	Prepare Output arrays
;

   flux = fltarr(nRowExtract, nTrace)
   finv = fltarr(nRowExtract, nTrace)

   whoppingct = 0
   if(whopping[0] NE -1) then $
	 whoppingct = n_elements(whopping)

   ma = nTrace*nCoeff + nPoly + whoppingct

   fullcovar = fltarr(ma,ma)

   squashprofile = 0
   if ARG_PRESENT(fitans) then squashprofile = 1
;
;	Now loop over each row specified in YROW 
;	and extract with rejection with a call to extract_row
;       Check to see if keywords are set to fill optional arrays
;

   ii = where(mask EQ 0, initiallyrejected)

   for iy=0, nRowExtract-1 do begin
     cur = yrow[iy]
     print, format='($, ".",i4.4,a5)',cur,string([8b,8b,8b,8b,8b])

;     xcencurrent = xcenuse[*,cur]
 
     if (sigmasize[0] EQ 2) then  sigmacur = sigma[*, cur] $
     else sigmacur = sigma1
     

     masktemp = mask[*,cur]

     whoppingct = 0
     if(whopping[0] NE -1) then begin
         whoppingcur = transpose(xcen[cur,whopping])
	 whoppingct = n_elements(whopping)
     endif
     
     if ARG_PRESENT(fitans) then begin
          inputans = fitans[0:nTrace*nCoeff-1,cur]
          iback = fitans[nTrace*nCoeff:nTrace*nCoeff+nPoly-1,cur]
     endif

     ansrow = extract_row(fimage[*,cur], invvar[*,cur], $
      xcen[cur,*], sigmacur, ymodel=ymodelrow, fscat=fscatrow, $
      proftype=proftype, iback=iback, $
      wfixed=wfixed, mask=masktemp, diagonal=prow, nPoly=nPoly, $
      oback=oback, niter=niter, squashprofile=squashprofile,inputans=inputans, $
      maxIter=maxIter, highrej=highrej, lowrej=lowrej, calcCovar=calcCovar, $
      whopping=whoppingcur, relative=relative, fullcovar=fullcovar)

     mask[*,cur] = masktemp
     if(ARG_PRESENT(ansimage)) then begin 
       ansimage[0:nTrace*nCoeff-1,iy] = ansrow
       ansimage[nTrace*nCoeff:nTrace*nCoeff+nPoly-1,iy] = oback
     endif

     if(ARG_PRESENT(ymodel)) then ymodel[*,cur] = ymodelrow
     if(ARG_PRESENT(fscat)) then fscat[iy,*] = fscatrow

     calcflux, ansrow, prow, fluxrow, finvrow, wfixed, proftype, lTrace,nCoeff,$
            squashprofile=squashprofile
     flux[iy,*] = fluxrow 
     finv[iy,*] = finvrow
   endfor	  

   ii = where(mask EQ 0, finallyrejected)

   print, 'I masked ', finallyrejected - initiallyrejected, ' pixels'

   ;
   ;	Clean up some memory 

   if (NOT ARG_PRESENT(mask)) then mask = 0
   fullcovar = 0
   inputans = 0
   iback = 0

   return
end
