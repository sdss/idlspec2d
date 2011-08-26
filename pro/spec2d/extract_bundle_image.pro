;+
; NAME:
;   extract_bundle_image
;
; PURPOSE:
;   Extract the fiber profile flux for an entire image, using
;   extract_bundle_row to chunk it up over bundles.
;
;   Modified from extract_image, which called extract_row.
;
; **SEE DETAILED DOCUMENTATION COMMENTS IN EXTRACT_BUNDLE_ROW.PRO**
; **FOR MORE INFO ON THE CONVERSION FROM EXTRACT_IMAGE/ROW**
;
; CALLING SEQUENCE:
;   extract_bundle_image, fimage, invvar, xcen, sigma, flux, [finv, yrow=,
;              ymodel=, fscat=, proftype=, ansimage=,
;              wfixed=, mask=mask, pixelmask=,  reject=, wsigma=, 
;              nPoly=, maxIter=, highrej=, lowrej=,
;              fitans=, whopping=, /relative, 
;              nband= ])
;
; INPUTS:
;   fimage     - Image [NCOL,NROW]
;   invvar     - Inverse variance [NCOL,NROW]
;   xcen       - Initial guesses for X centers [NROW,NFIBER]
;   sigma      - Input sigma of gaussian profile; default to 1.0 pixels.
;                This can be a scalar, an [NFIBER] vector, or
;                an [NROW,NFIBER] array.
;
; OPTIONAL KEYWORDS (retained):
;   yrow       - List of row numbers (0-indexed) to extract; default to all.
;   mask       - byte mask: 1 is good and 0 is bad [NCOL,NROW] 
;   pixelmask  - bits set due to extraction rejection [NROW,NFIBER]
;   reject     - Array setting rejection threshholds; defaults are set
;                in EXTRACT_BUNDLE_ROW().
;   nPoly      - order of polynomial background, default to 2
;   maxIter    - maximum number of profile fitting iterations; default to 20
;   highrej    - positive sigma deviation to be rejected (default 10.0)
;   lowrej     - negative sigma deviation to be rejected (default 10.0)
;   relative   - Scale rejection thresholds by reduced chi-squared (default 0)
;   oldreject  - ???
;
; OPTIONAL KEYWORDS (deprecated):
;   proftype   - currently, one can only use 1: Gaussian (scalar)
;              - or                          2: Exp Cubic
;              - or                          3: Double Gaussian
;              - or              4: Exp Cubic with doublewide Gaussian
;   wfixed     - array of 1's and zero's which set which parameters are fixed.
;                e.g. [1] just gaussian's with fixed width sigma
;                     [1, 1] fit gaussian + sigma correction
;                     [1, 0, 1] fit gaussian + center correction
;                     [1, 1, 1] fit gaussian + sigma and center corrections.   
;   nband      - band-width of full covariance fiber profile matrix;
;                default to 1.
;   fitans     - ratio of profiles to do in single profile fitting
;   whopping   - traces which have WHOPPINGingly high counts, and need extra
;                background terms
;   wsigma     - sigma width of whopping profile (exponential, default 25)
;
; OUTPUTS:
;   flux       - Total extracted flux in each profile [nRowExtract,NFIBER]
;
; OPTIONAL OUTPUTS (retained):
;   ansimage   - Coefficients of fit for each row [nCoeff,nRow]
;   mask       - Modified by setting the values of bad pixels to 0
;   finv       - Estimated inverse variance each profile [nRowExtract,NFIBER]
;   ymodel     - Model best fit of row [NCOL,NROW]
;   pimage     - ???
;   chisq      - Chi^2 of each row [NROW]
;
; OPTIONAL OUTPUTS (deprecated):
;   ansimage   - Coefficients of fit for each row [nCoeff,nRow]
;   fscat      - Scattered light contribution in each fiber [NROW,NFIBER]
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   calcflux
;   extract_row()
;   pixelmask_bits()
;   splog
;
; REVISION HISTORY:
;   08-Aug-1999  Written by Scott Burles, Chicago 
;   22-Aug-2000  Added banded-matrix possibility 
;   2010 Modified to _bundle_ form by A. Bolton, U. of Utah.
;   2011 Aug: documentation cleanup by A. Bolton, U. of Utah.
;-
;------------------------------------------------------------------------------
pro extract_bundle_image, fimage, invvar, xcen, sigma, flux, finv, yrow=yrow, $
               ymodel=ymodel, fscat=fscat,proftype=proftype,ansimage=ansimage, $
               wfixed=wfixed, mask=mask, pixelmask=pixelmask, reject=reject, $
               nPoly=nPoly, maxIter=maxIter, highrej=highrej, lowrej=lowrej, $
               fitans=fitans, whopping=whopping, oldreject=oldreject, $
               relative=relative, chisq=chisq, wsigma=wsigma, nband=nband, $
               pimage=pimage, nperbun=nperbun, buffsize=buffsize

   ; Need 5 parameters
   if (N_params() LT 5) then begin
      print, 'Syntax - extract_image(fimage, invvar, xcen, sigma, flux, [finv,'
      print, ' yrow=, ymodel=, fscat=, proftype=, '
      print, ' ansimage, fitans=, /relative'
      print, ' wfixed=, mask=, chisq=, reject=, wsigma='
      print, ' nPoly=, maxIter=, highrej=, lowrej= ])'
      return
   endif

;
; fimage should have [nCol,nRow]
;
   fimagesize = size(fimage)
   if (fimagesize[0] NE 2) then message, 'FIMAGE must be 2 dimensional'

   invvarsize = size(invvar)
   if (invvarsize[0] NE 2) then message, 'INVVAR must be 2 dimensional'

   xcensize = size(xcen)
   if (xcensize[0] NE 2) then message,'XCEN must be 2 dimensional [nRow,nTrace]'

   if (NOT keyword_set(wsigma)) then wsigma = 25.0

;
; Check dimensions
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
; Xcen should have dimensions [nRows, nTrace]
;
   nTrace = xcensize[2]

;
; For this procedure, we want to work with transposes:)
; That is [nTrace, nRow] since we work row by row.
; But all answers will be returned as [nRow, nTrace]
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

   if (N_elements(nPoly) EQ 0) then nPoly = 2      ; order of background, per bundle now!
   if (NOT keyword_set(nband)) then nband = 1L
   if (NOT keyword_set(maxIter)) then maxIter = 20
   if (NOT keyword_set(highrej)) then highrej = 15.0
   if (NOT keyword_set(lowrej)) then lowrej = 20.0 
   if (NOT keyword_set(wfixed)) then wfixed = [1]  ; Zeroth order term
   if (NOT keyword_set(proftype)) then proftype = 1  ; Gaussian
   if (NOT keyword_set(whopping)) then whopping = -1
   relative = keyword_set(relative)
; extract_bundle_row keywords:
   if (NOT keyword_set(nperbun)) then nperbun = 20L
   if (NOT keyword_set(buffsize)) then buffsize = 8L


   if (ARG_PRESENT(ymodel)) then ymodel = fltarr(nx,ny) 
   chisq = fltarr(ny) 

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
;   oldma = nPoly + nTrace*nCoeff
oldma = nTrace
   maxIter = LONG(maxIter)
   proftype = LONG(proftype)

   ; Allocate memory for C routines
   if (ARG_PRESENT(ansimage)) then begin
     if ((size(ansimage))[0] NE 2) then $
        ansimage = fltarr(oldma,nRowExtract)  $
     else if ((size(ansimage))[1] NE oldma OR $
        (size(ansimage))[2] NE nRowExtract OR (size(ansimage))[3] NE 4) then $
            ansimage = fltarr(oldma,nRowExtract)       ; parameter values
   endif

   if (ARG_PRESENT(pimage)) then pimage = fltarr(oldma,nRowExtract)
   

   ymodelrow = fltarr(nx)
   fscatrow = fltarr(nTrace)
   lTrace = lindgen(nTrace)

;
; Prepare Output arrays
;
    if ((size(flux))[0] NE 2) then flux = fltarr(nRowExtract, nTrace) $
     else if ((size(flux))[1] NE nRowExtract OR $
        (size(flux))[2] NE nTrace OR (size(flux))[3] NE 4) $
               then flux = fltarr(nRowExtract, nTrace)

    if ((size(finv))[0] NE 2) then finv = fltarr(nRowExtract, nTrace) $
     else if ((size(finv))[1] NE nRowExtract OR $
        (size(finv))[2] NE nTrace OR (size(finv))[3] NE 4) $
               then finv = fltarr(nRowExtract, nTrace)

   whoppingct = 0
   if(whopping[0] NE -1) then $
    whoppingct = n_elements(whopping)

   ma = nTrace*nCoeff + nPoly + whoppingct

   squashprofile = 0
   if ARG_PRESENT(fitans) then squashprofile = 1
;
; Now loop over each row specified in YROW 
; and extract with rejection with a call to extract_row
; Check to see if keywords are set to fill optional arrays
;

   ii = where(mask EQ 0, initiallyrejected)

   print, ' ROW NITER SIG(med) CHI^2'
   for iy=0, nRowExtract-1 do begin
     cur = yrow[iy]

     if (sigmasize[0] EQ 2) then  sigmacur = sigma[cur, *] $
     else sigmacur = sigma1
     

     masktemp = mask[*,cur]

     whoppingct = 0
     if(whopping[0] NE -1) then begin
         whoppingcur = transpose(xcen[cur,whopping])
         whoppingct = n_elements(whopping)
     endif
     
     if ARG_PRESENT(fitans) then begin
          inputans = fitans[0:nTrace*nCoeff-1,cur]
          if ((size(fitans))[1] GT nTrace*nCoeff) then $
            iback = fitans[nTrace*nCoeff:nTrace*nCoeff+nPoly-1,cur]
     endif

     contribution = 0.02 * (1.0 + 1.5*(cur/1200.0)^2)
     pixelmasktemp = 0
     ansrow = extract_bundle_row(fimage[*,cur], invvar[*,cur], $
      xcen[cur,*],sigmacur[*],ymodel=ymodelrow, fscat=fscatrow, $
      proftype=proftype, iback=iback, reject=reject, pixelmask=pixelmasktemp, $
      wfixed=wfixed, mask=masktemp, diagonal=prow, nPoly=nPoly, $
      niter=niter, squashprofile=squashprofile,inputans=inputans, $
      maxIter=maxIter, highrej=highrej, lowrej=lowrej, $
      whopping=whoppingcur, relative=relative, oldreject=oldreject, $
      reducedChi=chisqrow, nband=nband, contribution=contribution, $
      nperbun=nperbun, buffsize=buffsize, skew=skew, kurt=kurt)

     mask[*,cur] = masktemp

     if (total(finite(ansrow) EQ 0) GT 0) then $
       splog, 'ABORT! ansrow has NaNs at row', cur

     if (total(finite(ymodelrow) EQ 0) GT 0) then $
       splog, 'ABORT! ymodelrow has NaNs at row', cur

;     if (total(finite(fscatrow) EQ 0) GT 0) then $
;       splog, 'ABORT! fscatrow has NaNs at row', cur

     if(ARG_PRESENT(ymodel)) then ymodel[*,cur] = ymodelrow
 ;    if(ARG_PRESENT(fscat)) then fscat[iy,*] = fscatrow
     chisq[cur] = chisqrow

 ;    calcflux, ansrow, prow, fluxrow, finvrow, wfixed, proftype, lTrace,nCoeff,$
 ;           pixelmasktemp, squashprofile=squashprofile
 ;    flux[iy,*] = fluxrow 
 ;    finv[iy,*] = finvrow
; Simpler _bundle_ form for flux and its inverse variance:
     flux[iy,*] = ansrow
     finv[iy,*] = prow^2

     if(ARG_PRESENT(ansimage)) then ansimage[*,iy] = ansrow[0:oldma-1]
     if(ARG_PRESENT(pimage)) then pimage[*,iy] = prow[0:oldma-1]

     if(ARG_PRESENT(pixelmask)) then begin

       ;---------------------------------------------------
       ; Take care of extraction bits first
       ;
       pixelmask[cur,*] = pixelmask[cur,*] OR pixelmasktemp


       ;---------------------------------------------------
       ; Now attempt a cross-talk flag for whopping terms
       ;  do we need a flag for regular profiles, or just test the same??
       ;
       if (whoppingct GT 0) then begin

         if (squashprofile) then wflux = ansrow[nTrace+nPoly:nTrace+nPoly+whoppingct-1]$
         else wflux = ansrow[ma-whoppingct:ma-1]

         for ww = 0,whoppingct - 1 do begin
           guessdist = abs(xcen[cur,whopping[ww]] - xcen[cur,*])/wsigma
           guessflux = exp(-guessdist) * (wflux[ww]/wsigma) * (guessdist LT 5.0)
           crosstalk = (guessflux GT 0.5 * abs(fluxrow))
           crosstalk[ww] = 0
           pixelmask[cur,*] = pixelmask[cur,*] OR (pixelmask_bits('CROSSTALK') * crosstalk)
         endfor
       endif

     endif
;     print, cur, niter, djs_median(sigmacur), chisqrow, $
;      string(13b), format='($, ".",i4.4,i4,f8.2,f8.2,a1)'
   endfor

   if total(finite(chisq) EQ 0) GT 0 then $
      message, "There are infinities in extract_image, need to investigate, related to sdss-pr idlspec2d/2229"

   ii = where(mask EQ 0, finallyrejected)

   splog, 'masked ', finallyrejected - initiallyrejected, ' pixels'

   return
end
