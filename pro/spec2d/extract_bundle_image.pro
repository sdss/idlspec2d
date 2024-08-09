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
;   extract_bundle_image, fimage, invvar, rdnoise, xcen, sigma, flux,
;              [finv, yrow=,
;              ymodel=, fscat=, proftype=, ansimage=,
;              wfixed=, mask=mask, pixelmask=,  reject=, wsigma=, 
;              nPoly=, maxIter=, highrej=, lowrej=,
;              fitans=, whopping=, /relative, 
;              nband= ])
;
; INPUTS:
;   fimage     - Image [NCOL,NROW]
;   invvar     - Inverse variance [NCOL,NROW]
;   rdnoise    - Read noise image [NCOL,NROW]
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
;   ybkg       - Background term contributing to ymodel [NCOL,NROW]
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
pro extract_bundle_image, fimage, invvar, rdnoise, xcen, sigma, flux, finv, yrow=yrow, $
               ymodel=ymodel, ybkg=ybkg, fscat=fscat,proftype=proftype,ansimage=ansimage, $
               wfixed=wfixed, mask=mask, pixelmask=pixelmask, reject=reject, $
               nPoly=nPoly, maxIter=maxIter, highrej=highrej, lowrej=lowrej, $
               fitans=fitans, whopping=whopping, oldreject=oldreject, $
               relative=relative, chisq=chisq, wsigma=wsigma, nband=nband, $
               pimage=pimage, buffsize=buffsize, chi2pdf=chi2pdf, plottitle=plottitle,$
               use_image_ivar=use_image_ivar,$ ; JG more robust to trace offsets
               nbundles=nbundles, bundlefibers=bundlefibers, obs=obs, outname=outname, $
               debug=debug
               

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

   splog,'sigmasize: ',sigmasize[0]
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

  ; JG:
  ; no extra background per row by default
  ; because we subtract a smooth background in the CCD image 
  ; fitted in the gaps between fiber bundles  
  if (N_elements(nPoly) EQ 0) then nPoly = 0L 

   
   if (NOT keyword_set(nband)) then nband = 1L
   if (NOT keyword_set(maxIter)) then maxIter = 50
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
   
   ybkg = fltarr(nx,ny) 
   
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
   ybkgrow = fltarr(nx)
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


   fit_smooth_background = 1

   if ( fit_smooth_background NE 0 ) then begin
      
      print, "JG : Fit a smooth background in the gaps between the fiber bundles"
; fit a smooth background in the gaps between the fiber bundles
      nbun = nbundles   ; number of bundle
; there are nbun+1 bands (nbuns-1 between bundles + 2 the left and
; right CCD edges)
      mbkg=fltarr(nRowExtract, nbun+1)

; loop on gaps   
      for b=0, nbun do begin
; loop on CCD rows
         for r=0, nRowExtract-1 do begin

; find edges of each gap         
            if b gt 0 and b lt nbun then begin
               tl = total(bundlefibers[0:b-1],/int)-1 ; last trace of bundle on the left of this gap
               th=tl+1          ; first trace on the right of this gap
               xmin=xcen[r,tl]+3*sigma[r,tl]
               xmax=xcen[r,th]-3*sigma[r,th]
            endif else begin 
               if b eq 0 then begin
                  xmax=xcen[r,0]-3*sigma[r,0]
                  xmin=xmax-10
               endif else begin ; b=nbun
                  xmin=xcen[r,nTrace-1]+3*sigma[r,nTrace-1]
                  xmax=xmin+10
               endelse
            endelse
            xmin=ceil(xmin)
            xmax=floor(xmax)
            if xmin lt 0 then xmin=0
            if xmax gt nx-1 then xmax=nx-1
            
; compute weighted mean in gap and row
            yr=yrow[r]
            ;splog, xmin, xmax, nx, yr, total(invvar[*, yr])
            if total(invvar[*, yr]) EQ 0. then continue
            if xmax LT xmin then stop
            sw=total(invvar[xmin:xmax,yr])
            swf=total(invvar[xmin:xmax,yr]*fimage[xmin:xmax,yr])
            if sw gt 0 then mbkg[r,b] = swf/sw
                        
         endfor                 ; end of loop on rows

; median filtering along rows to reduce noise (and erase effet of
; cosmic rays)  
         nbkg = n_elements(mbkg[yrow[0]:yrow[nRowExtract-1],b])
         if nbkg gt 0 then $
             mbkg[yrow[0]:yrow[nRowExtract-1],b] = djs_median(mbkg[yrow[0]:yrow[nRowExtract-1],b], $
                    width=(50<nbkg), boundary='reflect')
      endfor                    ; and of loop on bands between trace bundles
      
; interpolate on fiber bundles (between the gaps), per row (after the
; median filtering)

; loop on bundles
      start = 0
      for b=0, nbun-1 do begin
         nperbun = bundlefibers[b]
         endf = start + nperbun
        ; loop on rows
         for r=0, nRowExtract-1 do begin
            yr=yrow[r]
            ; find edges (mid point between edges of bundles)
            
            if b eq 0 then begin
                x1 = xcen[r,0]-5-3
                x2 = (xcen[r,endf-1]+xcen[r,endf])/2.
            endif else begin
                if b eq nbun - 1 then begin
                    x1=(xcen[r,start-1]+xcen[r,start])/2.
                    x2=xcen[r,nTrace-1]+5+3
                endif else begin
                    x1 = (xcen[r, start-1] + xcen[r, start])/2.
                    x2 = (xcen[r, endf-1] + xcen[r, endf])/2.
                endelse
            endelse
            x1=floor(x1)
            x2=ceil(x2)
; saved values of background on edges
            v1=mbkg[yr,b]
            v2=mbkg[yr,b+1]
; linear interpolation
            x=x1+lindgen(x2-x1+1)
            ybkg[x1:x2,yr]=((x2-x)*v1+(x-x1)*v2)/(x2-x1)
         endfor
         start = endf
      endfor
         
; we simply subtract this background to the image before the row by
; row extraction. we ignore the increase of noise which is 
; small (we averaged about 10*50 pixel values)
   
      print, "JG : Subtract the smooth background"
      fimage_raw = fimage
      fimage = fimage - ybkg 
      
      if keyword_set(outname) and keyword_set(debug) then begin
	     bkgname = repstr(outname, 'spFrame', 'sdProc_bksub')
	     mwrfits_named, fimage,     bkgname, name='FIMG', /create
	     mwrfits_named, fimage_raw, bkgname, name='RAW_FIMG'
	     mwrfits_named, ybkg,       bkgname, name='YBKG'
	     mwrfits_named, mbkg,       bkgname, name='MBKG'
	     mwrfits_named, invvar,     bkgname, name='IVAR'
	     mwrfits_named, xcen,       bkgname, name='XCEN'
	     mwrfits_named, sigma,      bkgname, name='SIGMA'
      endif
; end of test on smoothe background fit
endif

   chi2pdf=fltarr(nRowExtract, nTrace) ; JG, I need this to discard outliers

; Now loop over each row specified in YROW 
; and extract with rejection with a call to extract_row
; Check to see if keywords are set to fill optional arrays
;

   ii = where(mask EQ 0, initiallyrejected)
   sigout = []

   print, ' ROW NITER SIG(med) CHI^2'
   for iy=0, nRowExtract-1 do begin
     cur = yrow[iy]
     
     if (iy MOD 100 ) eq 0 then print, "JG : extracting CCD row", yrow[iy]


;     print, ' JG DEBUG, go directly to one row of interest'
;     if ( cur LT 3018 ) then continue
;     if ( cur GT 3018 ) then message, 'JG STOP HERE'
     

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
     chi2pdf_of_row = 0.
     ansrow = extract_bundle_row(fimage[*,cur], invvar[*,cur], rdnoise[*,cur], $
      xcen[cur,*],sigmacur[*],ymodel=ymodelrow,ybkg=ybkgrow,fscat=fscatrow, $
      proftype=proftype, iback=iback, reject=reject, pixelmask=pixelmasktemp, $
      wfixed=wfixed, mask=masktemp, diagonal=prow, nPoly=nPoly, $
      niter=niter, squashprofile=squashprofile,inputans=inputans, $
      maxIter=maxIter, highrej=highrej, lowrej=lowrej, $
      whopping=whoppingcur, relative=relative, oldreject=oldreject, $
      reducedChi=chisqrow, nband=nband, contribution=contribution, $
      buffsize=buffsize, skew=skew, kurt=kurt, chi2pdf=chi2pdf_of_row, $
      use_image_ivar=use_image_ivar, nbundles=nbundles, bundlefibers=bundlefibers)

     if (iy eq fix(nRowExtract/2)) or (iy eq 1400) or (iy eq 2800) then begin       
       if keyword_set(plottitle) then begin
	       plot_extraction_profiles, xcen[cur, *], sigmacur, fimage[*,cur], ymodelrow, plottitle, nx, iy
       endif
     endif


     if keyword_set(outname) and keyword_set(debug) then begin 
	   FILE_MKDIR, 'extraction'
	   if iy mod 200 eq 0 then begin
   	    extname = repstr(repstr(outname, 'spFrame', 'extraction/extract_'), '.fits', '_'+strtrim(iy,2)+'.prt')
	    forprint, fimage[*,cur], ymodelrow,  /NOCOMMENT, textout=extname, /silent
	    forprint, xcen[cur,*],sigmacur[*], /NOCOMMENT, textout=repstr(extname, 'extraction/extract_', 'extraction/gauss_'), /silent
       endif
     endif
     sigout=[[sigout],[sigmacur[*]]]
     chi2pdf[iy,*]=chi2pdf_of_row[*]
     
     mask[*,cur] = masktemp

     if (total(finite(ansrow) EQ 0) GT 0) then $
       splog, 'ABORT! ansrow has NaNs at row', cur

     if (total(finite(ymodelrow) EQ 0) GT 0) then $
       splog, 'ABORT! ymodelrow has NaNs at row', cur

;     if (total(finite(fscatrow) EQ 0) GT 0) then $
;       splog, 'ABORT! fscatrow has NaNs at row', cur

     if(ARG_PRESENT(ymodel)) then ymodel[*,cur] = ymodelrow
     if(ARG_PRESENT(ybkg)) then ybkg[*,cur] = ybkgrow

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
;      string(13b), format='($, ".",i4.4,i4,f8.2,f8.2,a1)'  ; OK
   endfor
   if keyword_set(outname) and keyword_set(debug) then begin 
	   bkg_subname = repstr(outname, 'spFrame', 'spFrame_bksub')
	   mwrfits_named, flux, bkg_subname,  name='flux', '/create
	   mwrfits_named, finv, bkg_subname, name='ivar'
    
       extractname = repstr(outname, 'spFrame', 'extraction/sigma_')
       mwrfits_named, sigout, extractname, /create
   endif
; JG : look at chi2pdf to detect outliers and mask them out
; use same highrej sigma threshold as for CCD pixel rejection
; but here with normalize the chi2pdf
   for itrace=0, nTrace-1 do begin
      
      good=where(chi2pdf[*,itrace] gt 0,ngood)

      mean_chi2pdf=0.
      rms_chi2pdf=0.
      
      ; value to reject only terrible outliers
      nsig=(8>highrej) ; if highrej is set to a higher value than 8 use it
            
      ; JG :  clipping
      for loop=0,10 do begin
         if (ngood lt 2) then break
         mean_chi2pdf=mean(chi2pdf[good,itrace])
         rms_chi2pdf=sqrt(mean((chi2pdf[good,itrace]-mean_chi2pdf)^2))
         max_chi2pdf = (mean_chi2pdf>1)+nsig*rms_chi2pdf
         ; print,"DEBUG",loop," ",mean_chi2pdf," ", rms_chi2pdf," ",max_chi2pdf
         previous_ngood=ngood
         good=where((chi2pdf[*,itrace] gt 0) and (chi2pdf[*,itrace] lt max_chi2pdf),ngood)
         if previous_ngood eq ngood then break
      endfor
      if (ngood lt 2) then continue

      ; JG : look for outliers and mask them out
      max_chi2pdf = (mean_chi2pdf>1)+nsig*rms_chi2pdf
      bad=where(chi2pdf[*,itrace] gt max_chi2pdf, nbad)
      nfullrej=nbad
      npartrej=0
      if (nbad gt 0) then begin 
         chi2pdf[bad,itrace]=0.       
         finv[bad,itrace]=0.
         if(ARG_PRESENT(pixelmask)) then begin           
            pixelmask[bad,itrace] = pixelmask[bad,itrace] OR pixelmask_bits('FULLREJECT')
         endif
      endif
      if(ARG_PRESENT(pixelmask)) then begin      
         nfullrej = total((pixelmask[*,itrace] AND pixelmask_bits('FULLREJECT')) gt 0)
         npartrej = total((pixelmask[*,itrace] AND pixelmask_bits('PARTIALREJECT')) gt 0)
      endif
      if (ngood gt 1) then splog, "trace",itrace," chi2/ndf=",mean_chi2pdf," rms=",rms_chi2pdf," nbad=", nbad, " nfullrej=",long(nfullrej), " npartrej=",long(npartrej)
      
   endfor
   

   if total(finite(chisq) EQ 0) GT 0 then $
      message, "There are infinities in extract_image, need to investigate, related to sdss-pr idlspec2d/2229"

   ii = where(mask EQ 0, finallyrejected)

   splog, 'masked ', finallyrejected - initiallyrejected, ' pixels'

   return
end
