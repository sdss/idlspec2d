;+
; NAME:
;   trace320cen
;
; PURPOSE:
;   Find the 320 fiber positions for the central row of an image.
;
; CALLING SEQUENCE:
;   xfiber = trace320cen( fimage, [mthresh=, ystart=, nmed=, xgood=, $
;          bundlebreakrange=, nextpeakrange=, bundlegap=,
;          numberFibersPerSide=, deltax=, mthreshfactor= ] )
;
; INPUTS:
;   fimage     - Image
;
; OPTIONAL INPUTS:
;   mthresh    - Threshold for peak-finding in convolved row; default to 0.5
;                times the dispersion (found with djs_iterstat).
;   ystart     - Y position in image to search for initial X centers; default
;                to the central row
;   nmed       - Number of rows to median filter around YSTART;
;                default to 21
;   bundlebreakrange - range in which the bundle break is expected in
;                      units of deltax
;   nextpeakrange - range in which the next peak is expected in
;                      units of deltax
;   bundlegap   - approximate gap between bundles
;   numberFibersPerSide - default 320 if not set
;   deltax      - the expected spacing between fibers in pixel units
;   mtheshfactor      - noise factor for getting lines
;
; OUTPUTS:
;   xfiber     - Vector of 320 X centers
;
; OPTIONAL OUTPUTS:
;   xgood      - Set to 1 for fibers that were actually found, 0 otherwise
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_iterstat
;
; REVISION HISTORY:
;   13-Sep-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function trace320cen, fimage, mthresh=mthresh, ystart=ystart, nmed=nmed, $
 xgood=xgood, bundlebreakreange=bundlebreakrange, nextpeakrange=nextpeakrange, $
bundlegap=bundlegap, numberFibersPerSide=nfiber, deltax=deltax, mthreshfactor=mthreshfactor

  if (NOT keyword_set(bundlebreakrange)) then bundlebreakrange=[1,1.6]
  if (NOT keyword_set(nextpeakrange)) then nextpeakrange=[0.35,0.25]
  if (NOT keyword_set(bundlegap)) then bundlegap=0.28

   ; Need 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - xfiber = trace320cen( fimage, [ mthresh=, ystart=, nmed= ]'
      return, -1
   endif

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   if (NOT keyword_set(ystart)) then ystart = ny/2
   if (NOT keyword_set(nmed)) then nmed = 21

  if (NOT keyword_set(nfiber)) then nfiber = 320     ; Number of fibers per side
   npbundle = 20    ; Number of fibers per bundle
   nbun = nfiber / npbundle   ; Number of bundles
   if (NOT keyword_set(deltax)) then deltax = 6.25
   
   ; Make a sub-image copy of the image and error map
   ylo = ystart - (nmed-1)/2 > 0
   yhi = ystart + (nmed-1)/2 < ny-1
   subimg = fimage[*,ylo:yhi]

   ; Make a copy of the image and error map
   if (keyword_set(invvar)) then begin
      subivar = invvar[*,ylo:yhi]
   endif else begin
      subivar = 1.0 / (subimg > 1)
   endelse

   ; Median filter along each column of the subimage
   imrow = djs_median(subimg, 2)

   ; Convolve with a symmetric function
   kern = [-0.25,-0.25,0.25,0.50,0.25,-0.25,-0.25]
   conrow = convol(imrow, kern)

   ; Set threshold based upon the dispersion in the filtered row
   if (NOT keyword_set(mthreshfactor)) then mthreshfactor=0.5
   if (NOT keyword_set(mthresh)) then begin
      djs_iterstat, imrow, sigma=sig
      mthresh = mthreshfactor * sig
   endif

   ; Find all local peaks that are also above MTHRESH in the convolved row
   rderiv = conrow[1:nx-1] - conrow[0:nx-2]
   izero = where( rderiv[0:nx-3] GT 0 AND rderiv[1:nx-2] LE 0 $
    AND conrow[1:nx-2] GT mthresh)
   xpeak = izero + 0.5 + rderiv[izero] / (rderiv[izero] - rderiv[izero+1])

   ; Go through a first iteration of finding 320 centers, but just to
   ; identify the bundle gaps and where exactly the first fiber should be.
   xfiber = fltarr(nfiber)
   ibiggap = intarr(nfiber-1)
   nbiggap = 0
   xfiber[0] = xpeak[0]
   for j=1, nfiber-1 do begin

      ; Default case if no peak found
      xfiber[j] = xfiber[j-1] + deltax

      k = where(xpeak GT xfiber[j-1], ct)
      if (ct GT 0) then begin
         dx = min(xpeak[k] - xfiber[j-1], i)
         if (dx GT 0.75*deltax AND dx LT 1.20 * deltax) then begin
            xfiber[j] = xpeak[k[i]]
         endif else if (dx GT bundlebreakrange[0] * deltax AND dx LT bundlebreakrange[1] * deltax) then begin
            ; Try to find break between bundles of 20 by looking
            ; from [1,1.6]*deltax after last center
            xfiber[j] = xpeak[k[i]]
            ibiggap[nbiggap] = j
            nbiggap = nbiggap + 1
            splog, 'Big gap at fiber=', j, ', x=', xfiber[j], $
             ', dx=', xfiber[j]-xfiber[j-1], format='(a,i4,a,f7.2,a,f5.2)'
         endif else begin
            xfiber[j] = xfiber[j-1] + deltax
         endelse
      endif
   endfor

   ; Assume that we may have either missing or extra fibers in the beginning.
   ; Find where this first fiber should be relative to the bundle gaps.
  ;; ffiber = fix( median( (ibiggap[0:nbiggap-1]+20) MOD 20 ) )
   ;this identifies the last fiber before a gap
  gaps = shift(xfiber,-1)-xfiber
  gaps=gaps[0:n_elements(gaps)-2]
  normgaps=gaps/deltax
  bundlegap=median(normgaps(where(normgaps gt 2.2 and normgaps lt 3.8)))
  ffiber=min(where(abs(normgaps-bundlegap) lt 0.1))
   if (ffiber LE 10) then begin
      xcen = xfiber[ffiber]; up to 10 bogus peaks
      endif  else begin
      ;;trying a robust determination of the location of the "first" fiber
 ;used to be something wrong     xcen = xfiber[0] - (20-ffiber) * deltax ; missing up to 10 fibers
        gaps = shift(xfiber,-1)-xfiber
        gaps=gaps[0:ffiber]
        w=where(gaps/deltax LT 1.5)
        guessdeltax=median(gaps(w))
       xcenGuess = xfiber[ffiber] - 19 * guessdeltax
       xChoices=xfiber[0]+deltax*(findgen(41)-20)
       dum2=min(abs(xChoices-xcenGuess),in)
       xcen=xChoices[in]
    endelse
    
   ; Second iteration, where we insist upon the bundle gap positions
   xgood = lonarr(nfiber)
   for ibun=0, nbun-1 do begin
      for jbun=0, npbundle-1 do begin
         n = ibun * npbundle + jbun      ; 0-indexed fiber number

         dx = min(abs(xpeak - xcen), i)
         ; Be more lenient about finding the next peak where it should be
         ; if we are at a big gap.
         if ((jbun EQ 0 AND dx LT nextpeakrange[0] * deltax) OR $
             (jbun NE 0 AND dx LT nextpeakrange[1] * deltax)) then begin
            xfiber[n] = xpeak[i]
            xgood[n] = 1
         endif else begin
            xfiber[n] = xcen
            xgood[n] = 0
         endelse

         xcen = xfiber[n] + deltax
      endfor
      xcen = xcen + bundlegap * deltax ; Add approximate bundle gap
   endfor

;plot,imrow,xr=[0,300],yr=[-50,50],/xst
;djs_oplot,xpeak,xpeak*0+40,ps=1
;djs_oplot,xfiber,xfiber*0+50,ps=1,color='red'
;djs_oplot,xfiber[indgen(16)*20],intarr(16)+50,ps=1,color='blue'
;djs_oplot,conrow,color='green'

   return, xfiber
end
;------------------------------------------------------------------------------
