;+
; NAME:
;   trace320cen
;
; PURPOSE:
;   Find the 320 fiber positions for the central row of an image.
;
; CALLING SEQUENCE:
;   xfiber = trace320cen( fimage, [mthresh=, ystart=, nmed=, xgood= ] )
;
; INPUTS:
;   fimage     - Image
;
; OPTIONAL INPUTS:
;   mthresh    - Threshold for peak-finding in convolved row; default to 0.5
;                times the dispersion (found with djs_iterstat).
;   ystart     - Y position in image to search for initial X centers; default
;                to the central row
;   nmed       - Number of rows to median filter around YSTART; default to 21
;
; OUTPUTS:
;   xfiber       - Vector of 320 X centers
;
; OPTIONAL OUTPUTS:
;   xgood        - Set to 1 for fibers that were actually found, 0 otherwise
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
 xgood=xgood

   ; Need 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - xfiber = trace320cen( fimage, [ mthresh=, ystart=, nmed= ]'
      return, -1
   endif

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   if (NOT keyword_set(ystart)) then ystart = ny/2
   if (NOT keyword_set(nmed)) then nmed = 21

   nfiber = 320     ; Number of fibers per side
   npbundle = 20    ; Number of fibers per bundle
   nbun = nfiber / npbundle   ; Number of bundles
   deltax = 6.25

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
   if (NOT keyword_set(mthresh)) then begin
      djs_iterstat, imrow, sigma=sig
      mthresh = 0.5 * sig
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
         endif else if (dx GT 1.0 * deltax AND dx LT 1.6 * deltax) then begin
            ; Try to find break between bundles of 20 by looking
            ; from [1,1.6]*deltax after last center
            xfiber[j] = xpeak[k[i]]
            ibiggap[nbiggap] = j
            nbiggap = nbiggap + 1
print,'BIG GAP',j,xfiber[j], xfiber[j]-xfiber[j-1]
         endif else begin
            xfiber[j] = xfiber[j-1] + deltax
         endelse
      endif
   endfor

   ; Assume that we may have either missing or extra fibers in the beginning.
   ; Find where this first fiber should be relative to the bundle gaps.
   ffiber = fix( median( (ibiggap[0:nbiggap-1]+20) MOD 20 ) )
   if (ffiber LE 10) then xcen = xfiber[ffiber] $ ; up to 10 bogus peaks
    else xcen = xfiber[0] - (20-ffiber) * deltax ; missing up to 10 fibers

   ; Second iteration, where we insist upon the bundle gap positions
   xgood = lonarr(nfiber)
   for ibun=0, nbun-1 do begin
      for jbun=0, npbundle-1 do begin
         n = ibun * npbundle + jbun      ; 0-indexed fiber number

         dx = min(abs(xpeak - xcen), i)
         if (dx LT 0.25 * deltax) then begin
            xfiber[n] = xpeak[i]
            xgood[n] = 1
         endif else begin
            xfiber[n] = xcen
            xgood[n] = 0
         endelse

         xcen = xfiber[n] + deltax
      endfor
      xcen = xcen + 0.28 * deltax ; Add approximate bundle gap
   endfor

;plot,imrow,xr=[0,300],yr=[-50,50],/xst
;djs_oplot,xpeak,xpeak*0+40,ps=1
;djs_oplot,xfiber,xfiber*0+50,ps=1,color='red'
;djs_oplot,xfiber[indgen(16)*20],intarr(16)+50,ps=1,color='blue'
;djs_oplot,conrow,color='green'

   return, xfiber
end
;------------------------------------------------------------------------------
