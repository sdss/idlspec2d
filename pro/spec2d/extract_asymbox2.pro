;+
; NAME:
;   extract_asymbox2
;
; PURPOSE:
;   Extract the total flux within a boxcar window at many positions.
;   This routine will accept an asymmetric/variable window
;   Traces are expected to run vertically to be consistent with other
;    extract_  routines
;
;   MODIFIED MODEL WEIGHTING
;
; CALLING SEQUENCE:
;   fextract = extract_asymbox2( image, left, right, $
;                 [ycen, weight_image=weight_image, f_ivar=f_ivar, model=model])
;
; INPUTS:
;   image     - Image
;   left      - Lower boundary of boxcar window (given as floating pt pixels)
;   right     - Upper boundary of boxcar window (given as floating pt pixels)
;
; OPTIONAL KEYWORDS:
;   ycen      - Y positions corresponding to "left" (expected as integers)
;   weight_image -  Weights to be applied to image before boxcar
;
; OUTPUTS:
;   fextract - Extracted flux at positions specified by (left<-->right, ycen)
;
; OPTIONAL OUTPUTS:
;   f_ivar    - the boxcar summed weights, weights=1 if weight_image is missing
;   model     - A model image (of size image)
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:

; REVISION HISTORY:
;   24-Mar-1999  Written by David Schlegel, Princeton.
;   17-Feb-2003  Written with slow IDL routine, S. Burles, MIT
;   27-Jan-2004  Adopted to do asymmetric/varying boxcar
;-
;------------------------------------------------------------------------------
function extract_asymbox2, image, left, right, ycen, $
                  weight_image=weight_image, f_ivar=f_ivar, model=model

   ; Need 2 parameters
   if (N_params() LT 2) then begin
      print, 'Syntax - fextract = extract_idlcar( image, left, right, [ycen,  weight_image=weight_image, f_ivar=f_ivar, model=model] )'
      return, -1
   endif

   if NOT keyword_set(left) then message, 'LEFT is required'
   if NOT keyword_set(right) then message, 'RIGHT is required'
   if (N_elements(LEFT) NE N_elements(RIGHT)) then $
    message, 'Number of elements in LEFT and RIGHT must be equal'

   ndim = (size(left))[0]
   npix = (size(left))[1]
   nTrace =  (ndim EQ 1) ? 1L : (size(left))[2]

   if (NOT keyword_set(ycen)) then begin
     if ndim EQ 1 then ycen = findgen(N_elements(left))
     if ndim EQ 2 then ycen = findgen(npix) # replicate(1,nTrace)
     if ndim GT 2 then message, 'LEFT is not 1 or 2 dimensional'
   endif

   if (N_elements(LEFT) NE N_elements(YCEN)) then $
    message, 'Number of elements in LEFT and YCEN must be equal'

   nx = (size(image))[1]
   ny = (size(image))[2]
   ncen = N_elements(left)

   if ARG_PRESENT(model) then begin
     model = image * 0
     modelwt = image * 0
   endif

   maxwindow = max(right - left)

;   if (min(ycen) LT 0 OR max(ycen) GT ny-y) then $
;    message, 'YCEN contains values out of range'

   tempx = long(maxwindow + 3)
   bigleft =  (left[*]) # replicate(1,tempx)
   bigright =  (right[*]) # replicate(1,tempx)
   spot = dindgen(tempx)## replicate(1,npix*nTrace) + bigleft - 1
   bigy =  ycen[*] # replicate(1,tempx)
   fullspot = (((round(spot+1)-1) > 0) < (nx-1))
   fracleft =  ((fullspot - bigleft) < 0.5) > (-0.5)
   fracright =  ((bigright - fullspot) < 0.5) > (-0.5)
   bigleft = 0.
   bigright = 0.
   weight = (((fracleft+fracright) > 0) < 1)  * $
            ((spot GE -0.5) AND (spot LT (nx-0.5))) * $
            ((bigy GE 0) AND (bigy LE (ny-1)))
   spot = 0.
   fracleft = 0.
   fracright = 0.
   bigy = (temporary(bigy) > 0) < (ny-1)

   if keyword_set(weight_image) then begin
     fextract =total(reform(weight * weight_image[fullspot,bigy] * $
                             image[fullspot,bigy],npix,nTrace,tempx),3)
     f_ivar   =total(reform(weight * weight_image[fullspot,bigy], $
                              npix,nTrace,tempx),3)
     fextract = fextract/(f_ivar + (f_ivar EQ 0)) * (f_ivar GT 0)
   endif else $
     fextract =total(reform(weight * image[fullspot,bigy],npix,nTrace,tempx),3)

   if ARG_PRESENT(model) then begin
       for i=0,tempx-1 do begin
        model[fullspot[*,i],bigy[*,i]] = $
                  model[fullspot[*,i],bigy[*,i]] + fextract*weight[*,i]
        modelwt[fullspot[*,i],bigy[*,i]] = $
                  modelwt[fullspot[*,i],bigy[*,i]] + weight[*,i]
       endfor
     model = (model / (modelwt + (modelwt eq 0.))) * (modelwt gt 0)
   endif

   return, fextract
end
;------------------------------------------------------------------------------

