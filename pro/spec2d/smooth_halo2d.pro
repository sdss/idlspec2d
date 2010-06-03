;+
; NAME:
;   smooth_halo
;
; PURPOSE:
;   Given a model image of the CCD, smooth and return a halo image
;
; CALLING SEQUENCE:
;   smooth  = smooth_halo2d ( image, wset, [kernel_size= ] )
;
; INPUTS:
;   image   - Image to smooth, this should be something like ymodel
;                 form extract_image
;   wset    - Wavelength solution 
;
; OPTIONAL KEYWORDS:
;   kernel_size - size of box relative to r0, default = 3.
;                  a value of 5 or higher may be safer.
;
; OUTPUTS:
;   smooth  - Smoothed halo image, same size as image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;    2d version with surface brightness function given by J Gunn
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   traceset2xy
;
; REVISION HISTORY:
;   29-Sep-2000  Original Written by S. Burles, FNAL
;    3-Feb-2005  Finally attempted 2-d convolution, not too slow...
;-
;------------------------------------------------------------------------------
function smooth_halo2d, image, wset, kernel_size=kernel_size

   if NOT keyword_set(kernel_size) then kernel_size=3.

   nx = (size(image))[1]
   ny = (size(image))[2]
   traceset2xy, wset, xx, lam

;
;  convert to microns and just use median wavelength per row
;
   lam = 10^(djs_median(lam,2)-4)


   r0 = 50.0 * lam^2
   f = exp(11.51*(lam - 1.05))

   smooth_model = image * 0.

   for irow = 0,ny-1 do begin
      nr = long(kernel_size * r0[irow]) + 1
      qkernel = fltarr(nr, nr)
      x = findgen(nr) # replicate(1,nr)
      r = sqrt(x^2 + transpose(x^2))

      qBs = f[irow] / (1.0 - f[irow]) / 2.0 / !Pi/ r0[irow] * exp(-r/r0[irow])/ (r > 1)
   
      Bs = transpose([reverse(qbs[1:*,*]) , qbs])
      Bs = [reverse(Bs[1:*,*]), Bs] 

      inc = lindgen(2*nr-1) + irow - (nr-1)
      good = where(inc GE 0 AND inc LT ny, ng)
      for ic =0,ng-1 do $
        smooth_model[*,inc[good[ic]]] =  smooth_model[*,inc[good[ic]]] + $
          convol(image[*,irow], bs[*,good[ic]], /edge_trunc)

     ; Schlegel counter of step number...
;       print, format='("Step ",i5," of ",i5,a1,$)', $
;         irow, ny, string(13b)

   endfor
   return, smooth_model
end

