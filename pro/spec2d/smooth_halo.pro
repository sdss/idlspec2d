;+
; NAME:
;   smooth_halo
;
; PURPOSE:
;   Given a model image of the CCD, smooth and return a halo image
;
; CALLING SEQUENCE:
;   smooth  = smooth_halo ( image, wset )
;
; INPUTS:
;   image   - Image to smooth, this should be something like ymodel
;                 form extract_image
;   west    - Wavelength solution 
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   smooth  - Smoothed halo image, same size as image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Only does a 1d halo correction, does not spread the halo scattering
;    in the dispersion correction.  
;   Parameters for halo term are hardwired and only qualitatively tested
;    on a few frames.
;
; PROCEDURES CALLED:
;   traceset2xy
;
; REVISION HISTORY:
;   29-Sep-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function halo, image, sigma, exp=exp

   nx = (size(image))[1]
   ny = (size(image))[2]
   nexp  = n_elements(exp)

   ker = gauss_kernel(sigma)
   m = floor(n_elements(ker)/2.)

   imaget = [fltarr(m,nx), transpose(image), fltarr(m,nx)]
   smooth = transpose((convol(imaget, ker))[m:m+ny-1,*])


   for i=0,ny-1 do begin
     imaget = [fltarr(m), smooth[*,i], fltarr(m)]
     smooth[*,i] = (convol(imaget, ker))[m:m+nx-1]

     if keyword_set(exp) then begin
        kere = gauss_kernel(exp[(i < (nexp - 1))], exp=exp)
        me = floor(n_elements(kere)/2.)
        imaget = [fltarr(me), smooth[*,i], fltarr(me)]
        smooth[*,i] = (convol(imaget, kere))[me:me+nx-1]
     endif

   endfor
  
   return, smooth
end 

function smooth_halo, image, wset 

   nx = (size(image))[1]
   ny = (size(image))[2]
   traceset2xy, wset, xx, lam

;
;  convert to microns
;
   lam = 10^(djs_median(lam,2)-4)

   ; smooth = 0.03 * halo(image, 10.0)

   exp = 50.0 * lam^2
   correction = exp(9.0*lam - 9.7)
   smooth2 = halo(image, 4.0, exp=exp)

   ; return, smooth + smooth2 * (correction ## replicate(1,ny))
   return, smooth2 * (correction ## replicate(1,ny))
end

