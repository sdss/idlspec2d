;+
; NAME:
;   calcscatimage
;
; PURPOSE:
;   Just return smoothed scattered light image
;
; CALLING SEQUENCE:
;   scatfit = calcscatimage(ansimage, yrow, [ nscatbkpts= , nx=, ny= ] )
;
; INPUTS:
;     ansimage  -  Keyword Output from extract_image
;     yrow      -  Array of rows extracted in first pass 
;
; OPTIONAL KEYWORDS:
;     nscatbkpts- Number of break points in B-spline; default to 16
;     nx        - X dimension of output image
;     ny        - Y dimension of output image
;
; OUTPUTS:
;    scatfit    - Image of scattered light from smoothing polynomials
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Chebyshev background terms defined over a domain [0,NX-1]
;   mapped to [-1,1]
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   bspline_valu()
;   bspline_iterfit()
;   fchebyshev()
;
; REVISION HISTORY:
;   29-Sep-2000  Written by S. Burles, FNAL, adopted from fitansimage
;-
;------------------------------------------------------------------------------
function calcscatimage, ansimage, yrow, nscatbkpts=nscatbkpts, nx=nx, ny=ny

   if (N_params() LT 1) then begin
      print, 'Syntax - calcscatimage(ansimage, yrow, nscatbkpts=, nx=, ny= )'
      return, -1
   endif

   dims = size(ansimage, /dimens)
   npoly = dims[0]
   nrows = dims[1]
   if (nrows NE n_elements(yrow)) then $
    message, 'Inconsistent dimensions for ANSIMAGE and YROW'

   if(NOT keyword_set(nscatbkpts)) then nscatbkpts=16 ;scattered light

   fitans = fltarr(npoly,ny)
   xnorm = (2.0*findgen(nx)-nx)/(nx-1)
   ynorm = (2.0*findgen(ny)-ny)/(ny-1)
   ynorm_subsamp = (2.0*yrow-ny)/(ny-1)

   ;---------------------------------------------------
   ;	Now do background terms
   ;	First expand terms into NX x NY image
   ;    Without the step function at halfway

   scatfit = fltarr(nx,ny)
   fitans = fltarr(npoly, ny)
  
   for i=0, npoly-1 do begin
     sset = bspline_iterfit(ynorm_subsamp, ansimage[i,*], nbkpts=nscatbkpts, $
      requiren=2)
     fitans[i,*] = bspline_valu(ynorm, sset)
   endfor

   fullcheb = fchebyshev(xnorm, npoly)
   for i=0,ny-1 do $
	  scatfit[*,i] = fullcheb # fitans[*,i]  

   return, scatfit
end
;------------------------------------------------------------------------------
