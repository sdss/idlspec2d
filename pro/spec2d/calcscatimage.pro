;+
; NAME:
;   calcscatimage
;
; PURPOSE:
;   Just return smoothed scattered light image
;
; CALLING SEQUENCE:
;   scatfit = calcscatimage(ansimage, yrow, nscatbkpts=nscatbkpts, $
;            ymin=ymin, ymax=ymax, fullrows=fullrows)
;
; INPUTS:
;     ansimage  -  Keyword Output from extract_image
;     yrow      -  Array of rows extracted in first pass 
;
; OPTIONAL KEYWORDS:
;     ymin      -  lower limit for chebyshev background (default 0)
;     ymax      -  upper limit for chebyshev background (default 2047)
;     fullrows  -  number of rows in full image (default 2048)
;
; OUTPUTS:
;    scatfit    - Image of scattered light from smoothing polynomials
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
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
function calcscatimage, ansimage, yrow, nscatbkpts=nscatbkpts, $
        ymin=ymin, ymax=ymax, fullrows=fullrows

   if (N_params() LT 1) then begin
      print, 'Syntax - calcscatimage(ansimage, yrow, nscatbkpts=nscatbkpts, '
      print, '  ymin=ymin, ymax=ymax, fullrows=fullrows)'
      return, -1
   endif

   nrows = (size(ansimage))[2]
   if n_elements(yrow) EQ 0 then yrow = lindgen(nrows)
   if(NOT keyword_set(nscatbkpts)) then nscatbkpts=40 ;scattered light
   if(NOT keyword_set(ymin)) then ymin=0.0	 
   if(NOT keyword_set(ymax)) then ymax=2047.0	 
   if(NOT keyword_set(fullrows)) then fullrows=2048	 

   npoly = (size(ansimage))[1]
   ynorm = (2.0*yrow-(ymax+ymin))/(ymax-ymin) 
   yfnorm = (2.0*findgen(fullrows)-(ymax+ymin))/(ymax-ymin) 
   fitans = fltarr(npoly,fullrows)


   ;---------------------------------------------------
   ;	Now do background terms
   ;	First expand terms into nrows x nrows image
   ;    Without the step function at halfway

   scatfit = fltarr(fullrows,fullrows)
   fitans = fltarr(npoly, fullrows)
  
   for i=0, npoly-1 do begin
     sset = bspline_iterfit(ynorm, ansimage[i,*], nbkpts=nscatbkpts, $
      requiren=2)
     fitans[i,*] = bspline_valu(yfnorm, sset)
   endfor

   fullcheb = fchebyshev(yfnorm, nPoly)
   for i=0,fullrows-1 do $
	  scatfit[*,i] = fullcheb # fitans[*,i]  

   return, scatfit

end   
	  
