;+
; NAME:
;   make_skyinvvar
;
; PURPOSE:
;   To reconstruct the contribution of the sky and detector background
;    to the total variance.  Use spPlate*fits files.
;
; CALLING SEQUENCE:
;   skyinvvar = make_skyinvvar(filename)
;
; INPUTS:
;   filename   - Input spPlate file name 
;                     i.e. 'spPlate-0858-52316.fits'
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;   skyinnvar  - Inverse variance as derived from a polynomial fit of 
;                 variance vs. flux (size [npix])
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Return 0's on failure
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   18-Jun-2002   Written by Scott Burles, MIT
;-
;------------------------------------------------------------------------------
function make_skyinvvar, spplatefilename, ncoeff=ncoeff

   if NOT keyword_set(ncoeff) then ncoeff=3

   if (size(spplatefilename,/tname) NE 'STRING') then return, 0
   filename = spplatefilename[0]  

   flux = mrdfits(filename,0,hdr)
   invvar = mrdfits(filename,1)

   skyestimate = 0
   if (size(invvar))[0] NE 2 then return, 0


   npix   = (size(invvar))[1]
   nfiber = (size(invvar))[2]

   x = transpose(flux)
   it = transpose(invvar)
   mask = it GT 0
   y = 1.0/(it + (mask EQ 0))
   xy2traceset, x, y, tset, invvar=mask*it, ncoeff=ncoeff, yfit=yfit
  
   getzeropoint = transpose(fltarr(npix)) 
;
;   Evaluate tset at zero flux to get background variance estimate
;
   traceset2xy, tset, getzeropoint, skyvariance

;
;   Return inverse variance of sky background:
;
   skyinvvar = 1.0/(skyvariance + (skyvariance LE 0)) * (skyvariance GT 0)

   return, skyinvvar
end

