;+
; NAME:
;   fitvacset
;
; PURPOSE:
;   Re-fit the wavelength solution in vacuum wavelengths, and applying shifts
;   to the pixel positions and heliocentric corrections to wavelengths.
;
; CALLING SEQUENCE:
;   vacset = fitvacset(xpeak, lambda, wset, arcshift, helio=helio)
;
; INPUTS:
;   xpeak       - Arc line centroids 
;   lambda      - Corresponding wavelengths
;   wset        - Initial arc line wavelenthh solution; this is only used
;                 to determine the order of the fit (NCOEFF), and the range
;                 of the fit in pixel space (XMIN,XMAX).  The vectors XPEAK,
;                 LAMBDA are used to re-fit the wavelength solution itself.
;   arcshift    - Shifts to apply to arc lines in pix [NROW,NTRACE]
;
; OPTIONAL KEYWORDS:
;   helio       - Heliocentric correction to add to velocities in km/s.
;
; OUTPUTS:
;   vacset      - output wavelength solution which includes shift to
;                 sky lines and conversion to vacuum wavelengths
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   airtovac
;   splog
;   xy2traceset
;
; REVISION HISTORY:
;   20-Jan-2000  Written by S. Burles, Chicago
;   27-Jul-2011  Added CCD discontinuity handling: A. Bolton, Utah
;-
;------------------------------------------------------------------------------
function fitvacset, xpeak, lambda, wset, arcshift, helio=helio, airset=airset

   xmin = wset.xmin
   xmax = wset.xmax
   ncoeff = (size(wset.coeff, /dimens))[0]
   nfiber = (size(xpeak, /dimens))[0]

; Two-phase readout handling (ASB, jul2011):
   if tag_exist(wset, 'xjumpval') then begin
      xjumplo = wset.xjumplo
      xjumphi = wset.xjumphi
      xjumpval = wset.xjumpval
   endif else begin
      xjumplo = 0
      xjumphi = 0
      xjumpval = 0
   endelse

   if (NOT keyword_set(arcshift)) then arcshift = 0 $
    else splog, 'Tweaking to sky lines'

   ;----------
   ; First convert lambda, and skywaves to log10 vacuum

   splog, 'Converting wavelengths to vacuum'
   vaclambda = lambda
   airtovac, vaclambda

   ;----------
   ; Apply heliocentric correction

   if (keyword_set(helio)) then begin
      vaclambda = vaclambda / (1 + helio/299792.458)
   endif

   vacloglam = alog10(vaclambda)

   ;----------
   ; Re-fit the wavelength solution using LAMBDA converted to vacuum
   ; wavelengths, and pixels shifted by ARCSHIFT.

   xy2traceset, transpose(double(xpeak+arcshift)), $
    vacloglam # (dblarr(nfiber)+1), $
    vacset, ncoeff=ncoeff, xmin=xmin, xmax=xmax, $
    xjumplo=xjumplo, xjumphi=xjumphi, xjumpval=xjumpval

   if ARG_PRESENT(airset) then $
     xy2traceset, transpose(double(xpeak+arcshift)), $ 
                  alog10(lambda) # (dblarr(nfiber)+1), $
                  airset, ncoeff=ncoeff, xmin=xmin, xmax=xmax, $
                  xjumplo=xjumplo, xjumphi=xjumphi, xjumpval=xjumpval


   return, vacset
end
;------------------------------------------------------------------------------
