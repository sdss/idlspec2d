;+
; NAME:
;   fitvacset
;
; PURPOSE:
;   Use measured positions of arc lines and shift coefficients
;    passed through arcshift to produce a final vacuum wavelength solution 
;
; CALLING SEQUENCE:
;   vacset = fitvacset(xpeak, lambda, wset, arcshift, ncoeff=ncoeff)
;
; INPUTS:
;   xpeak  - Arc line centroids 
;   lambda - Corresponding wavelengths
;   wset   - Initial arc line solution coefficients
;   arcshift    - Shifts to apply to arc lines in pix [NROW,NTRACE]
;
; OPTIONAL KEYWORDS:
;   ncoeff - Number of coefficients to fit final wavelength solution (Default 5)
;
; OUTPUTS:
;   vacset - output wavelength solution which includes shift to
;                    sky lines and conversion to vacuum wavelengths
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   traceset2xy
;   xy2traceset
;
; REVISION HISTORY:
;   20-Jan-2000  Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
function fitvacset, xpeak, lambda, wset, arcshift, ncoeff=ncoeff

      if (NOT keyword_set(ncoeff)) then ncoeff=5

      ;------------------
      ; First convert lambda, and skywaves to log10 vacuum

      splog, 'Converting wavelengths to vacuum'
      vaclambda = lambda
      airtovac, vaclambda
      vacloglam = alog10(vaclambda)

      splog, 'Tweaking to sky lines'

      if (NOT keyword_set(arcshift)) then begin
         splog, 'WARNING: Sky lines are too noisy! No shifting!'
         arcshift = 0
      endif

      vacset = wset
      nfiber = (size(xpeak))[1]

      xy2traceset, transpose(double(xpeak+arcshift)), $
                   vacloglam # (dblarr(nfiber)+1), $
                   vacset, ncoeff=ncoeff, xmin=wset.xmin, xmax=wset.xmax

      return, vacset
end
