;+
; NAME:
;   fitvacset
;
; PURPOSE:
;   Use measured positions of arc lines and shift coefficients
;    passed through xset to produce a final vacuum wavelength solution 
;
; CALLING SEQUENCE:
;   vacset = fitvacset(xpeak, lambda, wset, xset, ncoeff=ncoeff)
;
; INPUTS:
;   xpeak  - Arc line centroids 
;   lambda - Corresponding wavelengths
;   wset   - Initial arc line solution coefficients
;   xset   - Coefficients specifying shift to coefficients 
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
function fitvacset, xpeak, lambda, wset, xset, ncoeff=ncoeff

      if (NOT keyword_set(ncoeff)) then ncoeff=5

      ;------------------
      ; First convert lambda, and skywaves to log10 vacuum

      splog, 'Converting wavelengths to vacuum'
      vaclambda = lambda
      airtovac, vaclambda
      vacloglam = alog10(vaclambda)

      splog, 'Tweaking to sky lines'

      if (size(xset,/tname) EQ 'UNDEFINED') then begin
         splog, 'WARNING: Sky lines are too noisy! No shifting!'
         xshift = xpeak*0.0
      endif else begin
	 traceset2xy, xset, transpose(xpeak), xshift 
         xshift = transpose(xshift)

         ; Move this QA plot elsewhere ???
         plot, xpeak, xshift, ps=3, xtitle='Arc line position', /ynozero, $
           ytitle = 'Offset to Sky Line [pix]', $
           title = 'Offset to Sky Lines From Wavelength-Solution'

      endelse

      vacset = wset
      nfiber = (size(xpeak))[1]

      xy2traceset, transpose(double(xpeak+xshift)), $
                   vacloglam # (dblarr(nfiber)+1), $
                   vacset, ncoeff=ncoeff, xmin=wset.xmin, xmax=wset.xmax

      return, vacset
end



