;+
; NAME:
;   fitflatwidth
;
; PURPOSE:
;   Fit a traceset to the first-order corrected width of the flat field
;
; CALLING SEQUENCE:
;   widthset = fitflatwidth(flux, fluxivar, ansimage, fibermask, $
;        ncoeff=ncoeff, sigma=sigma)
;
; INPUTS:
;   flux       - flat-field extracted flux
;   fluxivar   - corresponding inverse variance
;   ansimage   - output from extract image which contains parameter values
;   fibermask  - nTrace bit mask, which marks bad fibers
;
; OPTIONAL KEYWORDS:
;   ncoeff     - order of legendre polynomial to apply to width vs. row
;   sigma      - The zeroth order sigma of extracted profile
;
; OUTPUTS:
;   widthset   - a traceset structure containing fitted coefficients
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;    Used to fill flatstruct.widthset, which can then be applied
;     to object extraction (known profile widths)
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   1-Mar-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function fitflatwidth, flux, fluxivar, ansimage, fibermask, $
        ncoeff=ncoeff, sigma=sigma

        if (NOT keyword_set(ncoeff)) then ncoeff = 5
        if (NOT keyword_set(sigma)) then sigma = 1.0

        ntrace = (size(flux,/dimen))[1]
        nrow = (size(flux,/dimen))[0]

        mask = (flux GT 0) * (fluxivar GT 0) 

        good = where(mask)
        widthterm = transpose(ansimage[lindgen(ntrace)*2+1,*])

        width = flux * 0.0
        width[good] = widthterm[good]/flux[good]


        xy2traceset, findgen(nrow) # replicate(1,ntrace), $
             width + sigma, widthset, ncoeff=ncoeff, xmin=xmin, xmax=xmax

        if (keyword_set(fibermask)) then begin
          badflats = where(fibermask)
          if (badflats[0] NE -1) then begin
              widthset.coeff[*,badflats] = 0.0
              widthset.coeff[0,badflats] = sigma
          endif
        endif

        return, widthset
end


