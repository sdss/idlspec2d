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
;   xy2traceset
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

        if (keyword_set(fibermask)) then begin
          badflats = where(fibermask)
          if (badflats[0] NE -1) then mask[*,badflats] = 0
        endif

        good = where(mask)
        widthterm = transpose(ansimage[lindgen(ntrace)*2+1,*])

        width = flux * 0.0
        width[good] = widthterm[good]/flux[good]

        width = reform(width,nrow,20,16)
        mask = reform(mask,nrow,20,16)

        width_bundle = fltarr(nrow,16)

;
;       Perform median across bundles on good arclines only
;       somewhat tedious, but it works
;
        for i = 0, nrow - 1 do begin
          for j = 0, 15 do begin
             ss = where(mask[i,*,j])
             if (ss[0] NE -1) then $
               width_bundle[i,j] = djs_median(width[i,ss,j])
          endfor
        endfor

        expand_width = (1 + rebin(width_bundle,nrow,320,/sample)) * sigma

        xy2traceset, findgen(nrow) # replicate(1,ntrace), $
           expand_width, widthset, ncoeff=ncoeff, xmin=xmin, xmax=xmax

        return, widthset
end


