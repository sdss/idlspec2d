;+
; NAME:
;   fitflatwidth
;
; PURPOSE:
;   Fit a traceset to the first-order corrected width of the flat field
;
; CALLING SEQUENCE:
;   widthset = fitflatwidth(flux, fluxivar, ansimage, [ fibermask, $
;    ncoeff=, sigma= ])
;
; INPUTS:
;   flux       - flat-field extracted flux
;   fluxivar   - corresponding inverse variance
;   ansimage   - output from extract image which contains parameter values
;
; OPTIONAL INPUTS:
;   fibermask  - nTrace bit mask, which marks bad fibers
;   ncoeff     - Order of legendre polynomial to apply to width vs. row;
;                default to 5.
;   sigma      - The zeroth order sigma of extracted profile; default to 1.0.
;
; OUTPUTS:
;   widthset   - Traceset structure containing fitted coefficients
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The widths are forced to be the same as a function of row number
;   for all 16 fibers in each fiber bundle.
;
;   Used to fill flatstruct.widthset, which can then be applied
;   to object extraction (known profile widths).
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   xy2traceset
;
; REVISION HISTORY:
;   01-Mar-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function fitflatwidth, flux, fluxivar, ansimage, fibermask, $
 ncoeff=ncoeff, sigma=sigma

   if (NOT keyword_set(ncoeff)) then ncoeff = 5
   if (NOT keyword_set(sigma)) then sigma = 1.0

   ntrace = (size(flux,/dimen))[1]
   nrow = (size(flux,/dimen))[0]
   if (ntrace NE 320) then $
    message, 'Must have 320 traces'

   ;----------
   ; Generate a mask of good measurements based only upon the fibermask.

   mask = (flux GT 0) * (fluxivar GT 0) 

   if (keyword_set(fibermask)) then begin
      badflats = where(fibermask NE 0)
      if (badflats[0] NE -1) then mask[*,badflats] = 0
   endif

   ;----------
   ; Determine the widths from the output array from EXTRACT_IMAGE.

   igood = where(mask)
   widthterm = transpose(ansimage[lindgen(ntrace)*2+1,*])
   width = flux * 0.0
   width[igood] = widthterm[igood] / flux[igood]

   ;----------
   ; Trigger warning messages if widths are too large.

   medwidth = median(width[igood])
   splog, ((medwidth GT 1.2) ? 'WARNING: ' : '') $
    + 'Median spatial width term = ', medwidth

   ;----------
   ; Perform median across bundles on good arclines only
   ; somewhat tedious, but it works

   width = reform(width,nrow,20,16)
   mask = reform(mask,nrow,20,16)
   width_bundle = fltarr(nrow,16)

   for irow=0, nrow-1 do begin
      for j=0, 15 do begin
         ss = where(mask[irow,*,j])
         if (ss[0] NE -1) then $
          width_bundle[irow,j] = djs_median(width[irow,ss,j])
      endfor
   endfor

   expand_width = (1 + rebin(width_bundle,nrow,ntrace,/sample)) * sigma

   ;----------
   ; Turn the widths back into a traceset.

   xy2traceset, findgen(nrow) # replicate(1,ntrace), $
    expand_width, widthset, ncoeff=ncoeff, xmin=xmin, xmax=xmax

   return, widthset
end
;------------------------------------------------------------------------------
