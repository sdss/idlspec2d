;+
; NAME:
;   bin_spectra
;
; PURPOSE:
;   Make binned spectra for the purposes of fitting for the 
;   bandpasses 
;
; CALLING SEQUENCE:
;   bin_spectra, flux, binbounds
;
; INPUTS:
;   flux   - flux in a set of spectra [NPIX,NSPEC]
;   invvar   - inverse variance in a set of spectra [NPIX,NSPEC]
;   binbounds   - boundaries of desired bins
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   binflux  - binned flux
;
; COMMENTS:
;   Not currently binning the variances. Sorry.
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; DATA FILES:
;
; REVISION HISTORY:
;   05-APr-2000  Written by M. Blanton, Fermiland
;-
;------------------------------------------------------------------------------
pro bin_spectra, flux, invvar, binbounds, binflux=binflux

   for i=0, n_elements(binbounds)-2 do begin
     denominator= total(invvar[binbounds[i]:binbounds[i+1]-1,*],1,/double)
     denominator= denominator + (denominator LE 0.0)
     binflux[i,*]= $
       total((flux*invvar)[binbounds[i]:binbounds[i+1]-1,*],1,/double) / $
       denominator
   endfor

end
;------------------------------------------------------------------------------
