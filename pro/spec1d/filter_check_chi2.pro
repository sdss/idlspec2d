;+
; NAME:
;   filter_check_chi2
;
; PURPOSE:
;   
;   Calculate the chi^2 for a particular set of bandpass
;   residuals. 
;
; CALLING SEQUENCE:
;   chi2 = filter_check_chi2(binflux, binR, photocounts)
;
; INPUTS:
;   binflux  - binned flux
;   binR     - parameters to use for flux
;   photocounts  - psf_band + (fiber_r-psf_r) (not dereddened)
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   chi2  - resulting chi2
;
; COMMENTS:
;
; BUGS:
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
pro filter_check_chi2, binflux, binR, photocounts

  fullbinR=(fltarr((size(binflux))[1])+1)##binR
  chi2=total(binflux*fullbinR,/double)-total(photocounts,/double)
 
end
;------------------------------------------------------------------------------
