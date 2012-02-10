;+
; NAME:
;   filter_solve
;
; PURPOSE:
;   Minimize the chi^2 by varying the filter residuals. Basically solve Ax=b.
;
; CALLING SEQUENCE:
;   filter_solve,binflux,photocounts,binR
;
; INPUTS:
;   binflux  - binned flux
;   photocounts  - psf_band + (fiber_r-psf_r) (not dereddened)
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;
;   binR     - best fit filter parameters
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
pro filter_solve, binflux, photocounts, binR

  ;------------------
  ; Make A and invert
  nbins=(size(binflux))[1] 
  covar=fltarr(nbins,nbins)
  for i = 0, nbins-1 do begin
    for j = 0, nbins-1 do begin
      covar[i,j]=total(binflux[i,*]*binflux[j,*],/double)
    endfor
  endfor
  invcovar=invert(covar,/double)

  ;---------
  ; Make b
  rhs=fltarr(1,nbins)
  for i = 0, nbins-1 do begin
    rhs[0,i]=total(binflux[i,*]*photocounts,/double)
  endfor

  ;---------
  ; Multiply to get binR
  binR=invcovar##rhs

end
;------------------------------------------------------------------------------
