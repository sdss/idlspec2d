;+
; NAME:
;   fiberflat
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   fflat = fiberflat( flat_flux, flat_fluxivar, $
;    [ bkspace=, nord=, lower=, upper= ] )
;
; INPUTS:
;   flat_flux  - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   flat_fluxivar - Inverse variance map for FLAT_FLUX.
;
; PARAMTERS FOR SLATEC_SPLINEFIT:
;   bkspace
;   nord
;   lower
;   upper
;
; OUTPUTS:
;   fflat      - Array of flat-field flat-field vectors for each fiber
;                that remove relative flat-field variations as a function
;                of wavelength between fibers [Nrow, Ntrace]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The user should first "flat-field" the input array to take out
;   pixel-to-pixel variations.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   14-Oct-1999  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------

function fiberflat, flat_flux, flat_fluxivar, $
 bkspace=bkspace, nord=nord, lower=lower, upper=upper

   dims = size(flat_flux, /dimens)
   ny = dims[0]
   ntrace = dims[1]

   if (N_elements(bkspace) EQ 0) then bkspace = 15
   if (N_elements(nord) EQ 0) then nord = 4
   if (N_elements(lower) EQ 0) then lower = 10
   if (N_elements(upper) EQ 0) then upper = 10

   ; For each fiber, construct the spline fit through the data

   fflat = fltarr(ny,ntrace)
   yaxis = findgen(ny)
   for i=0, ntrace-1 do begin
      print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])
      fullbkpt = slatec_splinefit(yaxis, flat_flux[*,i], coeff, $
       invvar=flat_fluxivar[*,i], nord=4, bkspace=bkspace, $
       lower=lower, upper=upper, maxiter=3)
      fflat[*,i] = slatec_bvalu(yaxis, fullbkpt, coeff)
   endfor

   ; Divide FFLAT by a global average of all fibers
   fflat = fflat / mean(fflat)

   return, fflat
end
;------------------------------------------------------------------------------
