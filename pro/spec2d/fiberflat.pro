;+
; NAME:
;   fiberflat
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   fflat = fiberflat( flat_flux, flat_fluxivar, $
;    [ fibermask, plugmap=, bkspace=, nord=, lower=, upper= ] )
;
; INPUTS:
;   flat_flux  - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   flat_fluxivar - Inverse variance map for FLAT_FLUX.
;   fibermask  - Mask array from plugmap file
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

function fiberflat, flat_flux, flat_fluxivar, fibermask, plugmap=plugmap, $
 bkspace=bkspace, nord=nord, lower=lower, upper=upper

   dims = size(flat_flux, /dimens)
   ny = dims[0]
   ntrace = dims[1]

   if (N_elements(bkspace) EQ 0) then bkspace = 10
   if (N_elements(nord) EQ 0) then nord = 4
   if (N_elements(lower) EQ 0) then lower = 10
   if (N_elements(upper) EQ 0) then upper = 10
   if (N_elements(fibermask) NE ntrace) then fibermask = bytarr(ntrace) + 1

   ; For each fiber, construct the spline fit through the data

   fflat = fltarr(ny,ntrace)
   yaxis = findgen(ny)
   for i=0, ntrace-1 do begin
        print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])
        fullbkpt = slatec_splinefit(yaxis, flat_flux[*,i], coeff, $
         invvar=flat_fluxivar[*,i], nord=4, bkspace=bkspace, $
         lower=lower, upper=upper, maxiter=3)
        fflat[*,i] = slatec_bvalu(yaxis, fullbkpt, coeff)

	ugly = where(finite(fflat[*,i]) EQ 0 OR fflat[*,i] LE 0.0, uglyct)
        if (uglyct GT 0) then begin
          splog, 'bad spline in trace ', i, uglyct
          fflat[ugly,i] = 0.0
          if (uglyct GT 10) then fflat[*,i] = flat_flux[*,i]
        endif
   endfor

   ; Divide FFLAT by a global average of all fibers
   flatmean = total(fflat # fibermask)/(total(fibermask)*ny)
   fflat = fflat / flatmean

   fibermed = djs_median(fflat,1)

   plot,fibermed,ps=1
   if (keyword_set(plugmap)) then begin
   ;
   ; fitting out radial dependence
    
      r2 = (plugmap.xFocal^2 + plugmap.yFocal^2)*1.0e-6  ; units m^2
      plot,r2,fibermed,ps=1
     
      radialcoeff = polyfitw(r2,fibermed,fibermask,1,radialfit)

      for i=0,ntrace-1 do fflat[*,i] = fflat[*,i]/radialfit[i]

      splog, radialcoeff
    endif

   return, fflat
end
;------------------------------------------------------------------------------
