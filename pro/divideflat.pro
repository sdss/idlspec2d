;+
; NAME:
;   divideflat
;
; PURPOSE:
;   Divide an extracted image with a fiber-flat
;
; CALLING SEQUENCE:
;   divideflat, flux, fluxivar, fflat, [ fibermask=fibermask, minval=minval ]
;
; INPUTS:
;   flux       - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   fluxivar   - Inverse variance map for FLUX.
;   fflat      - Array of flat-field flat-field vectors [Nrow,Ntrace]
;
; OPTIONAL KEYWORDS:
;   fibermask  - Mask of 0 for bad fibers and 1 for good fibers [NFIBER]
;   minval     - Minimum value to consider good for flat-field;
;                default to 0.01.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Wherever the fiber is denoted bad in FIBERMASK, or wherever FFLAT is
;   <= MINVAL, we set FLUX=FLUXIVAR=0.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   splog
;
; REVISION HISTORY:
;   17-Nov-1999  Written by S. Burles, Chicago
;   23-Nov-1999  Modified by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

pro divideflat, flux, fluxivar, fflat, fibermask=fibermask, minval=minval

   dims = size(flux, /dimens)
   npix = dims[0] 
   ntrace = dims[1] 

   if (NOT keyword_set(minval)) then minval = 0.01
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber) + 1

   if (total(size(fluxivar,/dimens) NE dims) NE 0) then $
    message, 'FLUX and FLUXIVAR are not the same dimensions'

   if (total(size(fflux,/dimens) NE dims) NE 0) then $
    message, 'FLUX and FFLAT are not the same dimensions'

   if ((size(fibermask))[1] NE ntrace) then $
    message, 'FLUX and FIBERMASK have different number of fibers'
    
   for itrace=0, ntrace-1 do begin

      if (fibermask[i]) then begin ; GOOD FIBER

         ; Find where the flat field vector for this fiber is less than MINVAL
         qgood = fflat[*,itrace] GT minval
         igood = where(qgood, ngood)
         ibad = where(NOT qgood, nbad)

         if (ngood GT 0) then begin
            ; Only divide good fflat pixels
            flux[igood,itrace] = $
             flux[igood,itrace] / fflat[igood,itrace]
            fluxivar[igood,itrace] = $
             fluxivar[igood,itrace] * fflat[igood,itrace]^2
         endif else begin
            splog, 'No good flat field points in trace ', itrace
         endelse

         if (nbad GT 0) then begin
            flux[ibad,ifiber] = 0.0
            fluxivar[ibad,ifiber] = 0.0
            splog, 'Reject ', nbad, ' low points in trace ', itrace
         endif

      endif else begin ; BAD FIBER

         flux[*,itrace] = 0.0
         fluxivar[*,itrace] = 0.0

      endelse

    endfor
end
;------------------------------------------------------------------------------
