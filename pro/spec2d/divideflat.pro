;+
; NAME:
;   divideflat
;
; PURPOSE:
;   Divide an extracted image with a fiber-flat
;
; CALLING SEQUENCE:
;   divideflat, flux, fluxivar, fflat, [ fibermask=fibermask, minval=minval,
;    /quiet ]
;
; INPUTS:
;   flux       - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   fluxivar   - Inverse variance map for FLUX.
;   fflat      - Array of flat-field flat-field vectors [Nrow,Ntrace]
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   minval     - Minimum value to consider good for flat-field;
;                default to 0.03.
;   quiet      - don't print to splog?
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

pro divideflat, flux, fluxivar, fflat, fibermask=fibermask, minval=minval, $
   quiet=quiet

   dims = size(flux, /dimens)
   npix = dims[0] 
   ntrace = dims[1] 

   if (n_elements(minval) EQ 0) then minval = 0.03
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(ntrace) 

   if (total(size(fluxivar,/dimens) NE dims) NE 0) then $
    message, 'FLUX and FLUXIVAR are not the same dimensions'

   if (total(size(fflat,/dimens) NE dims) NE 0) then $
    message, 'FLUX and FFLAT are not the same dimensions'

   if ((size(fibermask))[1] NE ntrace) then $
    message, 'FLUX and FIBERMASK have different number of fibers'
    
   for itrace=0, ntrace-1 do begin

      ;  Do we really need to reject bad fibers here, does it hurt
      ;  to divide them out anyway???

;      if (fibermask[itrace] AND fibermask_bits('BADFLAT') NE 0) then begin
      if (1) then begin

         ; Find where the flat field vector for this fiber is less than MINVAL
         qgood = fflat[*,itrace] GT minval
         igood = where(qgood, ngood)
         ibad = where(qgood EQ 0, nbad)

         if (ngood GT 0) then begin
            ; Only divide good fflat pixels
            flux[igood,itrace] = $
             flux[igood,itrace] / fflat[igood,itrace]
            fluxivar[igood,itrace] = $
             fluxivar[igood,itrace] * fflat[igood,itrace]^2
         endif else begin
            if keyword_set(quiet) EQ 0 then $
              splog, 'No good flat field points in trace ', itrace
         endelse

         if (nbad GT 0) then begin
            flux[ibad,itrace] = 0.0
            fluxivar[ibad,itrace] = 0.0
            if (nbad GT 200 AND (keyword_set(quiet) EQ 0)) then $
             splog, 'WARNING: Reject ', nbad, ' low points in trace ', itrace
         endif

      endif else begin ; BAD FIBER

         flux[*,itrace] = 0.0
         fluxivar[*,itrace] = 0.0

      endelse

    endfor
end
;------------------------------------------------------------------------------
