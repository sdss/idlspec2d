;+
; NAME:
;   divideflat
;
; PURPOSE:
;   Divide an extracted image with a fiber-flat
;
; CALLING SEQUENCE:
;   divideflat, flux, fflat, [ invvar=, minval=, /quiet ]
;
; INPUTS:
;   flux       - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   fflat      - Array of flat-field flat-field vectors [Nrow,Ntrace]
;
; OPTIONAL KEYWORDS:
;   invvar     - Inverse variance map for FLUX [Nrow,Ntrace]
;   minval     - Minimum value to consider good for flat-field;
;                default to 0.03.
;   quiet      - don't print to splog?
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   We no longer utilize or set FIBERMASK 
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
pro divideflat, flux, fflat, invvar=invvar, minval=minval, quiet=quiet

   if (n_params() NE 2) then $
    message, 'Wrong number of arguments.'

   ndim = size(flux, /n_dimen)
   dims = size(flux, /dimens)
   npix = dims[0]
   if (ndim EQ 1) then ntrace = 1 $
    else ntrace = dims[1] 

   if (n_elements(minval) EQ 0) then minval = 0.03

   if (keyword_set(invvar) $
    AND (total(size(invvar,/dimens) NE dims) NE 0)) then $
    message, 'FLUX and FLUXIVAR are not the same dimensions'

   if (total(size(fflat,/dimens) NE dims) NE 0) then $
    message, 'FLUX and FFLAT are not the same dimensions'

   for itrace=0, ntrace-1 do begin

      ;  Do we really need to reject bad fibers here, does it hurt
      ;  to divide them out anyway???

      if (1) then begin

         ; Find where the flat field vector for this fiber is less than MINVAL
         qgood = fflat[*,itrace] GT minval
         igood = where(qgood, ngood)
         ibad = where(qgood EQ 0, nbad)

         if (ngood GT 0) then begin
            ; Only divide good fflat pixels
            flux[igood,itrace] = $
             flux[igood,itrace] / fflat[igood,itrace]
            if (keyword_set(invvar)) then $
             invvar[igood,itrace] = $
              invvar[igood,itrace] * fflat[igood,itrace]^2
         endif else begin
            if keyword_set(quiet) EQ 0 then $
              splog, 'No good flat field points in trace ', itrace
         endelse

         if (nbad GT 0) then begin
            flux[ibad,itrace] = 0.0
            if (keyword_set(invvar)) then $
             invvar[ibad,itrace] = 0.0
            if (nbad GT 200 AND (keyword_set(quiet) EQ 0)) then $
             splog, 'Warning: Reject ', nbad, ' low points in trace ', itrace
         endif

      endif else begin ; BAD FIBER

         flux[*,itrace] = 0.0
         if (keyword_set(invvar)) then $
          invvar[*,itrace] = 0.0

      endelse

    endfor
end
;------------------------------------------------------------------------------
