;+
; NAME:
;   skymask
;
; PURPOSE:
;   Mask regions in spectra where sky-subtraction errors are expected to
;   dominate.
;
; CALLING SEQUENCE:
;   newivar = skymask(objivar, objmask)
;
; INPUTS:
;   objivar    - Inverse variance array [NPIX,NOBJ]
;   objmask    - Mask from spectro-2D outputs, usually the AND-mask [NPIX,NOBJ]
;
; OUTPUTS:
;   newivar    - Modified OBJIVAR, where masked pixels are set to zero
;                [NPIX,NOBJ]
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   pixelmask_bits()
;
; REVISION HISTORY:
;   12-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function skymask, objivar, objmask

   ndim = size(objivar, /n_dimen)
   if (ndim EQ 1) then nobj = 1 $
    else nobj = (size(objivar, /dimens))[1]

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-subtraction
   ; residuals.  Grow that mask by 2 pixels in each direction.

   skymask = (objmask AND pixelmask_bits('BRIGHTSKY')) NE 0
   for iobj=0L, nobj-1 do $
    skymask[*,iobj] = smooth(float(skymask[*,iobj]),5) GT 0
   newivar = objivar * (1 - skymask)

   return, newivar
end
;------------------------------------------------------------------------------
