;+
; NAME:
;   skymask
;
; PURPOSE:
;   Mask regions in spectra where sky-subtraction errors are expected to
;   dominate.
;
; CALLING SEQUENCE:
;   newivar = skymask(objivar, andmask, [ormask, ngrow= ])
;
; INPUTS:
;   objivar    - Inverse variance array [NPIX,NOBJ]
;   andmask    - Mask from spectro-2D outputs, usually the AND-mask [NPIX,NOBJ]
;
; OPTIONAL INPUTS:
;   ormask     - Optional OR-mask [NPIX,NOBJ]; if specified, then also mask
;                wherever the REDMONSTER bit is set in this mask.
;   ngrow      - Numbers of pixels neighboring masked pixels to reject along
;                each spectrum; default to 2.
;
; OUTPUTS:
;   newivar    - Modified OBJIVAR, where masked pixels are set to zero
;                [NPIX,NOBJ]
;
; COMMENTS:
;   Grow the mask by 2 pixels in each direction for each spectrum.
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
function skymask, objivar, andmask, ormask, ngrow=ngrow

   ndim = size(objivar, /n_dimen)
   if (ndim EQ 1) then nobj = 1 $
    else nobj = (size(objivar, /dimens))[1]
   if (n_elements(ngrow) EQ 0) then ngrow = 2

   ;----------
   ; Construct a bad pixel mask that will be 0 for good points, 1 for bad.
   ; (Declare it of type long for the smoothing trick employed later.)

   badmask = make_array(size=size(objivar), /long)

   ;----------
   ; If the OR-mask is passed, mask wherever the BADSKYCHI or
   ; the REDMONSTER bit is set

   if (keyword_set(ormask)) then begin
      badmask = badmask OR ((ormask AND pixelmask_bits('BADSKYCHI')) NE 0)
      badmask = badmask OR ((ormask AND pixelmask_bits('REDMONSTER')) NE 0)
   endif

   ;----------
   ; If the AND-mask is passed, mask wherever the spectrum appears to
   ; be dominated by the sky in all exposures.

;   if (keyword_set(ormask)) then begin
;      badmask = badmask OR ((andmask AND pixelmask_bits('BRIGHTSKY')) NE 0)
;   endif

   ;----------
   ; Grow the bad pixel mask by NGROW pixels along each spectrum.

   if (ngrow GT 0) then begin
      width = 2*ngrow + 1
      for iobj=0L, nobj-1 do $
       badmask[*,iobj] = smooth(badmask[*,iobj]*width, width, /edge_truncate) GT 0
   endif

   return, objivar * (1 - badmask)
end
;------------------------------------------------------------------------------
