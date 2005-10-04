;+
; NAME:
;   plotspec_image
;
; PURPOSE:
;   Procedure to display an image, launched from PLOTSPEC
;
; CALLING SEQUENCE:
;   plotspec_image, ra, dec, [ _EXTRA= ]
;
; INPUTS:
;   ra         - Right ascension (degrees)
;   dec        - Declination (degrees)
;
; OPTIONAL INPUTS:
;   _EXTRA     - Keywords to pass to FPBIN_TO_FRAME
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   atv
;   fpbin_to_frame()
;
; REVISION HISTORY:
;   30-Sep-2005  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro plotspec_image, ra, dec, _EXTRA=EXTRA

   fpframe = fpbin_to_frame(ra=ra, dec=dec, rerun='*', $
    _EXTRA=EXTRA, hdr=hdr)
   atv, fpframe, head=hdr

   return
end
;------------------------------------------------------------------------------
