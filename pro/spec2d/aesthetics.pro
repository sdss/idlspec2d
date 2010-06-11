;+
; NAME:
;   aesthetics
;
; PURPOSE:
;   Fill in bad values in spectra with something that looks nice.
;
; CALLING SEQUENCE:
;   newflux = aesthetics(flux,invvar, [method=method])
;
; INPUTS:
;   flux   - A 1D spectrum
;   invvar - Inverse variance of flux
;
; OPTIONAL INPUTS:
;   method - A string describing the method for replacing bad values.  The
;            default value is 'traditional' (see COMMENTS).  method may
;            take the values:
;            'traditional' - the default
;            'noconst'     - like traditional but without setting the /const
;                            flag on djs_maskinterp()
;            'mean'        - bad values are replaced by the mean value of
;                            the good values of flux.
;            'nothing'     - Make no changes at all.
;
; OUTPUTS:
;   newflux - flux with bad values replaced by something more reasonable
;
; COMMENTS:
;   This function is intended to modularize the replacement of bad values in
;   spectra created by combine1fiber.  Traditionally bad values were replaced
;   by calling djs_maskinterp() with the /const flag.  This performs poorly on
;   the edges of the spectra because /const forces values on the edges to
;   take on the last non-bad value, which in some cases could be a large
;   fluctuation.
;
;   The default behaviour of this function is the traditional method
;   described above.
;
; PROCEDURES CALLED:
;   djs_maskinterp()
;
; REVISION HISTORY:
;   2010-06-11 Written by B. A. Weaver
;
; VERSION:
;   $Id$
;
;-
FUNCTION aesthetics, flux, invvar, method=method
    ;
    ; If there are no bad points, do nothing
    ;
    badpts = WHERE(invvar EQ 0,nbad)
    IF nbad EQ 0 THEN RETURN, flux
    ;
    ; Set default method
    ;
    IF ~KEYWORD_SET(method) THEN method='traditional'
    CASE method OF
        'traditional': $
            newflux = djs_maskinterp(flux,invvar EQ 0,/const)
        'noconst': $
            newflux = djs_maskinterp(flux,invvar EQ 0)
        'mean': BEGIN
            newflux = flux
            goodpts = WHERE(newivar GT 0, ngood)
            newflux[badpts] = TOTAL(newflux[goodpts])/ngood
            END
        'nothing': $
            newflux = flux
        ELSE: $
            MESSAGE, 'Unknown method: ' + STRING(method)
    ENDCASE
    RETURN, newflux
END
