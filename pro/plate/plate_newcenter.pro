;+
; NAME:
;   plate_newcenter
;
; PURPOSE:
;   Find the new plate center when trying to put an object at the same
;   position relative to the center.
;
; CALLING SEQUENCE:
;   plate_rotate, racen1, deccen1, racen2, deccen2, ra1, dec1, ra2, dec2, $
;    [maxerr= ]
;
; INPUTS:
;   racen1, deccen1  - Center of pointing 1
;   ra1, dec1        - positions on plate 1
;   ra2, dec2        - positions on plate 2
;
; OPTIONAL INPUTS:
;   maxerr           - Maximum error in arcsec; default to 0.01 arcsec
;
; OUTPUTS:
;   racen2, deccen2  - Center of pointing 2
;
; COMMENTS:
;   This routine is complementary to PLATE_ROTATE, solving for a
;   different unknown (the 2nd plate center).
;
; BUGS:
;
; REVISION HISTORY:
;   2001-Nov-27  Written by D. Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro plate_newcenter, racen1, deccen1, racen2, deccen2, ra1, dec1, ra2, dec2, $
 maxerr=maxerr

   if (n_params() EQ 0) then begin 
      print, 'plate_newcenter, racen1, deccen1, racen2, deccen2, ra1, dec1, ra2, dec2'
      return
   endif

   maxerr = 0.01 ; Max error in arcsec
   maxiter = 100

   ; Guess the initial solution -- not valid at high declination!!
   racen2 = (racen1 - ra1 + ra2)[0]
   deccen2 = (deccen1 - dec1 + dec2)[0]

   toterr = 9999.
   iiter = 0
   while (toterr GT maxerr AND iiter LT maxiter) do begin
      iiter = iiter + 1
      ; See how far we are off by
      plate_rotate, racen1[0], deccen1[0], racen2[0], deccen2[0], $
       ra1[0], dec1[0], tmpra, tmpdec
      raerr = ra2[0] - tmpra[0]
      decerr = dec2[0] - tmpdec[0]
      toterr = sqrt(raerr^2 + decerr^2) * 3600 ; in arcsec
      racen2 = racen2 + raerr
      deccen2 = deccen2 + decerr
   endwhile

   return
end
;------------------------------------------------------------------------------
