;+
; NAME:
;   tweaktrace
;
; PURPOSE:
;   Use fitans to tweak trace and sigma
;   This just perturbs the input values of xcen and sigma
;
; CALLING SEQUENCE:
;   tweaktrace, x, sigma, centershift, sigmashift, maxshift=maxshift
;
; INPUTS:
;   x            - the input trace positions
;   sigma        - the input sigma widths of profiles
;   centershift  - the pixel shift in x
;   sigmashift   - the fractional change of sigma
;
; OPTIONAL KEYWORDS:
;   maxshift     - the absolute allowed shift in x (default = 1.0)
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   19-Oct-1999  Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
pro tweaktrace, x, sigma, centershift, sigmashift, maxshift=maxshift

     if (n_params() LT 4) then begin
       print, 'Syntax - tweaktrace, x, sigma, centershift, sigmashift,'
       return
     endif

     if (NOT keyword_set(maxshift)) then maxshift=1.0

;
;	For proftypes 1 and 2, the following represents first order shifts
;	

     splog, 'median shift ', median(centershift)
     splog, 'max shift ', max(centershift)
     splog, 'min shift ', min(centershift)

     x = x - ((centershift < maxshift) > (-maxshift)) * sigma 
     sigma = (sigmashift + 1.0) * sigma 

    return
end
