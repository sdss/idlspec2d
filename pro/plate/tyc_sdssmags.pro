;------------------------------------------------------------------------------
; Approximate the SDSS ugriz magnitudes from the Tycho B+V magnitudes.
;
; One can very roughly transform from Johnson B+V to SDSS ugriz
; (at least for main sequence stars up to M0) with:
;   u = V - 0.024 + 2.317 * (B-V)  (Npts= 226)
;   g = V - 0.081 + 0.544 * (B-V)  (Npts= 184)
;   r = V + 0.140 - 0.477 * (B-V)  (Npts= 219)
;   i = V + 0.341 - 1.002 * (B-V)  (Npts= 198)
;   z = V + 0.509 - 1.360 * (B-V)  (Npts= 200)
; The above fits were done to the MT secondary standards.  I'll assume
; that Tycho B+V are approximately Johnson.
;------------------------------------------------------------------------------
function tyc_sdssmags, bminusv, vmag

; Horrible approximations for magnitudes !
;   umag = 0.0 + tycdat.vmag + 1.0 * tycdat.bmv
;   gmag = -0.149 + tycdat.vmag + 0.626 * tycdat.bmv

   umag = vmag - 0.024 + 2.317 * bminusv
   gmag = vmag - 0.081 + 0.544 * bminusv
   rmag = vmag + 0.140 - 0.477 * bminusv
   imag = vmag + 0.341 - 1.002 * bminusv
   zmag = vmag + 0.509 - 1.360 * bminusv
   magarr = transpose( [[umag],[gmag],[rmag],[imag],[zmag]] )

   return, magarr
end
;------------------------------------------------------------------------------
