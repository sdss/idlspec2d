; Approximate the SDSS ugriz magnitudes from the Tycho B+R magnitudes.
;------------------------------------------------------------------------------
function tyc_sdssmags, tycdat

   ; Horrible approximations for magnitudes !!!???
   umag = 0.0 + tycdat.vmag + 1.0 * tycdat.bmv
   gmag = -0.149 + tycdat.vmag + 0.626 * tycdat.bmv
   rmag = gmag
   imag = gmag
   zmag = gmag
   magarr = transpose( [[umag],[gmag],[rmag],[imag],[zmag]] )

   return, magarr
end
