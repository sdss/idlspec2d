;------------------------------------------------------------------------------
; Compute the frequencies for each bin for an IDL FFT using an array
; with NPIX elements.  Return them in the range (-0.5,0.5].

function fft_wavenums, npix

   if (npix MOD 2 EQ 0) then $
    return, [findgen(npix/2+1),-npix/2.0 + 1.0 + findgen(npix/2-1)]/npix

   return, [findgen(npix/2+1),-(npix- 1.0)/2.0 + findgen(npix/2)]/npix
end
;------------------------------------------------------------------------------
