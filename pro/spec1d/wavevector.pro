function wavevector, minfullwave, maxfullwave, binsz=binsz, wavemin=wavemin, $
 zeropoint=zeropoint

   if (NOT keyword_set(zeropoint)) then zeropoint = 3.5d
   if (NOT keyword_set(binsz)) then binsz = 1.0d-4

   if (NOT keyword_set(wavemin)) then begin
      spotmin = long((minfullwave - zeropoint)/binsz) + 1
      spotmax = long((maxfullwave - zeropoint)/binsz)
      wavemin = spotmin * binsz + zeropoint
      wavemax = spotmax * binsz + zeropoint
   endif else begin
      spotmin = 0
      if (NOT keyword_set(wavemax)) then begin
        spotmax = long((maxfullwave - wavemin)/binsz)
        wavemax = spotmax * binsz + wavemin
      endif else spotmax = long((wavemax - wavemin)/binsz)
   endelse

   nfinalpix = spotmax - spotmin + 1
   finalwave = dindgen(nfinalpix) * binsz + wavemin

   return, finalwave
end

