function poly_array, npix, npoly

   arr = fltarr(npix, npoly)
   xvec = findgen(npix) / npix
   for i=0, npoly-1 do $
    arr[*,i] = xvec^i

   return, arr
end
