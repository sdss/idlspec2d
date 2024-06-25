function quickboxcar, image, invvar, tset=tset, fluxivar=fluxivar

   if (NOT keyword_set(radius)) then radius = 3.0

   traceset2xy, tset, ytemp, xcen

   ; Boxcar extract image - no scattered light correction!
   flux = extract_boxcar(image, xcen, radius=radius)

   ; Boxcar extract for any masked pixels
   mask = extract_boxcar(float(invvar EQ 0), xcen, radius=radius) EQ 0

   ; Boxcar extract variance
   varimg = 0 * invvar
   ipos = where(invvar NE 0)
   if (ipos[0] NE -1) then varimg[ipos] = 1. / invvar[ipos]
   fluxvar = extract_boxcar(varimg, xcen, radius=radius)

   fluxivar = 0 * fluxvar
   ipos = where(fluxvar NE 0)
   fluxivar[ipos] = (1. / fluxvar[ipos]) * mask[ipos]

   return, flux
end
