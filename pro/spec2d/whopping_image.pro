;
;	return a "whopping" image from given extracted counts and xcen
;

;
;  We need a lookup table, otherwise it's way too slow!
;
function whopping_lookup_lorentz, sigma, fraclist

   if NOT keyword_set(fraclist) then fraclist = findgen(100)*0.01
   nfrac = n_elements(fraclist)

   x = findgen(long(10.0*sigma) + 1)
   nx = n_elements(x)
   
   fullx = x # replicate(1,nfrac) + fraclist ## replicate(1,nx)
   rx = reverse(x) # replicate(1,nfrac) + (1.0 - fraclist ## replicate(1,nx))

   lookup = 0.5*!Pi/sigma / ([rx,fullx]^2 / (sigma*sigma) + 1.0)
   return, reverse(lookup)
end

function whopping_lookup, sigma, fraclist

   if NOT keyword_set(fraclist) then fraclist = findgen(100)*0.01
   nfrac = n_elements(fraclist)

   x = findgen(long(5.0*sigma) + 1)
   nx = n_elements(x)
   
   fullx = x # replicate(1,nfrac) + fraclist ## replicate(1,nx)
   rx = reverse(x) # replicate(1,nfrac) + (1.0 - fraclist ## replicate(1,nx))

   lookup = 2.0 / sigma * exp(-[rx,fullx] / sigma)
   return, reverse(lookup)
end


function whopping_image, flux, xcen, sigma=sigma, xsize=xsize, ysize=ysize,$
   lorentz = lorentz 

   if NOT keyword_set(sigma) then sigma=25.0
   if NOT keyword_set(xsize) then xsize=2048
   if NOT keyword_set(ysize) then ysize=2048


   ; first get lookup array

   if keyword_set(lorentz) then lookup = whopping_lookup_lorentz(sigma) $
   else lookup = whopping_lookup(sigma)
   nl = (size(lookup))[1] / 2
   nf = (size(lookup))[2]

   center = long(xcen + xsize) - xsize
   frac = transpose(long((xcen - center) *100))
   ct = transpose(center)
   ft = transpose(flux)


   ntrace = (size(flux))[2]
   npix   = (size(flux))[1]

   b1 = ct - nl + 1 > 0 
   c1 = nl - ct  > 0
   c2 = (xsize + nl - 2 - ct) < 2 * nl - 1
   b2 = (c2 - c1) + b1

   image = fltarr(xsize,ysize)

;
;   I can't think of a better way to do this??

   for j=0,npix - 1 do begin

   print, format='($, ".",i4.4,a5)', j, string([8b,8b,8b,8b,8b])

     for i=0,ntrace - 1 do begin

         image[b1[i,j]:b2[i,j],j]= image[b1[i,j]:b2[i,j],j] + $
            lookup[c1[i,j]:c2[i,j],frac[i,j]] * ft[i,j] * 0.01
     endfor

   endfor


return, image
end

