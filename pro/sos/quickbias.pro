function quickbias, biasname, do_lock=do_lock

   if (n_elements(biasname) NE 1) then return, 0

   ;----------
   ; Read in image

   sdssproc, biasname, biasimg, biasivar, color=color, camname=camname, $
    do_lock=do_lock

   ;----------
   ; Test how much of the image was masked by SDSSPROC

   igood = where(biasivar NE 0, ngood)
   fracgood = float(ngood) / n_elements(biasivar)
   if (fracgood LT 0.6) then begin
      splog, 'ABORT: More than 40% of the image is rejected as bad!'
      ptiles = fltarr(101) - 9999.0
      return, create_struct('PERCENTILE', ptiles)
   endif

   ;----------
   ; Compute the percentiles

   ntile = 100
   isort = igood[ sort(biasimg[igood]) ]
   ptiles = biasimg[ isort[ngood*lindgen(ntile)/ntile] ]

   ;----------
   ; Add the 100-th percentile as the maximum unmasked value

   ptiles = [ptiles, biasimg[isort[ngood-1]]]

   ;----------
   ; Return a structure with the percentiles

   rstruct = create_struct('PERCENTILE', float(ptiles))

   return, rstruct
end
