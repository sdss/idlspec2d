function quickbias, biasname

   if (n_elements(biasname) NE 1) then return, 0

   ;----------
   ; Read in image

   sdssproc, biasname, biasimg, biasivar, color=color, camname=camname

   ;----------
   ; Test how much of the image was masked by SDSSPROC

   igood = where(biasivar NE 0, ngood)
   fracgood = float(ngood) / n_elements(biasivar)
   if (fracgood LT 0.9) then begin
      message, 'ABORT: More than 10% of the image is rejected as bad!'
      return, 0
   endif

   ;----------
   ; Compute the percentiles

   ntile = 100
   isort = igood[ sort(biasimg[igood]) ]
   ptiles = biasimg[ isort[ngood*lindgen(ntile)/ntile] ]

   ;----------
   ; Return a structure with the percentiles

   rstruct = create_struct('PERCENTILE', ptiles)

   ;----------
   ; Close splog file

   splog, 'Finished at ', systime()
   splog, /close

   return, rstruct
end
