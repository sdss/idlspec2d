pro bbspec_extract, image, invvar, xnow, flux, fluxivar, basisfile=basisfile

   stime0 = systime(1)

   nfibper = 3 ; number of fibers to extract in each call
   nsmallx = 40 ; number of columns to extract in each call
   nsmally = 80 ; number of rows to extract in each call
   npady = 15 ; number of rows to use as padding

   dims = size(image,/dimens)
   nx = dims[0]
   ny = dims[1]

   imgfile = 'tmp_img.fits'
   psffile = 'tmp_psf.fits'
   fluxfile = 'tmp_flux.fits'

   nhdu = 19 ; ??? should be dynamically determined
   bhdr = headfits(basisfile)
   nfiber = sxpar(bhdr,'NAXIS2')
   if (sxpar(bhdr,'NAXIS1') NE ny) then $
    message, 'Dimensions do not agree between image and PSF model!'
   flux = fltarr(ny,nfiber)
   fluxivar = fltarr(ny,nfiber)
   basis = dblarr(ny,nfiber,nhdu)
   bhdr = ptrarr(nhdu)
   for ihdu=0, nhdu-1 do begin
      basis1 = mrdfits(basisfile,ihdu,bhdr1)
; The below to replace NaNs shouldn't be necessary!???
      ibad = where(finite(basis1) EQ 0, nbad)
      igood = where(finite(basis1) EQ 1)
      if (nbad GT 0) then basis1[ibad] = median(basis1[igood]) ; replace NaNs
      basis[*,*,ihdu] = basis1
      bhdr[ihdu] = ptr_new(bhdr1)
   endfor
   if (total(1-finite(basis)) GT 0) then $
    message, 'NaN values in PSF file'
   ; Replace with the X centroids shifted, and trim to only the first entries
   ; if the PSF is only solved for the first fibers in the first rows
;   basis[*,*,0] = xnow[0:ny-1,0:nfiber-1] ; dimensions don't agree ???

   ; Loop through sub-images, solving for nsmall rows at a time on 1 fiber only
   nstepy = nsmally - 2*npady ; number of rows to step up in each call
   nchunk = ceil((ny - 2*npady)/nstepy)
   for ifiber=0, nfiber-1 do begin
      fib1 = ifiber - (nfibper-1)/2
      fib2 = fib1 + nfibper - 1
      fib1 = fib1 > 0
      fib2 = fib2 < (nfiber-1)
      for ichunk=0, nchunk-1 do begin
         y0 = ichunk * (nsmally - 2*npady)
         y1 = (y0 + nsmally - 1) < (ny-1)
         x0 = fix(median(xnow[y0:y1,ifiber])) - nsmallx/2
         x1 = x0 + nsmallx - 1
         x0 = x0 > 0 ; keep in bounds of image
         x1 = x1 < (nx-1) ; keep in bounds of image
print,ifiber,ichunk,x0,x1,y0,y1

         mwrfits, image[x0:x1,y0:y1], imgfile, /create
         mwrfits, invvar[x0:x1,y0:y1], imgfile

         for ihdu=0, nhdu-1 do begin
            basis1 = basis[y0:y1,fib1:fib2,ihdu]
            bhdr1 = *bhdr[ihdu]
            ; Replace the X and Y positions to refer to the subimage positions
            if (ihdu EQ 0) then basis1 -= x0
            if (ihdu EQ 1) then basis1 -= y0
            sxaddpar, bhdr1, 'NAXIS', 2 ; Does not work!!!???
            sxaddpar, bhdr1, 'NAXIS1', y1-y0+1
            sxaddpar, bhdr1, 'NAXIS2', fib2-fib1+1
            if (ihdu EQ 0) then begin
               sxaddpar, bhdr1, 'NPIX_X', x1-x0+1
               sxaddpar, bhdr1, 'NPIX_Y', y1-y0+1
               sxaddpar, bhdr1, 'NSPEC', fib2-fib1+1
               sxaddpar, bhdr1, 'NFLUX', y1-y0+1
            endif
            mwrfits, basis1, psffile, bhdr1, create=(ihdu EQ 0)
         endfor

         pyfile = djs_filepath('pix2spec.py', root_dir=getenv('BBSPEC_DIR'), $
          subdir='examples')
         spawn, 'python '+pyfile+' -i '+imgfile+' -p '+psffile+' -o '+fluxfile
         flux1 = mrdfits(fluxfile)
         fluxivar1 = mrdfits(fluxfile,1)
; The test for NaNs shouldn't be necessary!???
; This appears to happen if there are no good data points
         ibad = where(finite(flux1) EQ 0 OR finite(fluxivar1) EQ 0), nbad)
         if (nbad GT 0) then begin
            flux1[ibad] = 0
            fluxivar1[ibad] = 0
         endif
; The reform below shouldn't be necessary!???
         fluxivar1 = reform(mrdfits(fluxfile,1),size(flux1,/dimen))
         if (ichunk EQ 0) then trim1 = 0 else trim1 = npady
         if (ichunk EQ nchunk-1) then trim2 = 0 else trim2 = npady
         flux[y0+trim1:y1-trim2,ifiber] = flux1[trim1:y1-y0-trim2,ifiber-fib1]
         fluxivar[y0+trim1:y1-trim2,ifiber] = fluxivar1[trim1:y1-y0-trim2,ifiber-fib1]
      endfor
   endfor

   for ihdu=0, nhdu-1 do ptr_free, bhdr[ihdu]

   splog, 'Time to bbspec = ', systime(1)-stime0, ' seconds'
; Compute and return ymodel ???

   return
end
