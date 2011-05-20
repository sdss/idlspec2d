pro bbspec_extract, image, invvar, xnow, flux, fluxivar, basisfile=basisfile

   stime0 = systime(1)

   nsmallx = 40 ; number of columns to extract in each call
   nsmally = 200 ; number of rows to extract in each call
   npady = 20 ; number of rows to use as padding

   dims = size(image,/dimens)
   nx = dims[0]
   ny = dims[1]

   imgfile = 'tmp_img.fits'
   psffile = 'tmp_psf.fits'
   fluxfile = 'tmp_flux.fits'

   nhdu = 18 ; ??? should be dynamically determined
   bhdr = headfits(basisfile)
   nfiber = sxpar(bhdr,'NAXIS2')
   if (sxpar(bhdr,'NAXIS1') NE ny) then $
    message, 'Dimensions do not agree between image and PSF model!'
   flux = fltarr(ny,nfiber)
   fluxivar = fltarr(ny,nfiber)
   basis = dblarr(ny,nfiber,nhdu)
   for ihdu=0, nhdu-1 do $
    basis[*,*,ihdu] = mrdfits(basisfile,ihdu)
   ; Replace with the X centroids shifted, and trim to only the first entries
   ; if the PSF is only solved for the first fibers in the first rows
;   basis[*,*,0] = xnow[0:ny-1,0:nfiber-1] ; ???

   ; Loop through sub-images, solving for nsmall rows at a time on 1 fiber only
   nstepy = nsmally - 2*npady ; number of rows to step up in each call
   nchunk = ceil((ny - 2*npady)/nstepy)
   for ifiber=0, nfiber-1 do begin
      for ichunk=0, nchunk-1 do begin
         y0 = ichunk * (nsmally - 2*npady)
         y1 = (y0 + nsmally - 1) < (ny-1)
         x0 = median(xnow[y0:y1,ifiber]) - nsmallx/2
         x1 = x0 + nsmallx - 1
         x0 = x0 > 0 ; keep in bounds of image
         x1 = x1 < (nx-1) ; keep in bounds of image

         mwrfits, image[x0:x1,y0:y1], imgfile, /create
         mwrfits, invvar[x0:x1,y0:y1], imgfile

         for ihdu=0, nhdu-1 do begin
            basis1 = reform(basis[y0:y1,ifiber,ihdu],y1-y0+1,1)
            ; Replace the X and Y positions to refer to the subimage positions
            if (ihdu EQ 0) then basis1 -= x0
            if (ihdu EQ 1) then basis1 -= y0
            mwrfits, basis1, psffile, create=(ihdu EQ 0)
         endfor

         pyfile = djs_filepath('pix2spec.py', root_dir=getenv('BBSPEC_DIR'), $
          subdir='examples')
         spawn, 'python '+pyfile+' -i '+imgfile+' -p '+psffile+' -o '+fluxfile
         flux1 = mrdfits(fluxfile)
         fluxivar1 = mrdfits(fluxfile,1)
         if (ichunk EQ 0) then trim1 = 0 else trim1 = npady
         if (ichunk EQ nchunk-1) then trim2 = 0 else trim2 = npady
         flux[y0+trim1:y1-trim2,ifiber] = flux1[trim1:y1-y0-trim2]
         fluxivar[y0+trim1:y1-trim2,ifiber] = fluxivar1[trim1:y1-y0-trim2]
      endfor
   endfor

   splog, 'Time to bbspec = ', systime(1)-stime0, ' seconds'

   return
end
