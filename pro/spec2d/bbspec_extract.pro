;+
; NAME:
;   bbspec_extract
;
; PURPOSE:
;   Run bbspec 2D extraction code
;
; CALLING SEQUENCE:
;   bbspec_extract, image, invvar, flux, fluxivar, basisfile=, $
;    [ ximg=, frange=, yrange=, ymodel=, tmproot= ]
;
; INPUTS:
;   image      - Image [NX,NY]
;   innvar     - Inverse variance image corresponding to IMAGE [NX,NY]
;
; OPTIONAL INPUTS:
;   ximg       - X centroids of fibers on IMAGE, to replace those entries
;                in the PSF file [NY,NFIBER]
;   frange     - Fiber number range (0-indexed); default to all fibers
;                represented in the PSF file
;   yrange     - 0-indexed ange of rows to extract; default to all rows
;                that contain any unmasked pixels
;   tmproot    - Root file name for temporary files; default to 'tmp-';
;                necessary to be a unique string identifier if multiple
;                instances of this procedure running in the same directory,
;                for example those spawned by BBSPECT_TEST.
;
; OUTPUTS:
;   flux       - Extracted flux vectors [NY,NFIBER]
;   fluxivar   - Extracted inverse variance vectors [NY,NFIBER]
;   basisfile  - File with PSF model for bbspec
;
; OPTIONAL OUTPUTS:
;   ymodel     - Model-fit image
;
; COMMENTS:
;   Currently hard-wired to extract in chunks of 20 fibers by 100 rows,
;   with the 15 rows on the top and bottom to be used as padding and
;   discarded.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   splog
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   24-May-2011  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro bbspec_extract_readpsf, basisfile1, pbasis, phdr, nfiber=nfiber

   ; Determine the number of HDUs in the PSF file
   nhdu = 0
   basisfile = findfile(basisfile1, count=ct)
   if (ct EQ 0) then $
    message, 'PSF file not found '+string(basisfile1)
   while (size(headfits(basisfile,exten=nhdu,/silent),/tname) EQ 'STRING') do $
    nhdu++
   if (nhdu EQ 0) then $
    message, 'Error reading PSF file '+string(basisfile)

   phdr = ptrarr(nhdu)
   pbasis = ptrarr(nhdu)
   for ihdu=0, nhdu-1 do begin
      pbasis[ihdu] = ptr_new(mrdfits(basisfile, ihdu, hdr1))
      phdr[ihdu] = ptr_new(hdr1)
      if (ihdu EQ 0) then nfiber = sxpar(hdr1,'NAXIS2')
   endfor

   return
end
;------------------------------------------------------------------------------
pro bbspec_extract_trimpsf, pbasis, phdr, psffile, x0, x1, y0, y1, fib1, fib2, $
 ximg=ximg

   psftype = strtrim(sxpar(*phdr[0],'PSFTYPE'),2)
   ny = sxpar(*phdr[0],'NAXIS1')
   nfiber = sxpar(*phdr[0],'NAXIS2')

   for ihdu=0, n_elements(pbasis)-1 do begin
      basis1 = *pbasis[ihdu]
      hdr1 = *phdr[ihdu]

      ; Replace with the X centroids shifted, and trim to only the first entries
      ; if the PSF is only solved for the first fibers in the first rows
      if (ihdu EQ 0 AND keyword_set(ximg)) then basis1 = ximg[0:ny-1,0:nfiber-1]
 
     ; Replace the X and Y positions to refer to the subimage positions
      if (ihdu EQ 0) then basis1 -= x0
      if (ihdu EQ 1) then basis1 -= y0

      if (ihdu LE 2 OR psftype EQ 'GAUSS-HERMITE') then begin
         basis1 = basis1[y0:y1,fib1:fib2]
         sxaddpar, hdr1, 'NPIX_X', x1-x0+1
         sxaddpar, hdr1, 'NPIX_Y', y1-y0+1
         sxaddpar, hdr1, 'NSPEC', fib2-fib1+1
         sxaddpar, hdr1, 'NFLUX', y1-y0+1
         sxaddpar, hdr1, 'NAXIS', 2
         sxaddpar, hdr1, 'NAXIS1', y1-y0+1
         sxaddpar, hdr1, 'NAXIS2', fib2-fib1+1
      endif
      if (ihdu EQ 4 AND psftype EQ 'PCA-PIX') then begin
         basis1.x0 = basis1.x0 - x0
         basis1.y0 = basis1.y0 - y0
         basis1 = basis1[fib1:fib2]
      endif

      mwrfits, basis1, psffile, hdr1, create=(ihdu EQ 0), /silent
   endfor

   return
end
;------------------------------------------------------------------------------
pro bbspec_extract, image, invvar, flux, fluxivar, basisfile=basisfile, $
 ximg=ximg, frange=frange1, yrange=yrange1, tmproot=tmproot1, ymodel=ymodel

   if (n_params() NE 4 OR keyword_set(basisfile) EQ 0) then $
    message, 'Parameters not set'
   if (keyword_set(yrange1)) then yrange = yrange1 $
    else yrange = minmax(where(total(invvar,1) GT 0))
   if (keyword_set(tmproot1)) then tmproot = tmproot1 $
    else tmproot = 'tmp-'

   stime0 = systime(1)

   nfibper = 20 ; number of fibers to extract in each call
   npadx = 7 ; pad sub-image to left+right by this number of pixels
   nsmally = 100 ; number of rows to extract in each call
   npady = 15 ; number of rows to use as padding

   dims = size(image,/dimens)
   nx = dims[0]
   ny = dims[1]

   imgfile = tmproot+'img.fits'
   psffile = tmproot+'psf.fits'
   fluxfile = tmproot+'flux.fits'
   modfile = tmproot+'model.fits'

   ; Read the full PSF file
   bbspec_extract_readpsf, basisfile, pbasis, phdr, nfiber=nfiber

   if (keyword_set(frange1)) then begin
      if (frange1[0] LT 0 OR frange1[1] GE nfiber) then $
       message, 'FRANGE extends beyond the fiber numbers in PSF file'
      frange = frange1
   endif else begin
      frange = [0,nfiber-1]
   endelse

   flux = fltarr(ny,nfiber)
   fluxivar = fltarr(ny,nfiber)
   if (arg_present(ymodel)) then ymodel = fltarr(nx,ny)

   ; Loop through sub-images, solving for nsmall rows at a time
   ; on nfibper fibers
   nstepy = nsmally - 2*npady ; number of rows to step up in each call
   nchunk = ceil((yrange[1] - yrange[0] + 1 - 2*npady)/nstepy)
   for fib1=frange[0], frange[1], nfibper do begin
      fib2 = (fib1 + nfibper - 1) < frange[1]
      for ichunk=0, nchunk-1 do begin
         y0 = yrange[0] + ichunk * (nsmally - 2*npady)
         y1 = (y0 + nsmally - 1) < (ny-1)
         x0 = (fix(min(ximg[y0:y1,fib1:fib2])) - npadx) > 0
         x1 = (fix(max(ximg[y0:y1,fib1:fib2])) + npadx) < (nx-1)
         splog, 'Extracting FIBER=', fib1, fib2, ' Y=', y0, y1, ' '+systime()

         mwrfits, image[x0:x1,y0:y1], imgfile, /create, /silent
         mwrfits, invvar[x0:x1,y0:y1], imgfile, /silent

         bbspec_extract_trimpsf, pbasis, phdr, psffile, x0, x1, y0, y1, $
          fib1, fib2, ximg=ximg

         pyfile = djs_filepath('pix2spec.py', root_dir=getenv('BBSPEC_DIR'), $
          subdir='examples')
         cmd = 'python '+pyfile+' -i '+imgfile+' -p '+psffile+' -o '+fluxfile
         spawn, cmd, res, errcode
; Should not need to ignore some error messages below???
         if (keyword_set(errcode) $
          AND strmatch(errcode[0],'*LinAlgError*') EQ 0 $
          AND strmatch(errcode[0],'*Singular matrix*') EQ 0) then begin
            splog, errcode
            message, 'Error calling '+cmd
         endif
         flux1 = mrdfits(fluxfile, /silent)
         fluxivar1 = mrdfits(fluxfile, 1, /silent)

         pyfile = djs_filepath('spec2pix.py', root_dir=getenv('BBSPEC_DIR'), $
          subdir='examples')
         cmd = 'python '+pyfile+' -i '+fluxfile+' -p '+psffile+' -o '+modfile $
          + ' --hdu 3'
         spawn, cmd, res, errcode
         if (keyword_set(errcode)) then begin
            splog, errcode
            message, 'Error calling '+cmd
         endif
         if (arg_present(ymodel)) then ymodel1 = mrdfits(modfile, /silent)
; The test for NaNs shouldn't be necessary!???
         ibad = where(finite(flux1) EQ 0 OR finite(fluxivar1) EQ 0, nbad)
         if (nbad GT 0) then begin
            splog, 'Replacing ', nbad, ' NaN values in extracted flux'
            flux1[ibad] = 0
            fluxivar1[ibad] = 0
         endif
         if (ichunk EQ 0) then trim1 = 0 else trim1 = npady
         if (ichunk EQ nchunk-1) then trim2 = 0 else trim2 = npady
         flux[y0+trim1:y1-trim2,fib1:fib2] = flux1[trim1:y1-y0-trim2,*]
         fluxivar[y0+trim1:y1-trim2,fib1:fib2] = fluxivar1[trim1:y1-y0-trim2,*]
         ; Use the model image +/- 3.5 pix from the central fiber...
         if (arg_present(ymodel)) then $
          ymodel[x0:x1,y0+trim2:y1-trim2] = ymodel1[*,trim2:y1-y0-trim2]
      endfor
   endfor

   for i=0, n_elements(phdr)-1 do ptr_free, phdr[i]
   for i=0, n_elements(pbasis)-1 do ptr_free, pbasis[i]

   splog, 'Time to bbspec = ', systime(1)-stime0, ' seconds'

   return
end
;------------------------------------------------------------------------------
