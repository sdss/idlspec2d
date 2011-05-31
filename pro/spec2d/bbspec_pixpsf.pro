;------------------------------------------------------------------------------
; Generate pixelated PSF
;------------------------------------------------------------------------------
pro bbspec_pixpsf, arcstr, flatstr, pradius=pradius, rradius=rradius, $
 npoly=npoly, outfile=outfile1

arcstr='r1-00115982' ; ???
flatstr='r1-00115981' ; ???

   if (NOT keyword_set(npoly)) then npoly = [2,4]
   if (NOT keyword_set(pradius)) then pradius = 8. + fltarr(npoly[0]*npoly[1])
   if (NOT keyword_set(rradius)) then rradius = 0. + fltarr(npoly[0]*npoly[1])
   if (n_elements(pradius) NE n_elements(rradius)) then $
    message, 'Number of elements for PRADIUS,RRADIUS must agree'
   if (keyword_set(outfile1)) then outfile = outfile1 $
    else outfile ='spBasisPSF-'+arcstr+'.fits'

   ;----------
   ; Read the pipeline arc + flat files

   arcfile = 'spArc-'+arcstr+'.fits*'
   arcfile = (findfile(arcfile))[0]
   flatfile = 'spFlat-'+flatstr+'.fits*'
   flatfile = (findfile(flatfile))[0]
   archdr = headfits(arcfile)
   xset = mrdfits(flatfile,1)
   wset = mrdfits(arcfile,2)

   ;----------
   ; Read the raw arc image

   rawdata_dir = getenv('RAWDATA_DIR')
   mjdstr = string(sxpar(archdr, 'MJD'),format='(i5.5)')
   indir = concat_dir(rawdata_dir, mjdstr)
   arcname = 'sdR-'+arcstr+'.fit'
   sdssproc, arcname, image, invvar, indir=indir, $
    /applybias, /applypixflat, /applycrosstalk
;image = mrdfits('fakeimg.fits') ; ???
   if (NOT keyword_set(image)) then $
    message, 'Error reading file '+arcname
   dims = size(image, /dimens)
   nx = dims[0]
   ny = dims[1]

   ;----------
   ; Read the line list, and keep blends
   ; Note we are working in *air* not vacuum wavelengths here

   lampfile = filepath('lamplines.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='opfiles')
   splog, 'Reading lamp file ', lampfilename
   lamps = yanny_readone(lampfile)
   if (NOT keyword_set(lamps)) then $
    message, 'LAMPFILE not valid '+lampfile

   ;----------
   ; Determine the expected locations of arc lines

   ypix = traceset2pix(wset, alog10(lamps.lambda))
   traceset2xy, xset, ypix, xpix
   dims = size(ypix,/dimens)
   nlamp = dims[0]
   nfiber = dims[1]

   objs = reform(psolve_obj_struct(n_elements(xpix)),nlamp,nfiber)
   objs.xcen = xpix
   objs.ycen = ypix
;   objs.flux = rebin(lamps.intensity,nlamp,nfiber)
   objs[*].flux = djs_phot(objs[*].xcen, objs[*].ycen, 3., 0, image, $
    calg='none', salg='none', /quick)

   objs.goodmask = 1B
   ; Choose every 20th fiber with >3000 counts for PSF construction
   fibernum = djs_laxisgen([nlamp,nfiber], iaxis=1)
;   objs.bestmask = rebin(strmatch(lamps.use_psf,'*YES_PSF*'),nlamp,nfiber)
;   objs.bestmask = rebin(strmatch(lamps.use_wset,'*GOOD*'),nlamp,nfiber) $
;    AND (fibernum MOD 20) EQ 0 AND objs.flux GT 3000
   objs.bestmask = rebin(strmatch(lamps.use_wset,'*GOOD*'),nlamp,nfiber) $
    AND (fibernum LT 25) AND (objs.flux GT 3000)
objs = objs[where(fibernum LT 25)]
;invvar[415:*,*] = 0 ; zero-out everything after 20 fibers ???
;objs.bestmask = rebin(strmatch(lamps.use_psf,'*YES_PSF*'),nlamp,nfiber) $
; AND (fibernum MOD 5) EQ 0 ; every 5th fiber on these lines

   ; Trim to only objects on the image
   itrim = where(objs.xcen GE 0 AND objs.xcen LE nx-1 $
    AND objs.ycen GE 0 AND objs.ycen LE ny-1)
   objs = objs[itrim]
; Additional trimming for edge effects...
itrim = where(objs.xcen GE 10 AND objs.xcen LE nx-10 $
 AND objs.ycen GE 10 AND objs.ycen LE ny-10)
objs = objs[itrim]
;atv,image*(invvar ne 0)
;jj=where(objs.bestmask)
;atvplot,xpix[jj],ypix[jj],ps=1,syms=0.5,color='green'

   ;----------
   ; Solve for PSF, fixing which stars to use (nselect=0)
   ; and fixing the centers (maxshift=0)

   psfpix = psolve_pixelization(pradius=pradius, rradius=rradius)
;   skyimg = 0.*image
   psolve_iter, image, invvar, objs, psfpix, psfimg, skyimg, $
    xpad=0, ypad=0, npoly=npoly, niter=3, maxshift=0, fixpsf=0 ; shift???

   fakeimg = skyimg
   psolve_addstars, fakeimg, psfimg, objs

;stop
;jj=where(objs.bestmask)
;atv,image-fakeimg
;atvplot,objs.xcen,objs.ycen,ps=1,syms=0.5,color='red'
;atvplot,objs[jj].xcen,objs[jj].ycen,ps=4,syms=1,color='red'

   ;----------
   ; Write the output PSF file

   splog, 'Writing PSF file '+outfile
   traceset2xy, wset, yy, loglam
   traceset2xy, xset, ally, allx
   mkhdr, outhdr0, allx
   sxaddpar, outhdr0, 'PSFTYPE', 'PCA-PIX'
   sxaddpar, outhdr0, 'NPIX_X', nx
   sxaddpar, outhdr0, 'NPIX_Y', ny
   sxaddpar, outhdr0, 'NFLUX', ny
   sxaddpar, outhdr0, 'NSPEC', nfiber
   sxaddpar, outhdr0, 'PSFPARAM', 'X'
   mwrfits, allx, outfile, outhdr0, /create
   sxaddpar, outhdr0, 'PSFPARAM', 'Y'
   mwrfits, ally, outfile, outhdr0
   sxaddpar, outhdr0, 'PSFPARAM', 'LogLam'
   mwrfits, loglam, outfile, outhdr0
   mkhdr, outhdr1, psfimg[*,*,0]
   for iy=0, npoly[1]-1 do begin
      for ix=0, npoly[0]-1 do begin
         if (ix EQ 0 AND iy EQ 0) then begin
            psfparam = 'const'
         endif else begin
            psfparam = ''
            for j=1, ix do psfparam += 'x'
            for j=1, iy do psfparam += 'y'
         endelse
         sxaddpar, outhdr1, 'PSFPARAM', psfparam
         mwrfits, psfimg[*,*,ix+iy*npoly[0]], outfile, outhdr1
      endfor
   endfor

stop
   return
end
;------------------------------------------------------------------------------
