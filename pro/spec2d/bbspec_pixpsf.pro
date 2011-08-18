;------------------------------------------------------------------------------
; Generate pixelated PSF
;------------------------------------------------------------------------------
pro bbspec_pixpsf, arcstr, flatstr, pradius=pradius, rradius=rradius, $
 npoly=npoly, outfile=outfile1, batch=batch

;arcstr='r1-00115982' ; ???
;flatstr='r1-00115981' ; ???

   if (NOT keyword_set(npoly)) then npoly = [2,3]
   if (NOT keyword_set(pradius)) then pradius = 8. + fltarr(npoly[0]*npoly[1])
   if (NOT keyword_set(rradius)) then rradius = 0. + fltarr(npoly[0]*npoly[1])
   if (n_elements(pradius) NE n_elements(rradius)) then $
    message, 'Number of elements for PRADIUS,RRADIUS must agree'
   if (keyword_set(outfile1)) then outfile = outfile1 $
    else outfile ='spBasisPSF-PIX-'+arcstr+'.fits'

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

   rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   mjdstr = string(sxpar(archdr, 'MJD'),format='(i5.5)')
   indir = concat_dir(rawdata_dir, mjdstr)
   arcname = 'sdR-'+arcstr+'.fit'
   sdssproc, arcname, image, ivar, indir=indir, $
    /applybias, /applypixflat, /applycrosstalk
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

   objs_all = reform(psolve_obj_struct(n_elements(xpix)),nlamp,nfiber)
   objs_all.xcen = xpix
   objs_all.ycen = ypix
;   objs_all.flux = rebin(lamps.intensity,nlamp,nfiber)
   objs_all[*].flux = djs_phot(objs_all[*].xcen, objs_all[*].ycen, 3., 0, image, $
    calg='none', salg='none', /quick)

   objs_all.goodmask = 1B
   fibernum = djs_laxisgen([nlamp,nfiber], iaxis=1)

;atv,image*(ivar ne 0)
;jj=where(objs.bestmask)
;atvplot,xpix[jj],ypix[jj],ps=1,syms=0.5,color='green'

   ;----------
   ; Solve for PSF, fixing which stars to use (nselect=0)
   ; and fixing the centers (maxshift=0)

   psfpix = psolve_pixelization(pradius=pradius, rradius=rradius)
   xpad = 0
   ypad = 0
   niter = 3
   maxshift = 0
   fixpsf = 0

   ngroup = nfiber/20
   fakeimg = 0
   filename = 'tmppsf-'+arcstr+'-'+string(lindgen(ngroup),format='(i2.2)')+'.ss'
   pbsfile = 'pbs-tmppsf-'+arcstr+'-'+string(lindgen(ngroup),format='(i2.2)')
   jobid = lonarr(ngroup)
   for igroup=0, ngroup-1 do begin
      splog, 'Generating PSF for group ', igroup, ngroup
      objs = objs_all

      ; Choose every fiber with >1000 counts for PSF construction
      objs.bestmask = rebin(strmatch(lamps.use_wset,'*GOOD*'),nlamp,nfiber) $
       AND (fibernum GE igroup*20 AND fibernum LT (igroup+1)*20) $
       AND (objs.flux GT 1000)
      ; Additional trimming for edge effects...???
      ; Include stars in the neighboring bundles to left + right for blending purposes.
      itrim = where(fibernum GE igroup*20-20 AND fibernum LT (igroup+1)*20+20 $
       AND objs.xcen GE 10 AND objs.xcen LE nx-10 $
       AND objs.ycen GE 10 AND objs.ycen LE ny-10)
      objs = objs[itrim]

      skyimg = 0 ; re-fit the sky image on each group of fibers

      if (NOT keyword_set(batch)) then begin
         psolve_iter, image, ivar, objs, psfpix, psfimg, skyimg, $
          xpad=xpad, ypad=ypad, npoly=npoly, niter=niter, maxshift=maxshift, fixpsf=fixpsf
      endif else begin
         fq = "'"
         save, filename=filename[igroup], image, ivar, objs, psfpix, skyimg, $
          xpad, ypad, npoly, niter, maxshift, fixpsf
         jobid[igroup] = bbspec_batch('psolve_iter,filename='+fq+filename[igroup]+fq, $
          pbsfile=pbsfile[igroup])
      endelse

      if (NOT keyword_set(batch)) then begin
         if (NOT keyword_set(psfimg_all)) then $
          psfimg_all = fltarr([size(psfimg,/dimens),ngroup])
         psfimg_all[*,*,*,igroup] = psfimg
         fakeimg += skyimg
         psolve_addstars, fakeimg, psfimg, objs
      endif
   endfor

   if (keyword_set(batch)) then begin
       bbspec_batch_wait, jobid
      for igroup=0, ngroup-1 do begin
         restore, filename=filename[igroup]
         if (NOT keyword_set(psfimg_all)) then begin
            psfimg_all = fltarr([size(psfimg,/dimens),ngroup])
            fakeimg += skyimg
         endif
         psfimg_all[*,*,*,igroup] = psfimg
         indx = where(objs.bestmask, ct)
         if (ct GT 0) then $
          psolve_addstars, fakeimg, psfimg, objs
      endfor
   endif

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
   sxdelpar, outhdr0, 'SIMPLE'  ; Can't have SIMPLE in extensions
   mwrfits, ally, outfile, outhdr0
   sxaddpar, outhdr0, 'PSFPARAM', 'LogLam'
   mwrfits, loglam, outfile, outhdr0

   ; HDU #3 has list of polynomial exponents
   polydat = replicate(create_struct('IMODEL', 0L, 'XEXP', 0L, 'YEXP', 0L), $
    npoly[0]*npoly[1])
   polydat.imodel = lindgen(npoly[0]*npoly[1])
   polydat.xexp = (djs_laxisgen(npoly, iaxis=0))[*]
   polydat.yexp = (djs_laxisgen(npoly, iaxis=1))[*]
   mwrfits, polydat, outfile

   ; HDU #4 has the indexing for each fiber ID
   fibdat = replicate(create_struct('IGROUP', 0L, 'X0', 0.0, 'XSCALE', 0.0, $
    'Y0', 0.0, 'YSCALE', 0.0), nfiber)
   fibdat.igroup = lindgen(nfiber) / 20L
   fibdat.x0 = 0
   fibdat.xscale = 0.001
   fibdat.y0 = 0
   fibdat.yscale = 0.001
   mwrfits, fibdat, outfile

   ; HDU #5 has the PSF images indexed [X,Y,IMODEL,IGROUP]
   mkhdr, outhdr1, psfimg_all
   mwrfits, psfimg_all, outfile

   return
end
;------------------------------------------------------------------------------
