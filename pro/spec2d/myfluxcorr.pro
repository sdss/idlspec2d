
;------------------------------------------------------------------------------
; Mask the same low S/N points in both a smear image and a science image.
; Low S/N points are those with a S/N less than 10% of the median S/N.
; Mask the points be setting the inverse variance to zero.

pro flux_mask_low_sn, smearflux, smearivar, sciflux, sciivar

   nfiber = (size(smearflux,/dimens))[1]

   for ifiber=0, nfiber-1 do begin
      snvec1 = smearflux[*,ifiber] * sqrt(smearivar[*,ifiber])
      snvec2 = sciflux[*,ifiber] * sqrt(sciivar[*,ifiber])
      indx1 = where(snvec1 NE 0, ct1)
      indx2 = where(snvec2 NE 0, ct2)
      if (ct1 LT 100 OR ct2 LT 100) then begin
         smearivar[*,ifiber] = 0
         sciivar[*,ifiber] = 0
      endif else begin
         snmed1 = median( snvec1[indx1] )
         snmed2 = median( snvec2[indx2] )
         ibad = where(snvec1 LE 0.10 * snmed1 OR snvec2 LE 0.1 * snmed2)
         if (ibad[0] NE -1) then begin
            smearivar[ibad,ifiber] = 0
            sciivar[ibad,ifiber] = 0
         endif
      endelse
   endfor

   return
end

;------------------------------------------------------------------------------
pro fluxrebin, calibset, filename, newloglam, newflux, newivar

   ;----------
   ; Read in the flux, errors, mask, wavelengths

;   objflux = mrdfits(filename,0)
objflux = mrdfits(filename,0) * mrdfits(filename,6) ; ???
;   objivar = mrdfits(filename,1)
objivar = mrdfits(filename,1) * mrdfits(filename,6) ; ???
   objmask = mrdfits(filename,2)
   wset = mrdfits(filename,3)
   traceset2xy, wset, junk, objloglam

   dims = size(objloglam, /dimens)
   npix = dims[0]
   nfiber = dims[1]

   ;----------
   ; Mask where the spectrum may be dominated by sky-subtraction
   ; residuals.  Grow that mask by 2 pixels in each direction.

   skymask = (objmask AND pixelmask_bits('BRIGHTSKY')) NE 0
   for ifiber=0L, nfiber-1 do $
    skymask[*,ifiber] = smooth(float(skymask[*,ifiber]),5) GT 0
   ibad = where(skymask)
   if (ibad[0] NE -1) then objivar[ibad] = 0

   ;----------
   ; Apply the flux-calibration vector
   ; In this case, don't trust points with a low calibration factor.

; ???
;   calibfac = bspline_valu(objloglam, calibset)
;   divideflat, objflux, objivar, calibfac, minval=0.20*mean(calibfac)

   ;----------
   ; Make a map of the size of each pixel in delta-(log10-Angstroms),
   ; and re-normalize the flux to ADU/(dloglam)

   dlogimg = [ objloglam[1,*] - objloglam[0,*], $
    0.5 * (objloglam[2:npix-1,*] - objloglam[0:npix-3,*]), $
    objloglam[npix-1,*] - objloglam[npix-2,*] ]
   dlogimg = abs(dlogimg)

   dloglam = newloglam[1] - newloglam[0]
   divideflat, objflux, objivar, (dlogimg/dloglam), minval=0

   ;----------
   ; Linearly interpolate the data onto the new wavelengths

   nnewpix = n_elements(newloglam)
   newflux = fltarr(nnewpix, nfiber)
   newivar = fltarr(nnewpix, nfiber)
   newmask = bytarr(nnewpix, nfiber)

   for ifiber=0, nfiber-1 do begin
      w0 = 3 * objloglam[0,ifiber] - 2 * objloglam[1,ifiber]
      w1 = 2 * objloglam[0,ifiber] - 1 * objloglam[1,ifiber]
      w2 = 2 * objloglam[npix-1,ifiber] - 1 * objloglam[npix-2,ifiber]
      w3 = 3 * objloglam[npix-1,ifiber] - 2 * objloglam[npix-2,ifiber]
      newflux[*,ifiber] = interpol([0,0,objflux[*,ifiber],0,0], $
       [w0,w1,objloglam[*,ifiber],w2,w3], newloglam)
      ; Beware that round-off could put the inverse variance less than zero.
      newivar[*,ifiber] = interpol([0,0,objivar[*,ifiber],0,0], $
       [w0,w1,objloglam[*,ifiber],w2,w3], newloglam) > 0
      newmask[*,ifiber] = floor( $
       interpol([0,0,float(objivar[*,ifiber] NE 0),0,0], $
       [w0,w1,objloglam[*,ifiber],w2,w3], newloglam) )
   endfor

   newivar = newivar * newmask

   return
end

;------------------------------------------------------------------------------
pro myfluxcorr, bsmearfile, rsmearfile, bcalibset, rcalibset, $ 
 bsciname, rsciname, corrfile

   binsz = 100 ; in pixels for rebinned images

   ;----------
   ; Read the plug-map file (for identifying sky fibers)

   plugmap = mrdfits(bsmearfile,5)

   ;----------
   ; Re-bin all the spectra to the same, uniform log-lambda scale.
   ; Consistently mask exactly the same pixels in the smear images
   ; as in the science images.
   ; The new binning will have an integral number of BINSZ pixels

   minlog = alog10(3700)
   maxlog = alog10(9300)
   dloglam = 1.0d-4
   nnewpix = fix( (maxlog - minlog) / (dloglam * binsz) ) * binsz
   newloglam = minlog + findgen(nnewpix) * dloglam

   fluxrebin, bcalibset, bsmearfile, newloglam, bsmearflux, bsmearivar
   fluxrebin, rcalibset, rsmearfile, newloglam, rsmearflux, rsmearivar

   fluxrebin, bcalibset, bsciname, newloglam, bsciflux, bsciivar
   fluxrebin, rcalibset, rsciname, newloglam, rsciflux, rsciivar

   ;----------
   ; Mask the lowest S/N points in both smear + science

   flux_mask_low_sn, bsmearflux, bsmearivar, bsciflux, bsciivar
   flux_mask_low_sn, rsmearflux, rsmearivar, rsciflux, rsciivar

   ;----------
   ; Add the blue + red flux images

   ssmearflux = bsmearflux * (bsmearivar NE 0) + rsmearflux * (rsmearivar NE 0)
   ssmearivar = bsmearivar + rsmearivar
   ssciflux = bsciflux * (bsciivar NE 0) + rsciflux * (rsciivar NE 0)
   ssciivar = bsciivar + rsciivar

   ;----------
   ; Rebin these images into only NBIN wavelength bins.

   nbin = nnewpix / binsz
   nfiber = (size(bsmearflux,/dimens))[1]

   binloglam = rebin(newloglam, nbin)
   binlogarr = binloglam # (bytarr(nfiber)+1)
   binsmearflux = rebin(ssmearflux, nbin, nfiber)
   binsmearivar = rebin(ssmearivar, nbin, nfiber)
   binsciflux = rebin(ssciflux, nbin, nfiber)
   binsciivar = rebin(ssciivar, nbin, nfiber)

   ;----------
   ; Take the ratio smear/science.
   ; Construct a mask of the high S/N bins in the binned science exposure.
   ; Assume that we can deal with low S/N bins in the smear exposure,
   ; since we're taking the ratio smear/science.

   minflux = 1.0 ; in ADU/pix
   binratio = (binsmearflux > minflux) / (binsciflux > minflux)
   binmask = (binsciflux * sqrt(binsciivar) GT 20) $
    AND binsmearflux GT minflux $
    AND binsciflux GT minflux

   ;----------
   ; Special case: If the science image and smear image is the same,
   ; then force their ratio to be unity.

   if (bsmearfile EQ bsciname AND rsmearfile EQ rsciname) then begin
      binratio[*] = 1 ; Should already be true
      binmask[*] = 1
   endif

   ;----------
   ; Fit the ratio image on a fiber-by-fiber basis.
   ; Reject points whose fit values differ by more than 50% of the median
   ;   BINRATIO value.  For ex, if the typical value of BINRATIO is 0.2,
   ;   then reject any points that outlie by more than 0.1 from the fit.

ncoeff = 3 ; ???
   xy2traceset, binlogarr, binratio, corrset, ncoeff=ncoeff, $
    maxdev=0.50*median(binratio), inmask=binmask, outmask=outbinmask

   ;----------
   ; Which fits are considered bad?
   ; Require the median S/N per pix for the smear exposure to be > 0.1.
   ; Require the median S/N per pix for the science exposure to be > 1.0.
   ; Require the fit to have retained at least 50% of the binned points.
   ; Require the fits for a fiber are never less than 0.2 * median.
   ; Require the fits for a fiber are never greater than 5 * median.
   ; Require the fiber isn't SKY.

   ssmearmask = ssmearivar NE 0
   ssmearsnmed = djs_median(ssmearflux * sqrt(ssmearivar>0) * ssmearmask, 1)

   sscimask = ssciivar NE 0
   scisnmed = djs_median(ssciflux * sqrt(ssciivar>0) * sscimask, 1)

   traceset2xy, corrset, binlogarr, fitimg
   medfit = median(fitimg)
   splog, 'Median flux-correction factor = ', medfit

   qgood = scisnmed GT 1.0 $
    AND ssmearsnmed GT 0.1 $
    AND total(outbinmask,1) GT 0.50*nbin $
    AND total(fitimg LT 0.20 * medfit, 1) EQ 0 $
    AND total(fitimg GT 5.0 * medfit, 1) EQ 0 $
    AND strtrim(plugmap.objtype,2) NE 'SKY'

   igood = where(qgood EQ 1, ngood)
   splog, 'Total of', ngood, ' good flux-correction vectors'
   if (igood[0] EQ -1) then $
    message, 'Trouble in flux-correction fits'

   ;----------
   ; Make a mean fit for the good fits
 
   meanfit = total(fitimg[*,igood],2) / ngood
   xy2traceset, binloglam, meanfit, cset1, ncoeff=ncoeff

   ;----------
   ; Replace bad fits with the mean fit

   ibad = where(qgood EQ 0, nbad)
   splog, 'Replacing', nbad, ' bad flux-correction vectors'
   if (ibad[0] NE -1) then begin
      for i=0, nbad-1 do $
       corrset.coeff[*,ibad[i]] = cset1.coeff
   endif

   mwrfits, corrset, corrfile, /create

;stop ; ???
;   traceset2xy, corrset, binlogarr, fitimg2
;   plothist, fitimg2, bin=0.01

   return
end
