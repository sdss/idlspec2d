
;------------------------------------------------------------------------------
; Compute the polynomial function to multiply by YDATA1 to get YDATA2.
; ??? Need to iterate with the values of COEFF, and applying those to
;     the errors in YDATA1.  Also, should iterate each fit and reject points.
function myratiofit, xdata, ydata1, yivar1, ydata2, yivar2, ncoeff=ncoeff

   ndata = n_elements(xdata)
   acoeff = dblarr(ncoeff)
   acoeff[0] = 1.0

   ; Renormalize X to [-1,1]
   xnorm = 2.0 * (xdata - min(xdata)) / (max(xdata) - min(xdata)) -1.0

   ; Construct the matrix A, such that the first row is YDATA1,
   ; the second row is YDATA1 * XNORM, the third is YDATA1 * XNORM^2
   amatrix = fltarr(ndata,ncoeff)
   for icoeff=0, ncoeff-1 do $
    amatrix[*,icoeff] = ydata1 * xnorm^icoeff

   invsig = sqrt(yivar1 + yivar2)

   mmatrix = amatrix
   for icoeff=0, ncoeff-1 do $
    mmatrix[*,icoeff] = mmatrix[*,icoeff] * invsig
   bvec = ydata2 * invsig

   mmatrixt = transpose( mmatrix )
   mm = mmatrixt # mmatrix

   ; Use SVD to invert the matrix
;   mmi = invert(mm, /double)
   if (ncoeff EQ 1) then begin
      mmi = 1.0 / mm
   endif else begin
      svdc, mm, ww, uu, vv, /double
      mmi = 0 * vv
      for i=0L, ncoeff-1 do mmi[i,*] = vv[i,*] / ww[i]
      mmi = mmi ## transpose(uu)
   endelse

   acoeff = mmi # (mmatrixt # bvec)
   yfit = fltarr(ndata) + acoeff[0]
   for icoeff=1, ncoeff-1 do $
    yfit = yfit + acoeff[icoeff] * xnorm^icoeff

   return, yfit
end

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

   ;----------
   ; Read the plug-map file (for identifying sky fibers)

   plugmap = mrdfits(bsmearfile,5)

   ;----------
   ; Re-bin all the spectra to the same, uniform log-lambda scale.
   ; Consistently mask exactly the same pixels in the smear images
   ; as in the science images.

   minlog = alog10(3700)
   maxlog = alog10(9300)
   dloglam = 1.0d-4
   nnewpix = fix( (maxlog - minlog) / dloglam )
   newloglam = minlog + findgen(nnewpix) * dloglam

   fluxrebin, bcalibset, bsmearfile, newloglam, bsmearflux, bsmearivar
   fluxrebin, rcalibset, rsmearfile, newloglam, rsmearflux, rsmearivar

   fluxrebin, bcalibset, bsciname, newloglam, bsciflux, bsciivar
   fluxrebin, rcalibset, rsciname, newloglam, rsciflux, rsciivar

   ;----------
   ; Mask the lowest S/N points in both smear + science

;   flux_mask_low_sn, bsmearflux, bsmearivar, bsciflux, bsciivar
;   flux_mask_low_sn, rsmearflux, rsmearivar, rsciflux, rsciivar

   ;----------
   ; Add the blue + red flux images

   ssmearflux = bsmearflux * (bsmearivar NE 0) + rsmearflux * (rsmearivar NE 0)
   ssmearivar = bsmearivar + rsmearivar
   ssciflux = bsciflux * (bsciivar NE 0) + rsciflux * (rsciivar NE 0)
   ssciivar = bsciivar + rsciivar

   ;----------

ncoeff = 3 ; ???
   nfiber = (size(bsmearflux,/dimens))[1]
   fitimg = dblarr(nnewpix,nfiber)
   for ifiber=0, nfiber-1 do begin
      fitimg[*,ifiber] = $
       myratiofit(newloglam, $
        ssciflux[*,ifiber], ssciivar[*,ifiber], $
        ssmearflux[*,ifiber], ssmearivar[*,ifiber], $
        ncoeff=ncoeff)
   endfor

   ;----------
   ; Special case: If the science image and smear image is the same,
   ; then force their ratio to be unity.

   if (bsmearfile EQ bsciname AND rsmearfile EQ rsciname) then begin
      fitimg[*] = 1 ; Should already be true
   endif

   ;----------
   ; Convert these fits into a trace set

   xy2traceset, newloglam # (bytarr(nfiber)+1), fitimg, corrset, ncoeff=ncoeff

   ;----------
   ; Which fits are considered bad?
   ; Require the median S/N per pix for the science exposure to be > 1.0.
   ; Require the median S/N per pix for the smear exposure to be > 0.5.
   ; Require the fits for a fiber are never less than 0.2 * median.
   ; Require the fits for a fiber are never greater than 5 * median.
   ; Require the fiber isn't SKY.

   ssmearsnmed = djs_median(ssmearflux * sqrt(ssmearivar), 1)

   scisnmed = djs_median(ssciflux * sqrt(ssciivar), 1)

   medfit = median(fitimg)
   splog, 'Median flux-correction factor = ', medfit

   qgood = scisnmed GT 1.0 $
    AND ssmearsnmed GT 0.5 $
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
   xy2traceset, newloglam, meanfit, cset1, ncoeff=ncoeff

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
