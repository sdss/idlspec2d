
;------------------------------------------------------------------------------
pro fluxrebin, filename, newloglam, newflux, newivar, adderr=adderr

; Should we actually spline fit for the re-binning???

   ;----------
   ; Case where file does not exist

   if (NOT keyword_set(filename)) then begin
      newflux = 0
      newivar = 0
      return
   endif

   ;----------
   ; Read in the flux, errors, mask, wavelengths

   objflux = mrdfits(filename,0)
   objivar = mrdfits(filename,1)
   objmask = mrdfits(filename,2)
   wset = mrdfits(filename,3)
   traceset2xy, wset, xx, objloglam

   dims = size(objloglam, /dimens)
   npix = dims[0]
   nfiber = dims[1]

   ;----------
   ; Add an additional error term equal to ADDERR of the flux.

   if (keyword_set(adderr)) then begin
      gmask = objivar NE 0 ; =1 for good points
      objivar = 1.0 / ( 1.0/(objivar + (1-gmask)) $
       + (adderr * (objflux>0))^2 ) * gmask
   endif

   ;----------
   ; Mask where the spectrum may be dominated by sky-subtraction
   ; residuals.  Grow that mask by 2 pixels in each direction.

   smask = (objmask AND pixelmask_bits('BRIGHTSKY')) NE 0
   for ifiber=0L, nfiber-1 do $
    smask[*,ifiber] = smooth(float(smask[*,ifiber]),5) GT 0
   ibad = where(smask)
   if (ibad[0] NE -1) then objivar[ibad] = 0

   ;----------
   ; Make a map of the size of each pixel in delta-(log10-Angstroms),
   ; and re-normalize the flux to ADU/(dloglam)

   dlogimg = [ objloglam[1,*] - objloglam[0,*], $
    0.5 * (objloglam[2:npix-1,*] - objloglam[0:npix-3,*]), $
    objloglam[npix-1,*] - objloglam[npix-2,*] ]
   dlogimg = abs(dlogimg)

   dloglam = newloglam[1] - newloglam[0]
   divideflat, objflux, invvar=objivar, (dlogimg/dloglam), minval=0

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
pro myfluxcorr, bsmearfile, rsmearfile, bsciname, rsciname, corrfile, $
 adderr=adderr

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

   fluxrebin, bsmearfile, newloglam, bsmearflux, bsmearivar, adderr=adderr
   fluxrebin, rsmearfile, newloglam, rsmearflux, rsmearivar, adderr=adderr

   fluxrebin, bsciname, newloglam, bsciflux, bsciivar, adderr=adderr
   fluxrebin, rsciname, newloglam, rsciflux, rsciivar, adderr=adderr

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
        ssmearflux[*,ifiber], ssmearivar[*,ifiber], $
        ssciflux[*,ifiber], ssciivar[*,ifiber], $
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
   ; Require the median S/N per pix for the science exposure to be > 2.0.
   ; Require the median S/N per pix for the smear exposure to be > 2.0.
   ; Require the fits for a fiber are never less than 0.2 * median.
   ; Require the fits for a fiber are never greater than 5 * median.
   ; Require the fiber isn't SKY.

   smearsnmed = djs_median(ssmearflux * sqrt(ssmearivar), 1)

   scisnmed = djs_median(ssciflux * sqrt(ssciivar), 1)

   medfit = median(fitimg)
   splog, 'Median flux-correction factor = ', medfit

   qgood = scisnmed GT 2.0 $
    AND smearsnmed GT 2.0 $
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

;stop ; ???
;i=0
;splot,ssmearflux[*,i]
;soplot,ssciflux[*,i],color='red'
;soplot,ssciflux[*,i]*fitimg[*,i],color='green'

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
