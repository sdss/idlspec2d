;------------------------------------------------------------------------------
;+
; NAME:
;   fitfineconst
;
; PURPOSE:
;   Steinhardt's project to fit the fine structure constant from O_III lines
;
; CALLING SEQUENCE:
;   fitfineconst, plate, fiber, mjd=, [ /doplot, $
;    bestshift=, besterr= ]
;
; INPUTS:
;   plate      - Plate number
;   fiber      - Fiber number
;
; REQUIRED KEYWORDS:
;   mjd        - MJD number for this PLATE
;
; OPTIONAL INPUTS:
;   doplot     - If set, then output a bunch of plots
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   bestshift  - Best shift between lines in log10(Ang)
;   besterr    - Formal error in BESTSHIFT
;
; COMMENTS:
;
; EXAMPLES:
;   IDL> fitfineconst, 276, 251, mjd=51909, /doplot, $
;        bestshift=bestshift, besterr=besterr
;
; PROCEDURES CALLED:
;   combine1fiber
;   djs_oplot
;   find_nminima
;   readspec
;
; REVISION HISTORY:
;   09-Jul-2002  Written by D. Schlegel & C. Steinhardt, Princeton.
;------------------------------------------------------------------------------
pro fitfineconst, plate, fiber, mjd=mjd, doplot=doplot, $
 bestshift=bestshift, besterr=besterr

   if (n_elements(plate) NE 1 OR n_elements(fiber) NE 1 $
    OR n_elements(mjd) NE 1) then begin
      doc_library, 'fitfineconst'
      return
   endif

   plottitle = string(plate, fiber, mjd, $
    format='("Plate ",i4," Fiber ",i3," MJD ", i5)')
   pixscale = 1.d-4 ; pixel scale in log10(lambda)
   hwidth = 10 * pixscale ; 10 pixels half-width for fitting domain
   npoly = 2 ; Number of polynomial terms for continuum-fitting
   maxshift = 2.0 ; Maximum pixel shift to search
   nsamp = 20 ; Sub-sampling factor for pixels
   wtime = 0.05 ; Pause time in seconds between plots

   ;----------
   ; Read the spectrum and its redshift

   readspec, plate, fiber, mjd=mjd, $
    loglam=loglam, flux=objflux, invvar=objivar, zans=zans

   ; Compute the rest-frame wavelengths...
   rloglam = loglam - alog10(1. + zans.z)

   ;----------
   ; Set the wavelengths to fit -- the O_III lines at 4959,5007

   wave1 = 4958.911d0
   wave2 = 5006.843d0
   airtovac, wave1
   airtovac, wave2
   llam1 = alog10(wave1)
   llam2 = alog10(wave2)
   linesep = (llam2 - llam1)

   ;----------
   ; Oversample the object spectrum + its errors by a factor of NSAMP

   bigloglam = loglam[0] + (pixscale/nsamp) $
    * lindgen((n_elements(loglam)+1)*nsamp)
   combine1fiber, loglam, objflux, objivar, $
    newloglam=bigloglam, newflux=bigflux, newivar=bigivar, maxiter=0

   ;----------
   ; Loop over redshifts

   ; Pixel shifts from [-4,+4] pixels spaced every (1/nsamp) pixels
   nshift = 2 * long(maxshift * nsamp) + 1
   bigpixshift = (lindgen(nshift)-(nshift-1)/2)

   llamshift = fltarr(nshift)
   chi2vec = fltarr(nshift)

   ; Get the 1st line from the original spectrum
   ipix1 = where(rloglam GE llam1 - hwidth $
    AND rloglam LT llam1 + hwidth, npix1)
   ipix2 = where(rloglam GE llam2 - hwidth $
    AND rloglam LT llam2 + hwidth, npix2)

   for ishift=0, nshift-1 do begin

      ; Get the 2nd line from the over-sampled spectrum that most
      ; closely matches the wavelengths of the 1st line
      ibigpix = ipix2 * nsamp + bigpixshift[ishift]
      llamshift[ishift] = bigloglam[ibigpix[0]] - loglam[ipix1[0]]

      ; Trim to only good data points
      igood = where(objivar[ipix1] NE 0 AND bigivar[ibigpix] NE 0, ngood)
      itrim = ipix1[igood]
      ibigpix = ibigpix[igood]

      ; THE "MODEL" SPECTRUM: bigloglam[ibigpix], bigflux[ibigpix]
      ; THE "DATA" SPECTRUM: loglam[itrim], objflux[itrim]
      ; Approximate the variance in the "DATA" spectrum as:
      ;                      1/objivar[itrim]^2 + 0.5/bigivar[ibigpix]^2
      sigma2 = 1.0 / (objivar[itrim])^2 + 0.5/(bigivar[ibigpix])^2
      weight = 1.0 / sqrt(sigma2)

      ; Construct the basis vectors, which are polynomial terms
      ; and the shifted 5007 line.
      basisvec = [[poly_array(npix1,npoly)], [bigflux[ibigpix]]]
      mmatrix = fltarr(npix1,npoly+1)
      for i=0, npoly do mmatrix[*,i] = basisvec[*,i] * weight
      mtrans = transpose(mmatrix)
      bvec = objflux[itrim] * weight
      theta = invert(mtrans # mmatrix, /double) # mtrans # bvec
      thisfit = theta ## basisvec

      ; The following two lines are equivalent evaluations of chi^2
      chi2vec[ishift] = total( (objflux[itrim] - thisfit)^2 / sigma2)
;      chi2vec[ishift] = total( (mmatrix#theta - bvec)^2 )

      if (keyword_set(doplot)) then begin
         plot, loglam[itrim], objflux[itrim], $
          xtitle='Log10(Wavelength [Ang])', ytitle='Flux', psym=-4
         djs_oplot, loglam[itrim], thisfit, color='green'
         wait, wtime
      endif
   endfor

   ;----------
   ; Now fit for the best redshift using the chi^2

   ; Fit around the 4 points nearest to the minima.
   bestshift = find_nminima(chi2vec, llamshift, width=4*pixscale/nsamp, $
    xerr=besterr, errcode=errcode, doplot=doplot, $
    xtitle='Line separation in log10(Ang)', $
    plottitle=plottitle)

end
;------------------------------------------------------------------------------
