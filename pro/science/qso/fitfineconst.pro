;------------------------------------------------------------------------------
;+
; NAME:
;   fitfineconst
;
; PURPOSE:
;   Steinhardt's project to fit the fine structure constant from O_III lines
;
; CALLING SEQUENCE:
;   fitfineconst, plate, fiber, mjd=, [ linewave=, lineratio=, $
;    bestshift=, besterr=, bestratio=, bestrchi2=, /debug, /doplot ]
;
; INPUTS:
;   plate      - Plate number
;   fiber      - Fiber number
;
; REQUIRED KEYWORDS:
;   mjd        - MJD number for this PLATE
;
; OPTIONAL INPUTS:
;   linewave   - Two-element array with the line wavelengths (in air).
;                The fitting is done by shifting the 2nd line to fit the 1st.
;                Default to [4958.911d0, 5006.843d0] Ang.
;   lineratio  - Input flux ratio of (line 2):(line 1) for scaling the
;                errors; default to 2.0.
;   debug      - If set, then send debugging plots to the screen
;   doplot     - If set, then output a bunch of plots
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   bestshift  - Best shift between lines in log10(Ang)
;   besterr    - Formal error in BESTSHIFT, but rescaled by BESTCHI2
;   bestratio  - Fit to flux ratio (line 2):(line 1)
;   bestrchi2  - Reduced chi^2 of the best fit
;
; COMMENTS:
;
; EXAMPLES:
;   IDL> fitfineconst, 276, 251, mjd=51909, /debug, /doplot, $
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
function fitfineconst1, ishift, fitflux=fitflux, theta=theta

   common com_fitfine, npix1, npoly, lineratio, $
    nsamp, ipix1, ipix2, itrim, sigma2, bigpixshift, llamshift, $
    loglam, objflux, objivar, $
    bigloglam, bigflux, bigivar

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
;   sigma2 = 1. / (objivar[itrim])^2 + (1./lineratio^2) / (bigivar[ibigpix])^2
   sigma2 = 1. / (objivar[itrim]) + (1./lineratio^2) / (bigivar[ibigpix])
   weight = 1. / sqrt(sigma2)

   ; Construct the basis vectors, which are polynomial terms
   ; and the shifted 5007 line.
;   basisvec = [[poly_array(npix1,npoly)], [bigflux[ibigpix]]]
   basisvec = [[bigflux[ibigpix]], [poly_array(npix1,npoly)]]
   mmatrix = fltarr(npix1,npoly+1)
   for i=0, npoly do mmatrix[*,i] = basisvec[*,i] * weight
   mtrans = transpose(mmatrix)
   bvec = objflux[itrim] * weight
   theta = invert(mtrans # mmatrix, /double) # mtrans # bvec
   fitflux = theta ## basisvec

   ; The following two lines are equivalent evaluations of chi^2
   thischi2 = total( (objflux[itrim] - fitflux)^2 / sigma2)
;   thischi2 = total( (mmatrix#theta - bvec)^2 )

   return, thischi2
end
;------------------------------------------------------------------------------
pro fitfineconst, plate, fiber, mjd=mjd, debug=debug, doplot=doplot, $
 linewave=linewave, lineratio=lineratio1, $
 bestshift=bestshift, besterr=besterr, bestratio=bestratio, bestrchi2=bestrchi2

   common com_fitfine, npix1, npoly, lineratio, $
    nsamp, ipix1, ipix2, itrim, sigma2, bigpixshift, llamshift, $
    loglam, objflux, objivar, $
    bigloglam, bigflux, bigivar

   if (n_elements(plate) NE 1 OR n_elements(fiber) NE 1 $
    OR n_elements(mjd) NE 1) then begin
      doc_library, 'fitfineconst'
      return
   endif

   if (NOT keyword_set(linewave)) then linewave = [4958.911d0, 5006.843d0]
   if (keyword_set(lineratio1)) then lineratio = lineratio1 $
    else lineratio = 2.0

   plottitle = string(plate, fiber, mjd, $
    format='("Plate ",i4," Fiber ",i3," MJD ", i5)')
   pixscale = 1.d-4 ; pixel scale in log10(lambda)
   hwidth = 10 * pixscale ; 10 pixels half-width for fitting domain
   npoly = 2 ; Number of polynomial terms for continuum-fitting
   maxshift = 2.0 ; Maximum pixel shift to search
   nsamp = 20 ; Sub-sampling factor for pixels
   wtime = 0.05 ; Pause time in seconds between plots

   ;----------
   ; Set the wavelengths to fit -- the O_III lines at 4959,5007

   wave1 = linewave[0]
   wave2 = linewave[1]
   airtovac, wave1
   airtovac, wave2
   llam1 = alog10(wave1)
   llam2 = alog10(wave2)
   linesep = (llam2 - llam1)

   ;----------
   ; Read the spectrum and its redshift

   readspec, plate, fiber, mjd=mjd, $
    loglam=loglam, flux=objflux, invvar=objivar, zans=zans, zline=zline

   ; If ZLINE is set, then use the redshift of the line
   thisz = zans.z
   if (keyword_set(zline)) then begin
      ii = where(abs(zline.linewave - wave1) LT 3)
      if (ii[0] NE -1) then thisz = zline[ii[0]].linez
   endif
   print, 'Using z=', thisz

   ; Compute the rest-frame wavelengths...
   rloglam = loglam - alog10(1. + thisz)

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
      chi2vec[ishift] = fitfineconst1(ishift, fitflux=fitflux)

      if (keyword_set(debug)) then begin
         plot, loglam[itrim], objflux[itrim], $
          xtitle='Log10(Wavelength [Ang])', ytitle='Flux', psym=-4
         djs_oplot, loglam[itrim], fitflux, color='green'
         wait, wtime
      endif
   endfor

   ;----------
   ; Now fit for the best redshift using the chi^2

   ; Fit around the 4 points nearest to the minima.
   bestshift = find_nminima(chi2vec, llamshift, width=4*pixscale/nsamp, $
    xerr=besterr, ypeak=bestchi2, errcode=errcode, doplot=doplot, $
    xtitle='Line separation in log10(Ang)', $
    plottitle=plottitle)
   if (keyword_set(debug)) then begin
      print, 'Press any key...'
      cc = get_kbrd(1)
   endif

   ; Re-scale the errors by the chi^2 at the best fit
   bestrchi2 = bestchi2 / (npix1 - 3)
   besterr = besterr * sqrt(bestrchi2)

   junk = min( abs(llamshift - bestshift), ishift )
   junk = fitfineconst1(ishift, fitflux=fitflux, theta=theta)
   bestratio = 1. / theta[0]
   if (keyword_set(debug) OR keyword_set(doplot)) then begin
      plot, loglam[itrim], objflux[itrim], $
       title=plottitle, $
       xtitle='Log10(Wavelength [Ang])', ytitle='Flux', psym=-4
      djs_oploterr, loglam[itrim], objflux[itrim], yerr=sqrt(sigma2)
      djs_oplot, loglam[itrim], fitflux, color='green'
      if (keyword_set(debug)) then begin
         print, 'Press any key...'
         cc = get_kbrd(1)
      endif
   endif

   print, 'BESTSHIFT = ', bestshift
   print, 'BESTERR   = ', besterr
   print, 'BESTRATIO = ', bestratio
   print, 'BESTRCHI2 = ', bestrchi2
end
;------------------------------------------------------------------------------
