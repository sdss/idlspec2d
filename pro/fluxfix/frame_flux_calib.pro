;+
; NAME:
;   frame_flux_calib
;
; PURPOSE:
;   To derive the flux calibration appropriate to one frame.  This routine
;   takes as input a set of "flux correction vectors" which are formed by 
;   ratioing each standard with it's best fit model and smoothing the result.
;   The average correction vector, created from a PCA of all the frames (of
;   one of the CCDs b1/r1/b2/r2) is first divided out.  This gets rid of the
;   low order shape.  The residuals are then fit with a 4th order legendere.
;   The average low order residual is found by medianing the Legendere 
;   coefficients.  High order terms are fit by first dividing out the low order
;   terms and then medianing the result.  These are included in the final fit
;   if "fit_wiggles" is set.  The flux calibration vecotor for the frame is
;   constructed by multiplying the average for all frames, the lower order 
;   residual fit to this frame, and the high order residual fit (if desired).  
;   This vector is then fit with a bspline and the coefficients are returned.
;
; CALLING SEQUENCE:
;   calibset = frame_flux_calib(loglam, corvector, corvivar, avgcorvset, $
;     cormed, framename, fit_wiggles=, fsig=, fcor=, final=, noplot=)
;
; INPUTS:
;   loglam     - wavelength array (log10(Angstroms)) [npix]
;   corvector  - array of "flux correction vectors" -- these are
;                formed by ratioing each standard with it's best fit model
;                and smoothing the result [npix, nstar]
;   corvivar   - inverse variance of correction vectors -- used in rejecting
;                points from low order fit [npix, nstar]
;   avgcorvset - a set of bspline coefficients which describes the average
;                flux correction of all frames. Not needed if "final" is set 
;   cormed     - The median value of each raw correction vector.  Before being
;                fed into this code they are all normalized at 5700 -- 6300 A 
;                (see spdata2model) [nstar]
;
; OUTPUT:
;   A structure holding the coefficients of a bspline fit to the average
;   flux correction vector for the frame is returned.  Unless "noplot" is set,
;   diagnostic plots are created.  If "outfile" is set the structure is stored
;   as a FITS binary file.
;
; KEYWORDS:
;   framename   - Name of frame (to be used in plot titles)
;   outfile     - Name of output FITS binary file -- if not set no file is 
;                 written
;   fit_wiggles - if set a fit is done to the high order (as well as low
;                 order) residuals.  This should only be set when the S/N
;                 is good.  Note that the high order fit is plotted even if 
;                 it is not used.
;   fsig        - the variance between each of the flux correction vectors
;                 and the average -- this is a good measure of the spectro-
;                 photometric quality
;   fcor        - average flux correction for the frame (not the spline fit)
;   final       - if set, the average flux correction from a PCA of all the
;                 frames is not used -- this is useful for a final evaluation
;                 of the spectrophotometric quality of the combined exposures
;   noplot      - toggles plotting
;
; COMMENTS:
;
; BUGS:
;   Not sure if it is mathematically OK to median the Legendre coeffs! (But
;   it seems to work well ...)
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   bspline_bkpts
;   bspline_iterfit
;   bspline_valu
;   divideflat
;   djs_icolor()
;   djs_iterstat()
;   djs_median()
;   djs_oplot
;   djs_plot
;   mwrfits
;   splot
;   sxaddpar
;   traceset2xy
;   xytraceset
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   12-Aug-2003 Created by C. Tremonti Steward Observatory
;-
;------------------------------------------------------------------------------

function frame_flux_calib, loglam, corvector, corvivar, avgcorvset, cormed, $
         framename = framename, outfile = outfile, noplot = noplot, $
         fit_wiggles = fit_wiggles, fsig = fsig, fcor = fcor, final = final

  npix = n_elements(corvector[*,0])
  nstd = n_elements(corvector[0,*])
  loglam2d = rebin(loglam, npix, nstd)
  wave = 10.0^loglam 

  ;-----------------
  ; Divide each vector by the average and fit a 4th order legendre to the
  ; result -- this is the "low frequency" part

  residcorv = corvector
  residivar = corvivar
 
  if not keyword_set(final) then begin
    avgcorv = bspline_valu(loglam2d, avgcorvset)  ; average from PCA
    divideflat, residcorv, invvar=residivar, avgcorv, $
                minval=0.001*median(avgcorv)
  endif else avgcorv = corvector * 0 + 1

  xy2traceset, loglam2d, residcorv, polyset, invvar=residivar, $
    ncoeff=4, lower = 2.5, upper = 2.5, maxiter = 15

  ;----------------
  ; Create an average low-frequency vector
  traceset2xy, polyset, loglam2d, lowfvect
  lowfset = {func: 'legendre', xmin: polyset.xmin, xmax: polyset.xmax, $
             coeff: djs_median(polyset.coeff, 2)}
  traceset2xy, lowfset, loglam, lowfmed

  ;----------------
  ; Create an average hi-frequency vector
  hifvect = residcorv
  divideflat, hifvect, lowfvect, minval=0.01
  hifmed = djs_median(hifvect, 2)

  ;--------------
  ; Recombine high & low-frequencey parts (but only use high-f if allowed  
  ; by keyword)  

  fcor = lowfmed * avgcorv[*,0]
  if fit_wiggles then fcor = fcor * hifmed

  ;--------------
  ; Measure the variance between the fluxcalib vectors derived
  ; for each of the standard stars -- this is an indicator of the
  ; spectrophotometric quality.

  djs_iterstat, corvector / rebin(fcor, npix, nstd), mean=fmean, sigma=fsig, $
    maxiter=5, sigrej=5

  ;---------------
  ; Restore to original scale

  djs_iterstat, cormed, mean=scalefactor, sigma=cormedsig, sigrej=2.5

  ;splog, 'Zeropoint of corvector: ', string(cormed, format='(F6.2)')
  splog, 'Average zeropoint & stdev: ' + $
    string(scalefactor, format='(F6.2)') + ' +/- ' + $
    string(cormedsig, format='(F5.2)')

  fcor = fcor * scalefactor 

  ;--------------
  ; Do the spline fit

  if keyword_set(final) then nbkpts = 100 else nbkpts = 50
  padbkpts = bspline_bkpts(loglam, nord=4, nbkpts=nbkpts, bkpt=bkpts, /silent)

  calibset = bspline_iterfit(loglam, fcor, bkpt = bkpts, $
             nord=4, upper=3, lower=3, maxrej=ceil(0.10*n_elements(fcor)))

  calibfac = bspline_valu(loglam, calibset)

  ;--------------
  ; Create header cards describing the data range and write to FITS
  ; The keyword 'SPHOTERR' holds the standard deviation of the
  ; correction vectors for the individual stars -- this is a good measure
  ; of the quality of the spectrophotometry

  if keyword_set(outfile) then begin
    hdr = ['']
    sxaddpar, hdr, 'WAVEMIN', min(wave)
    sxaddpar, hdr, 'WAVEMAX', max(wave)
    sxaddpar, hdr, 'SPHOTERR', fsig
    mwrfits, 0, outfile, hdr, /create
    mwrfits, calibset, outfile
  endif

  ;------------
  ; QA plots

  if keyword_set(noplot) then return, calibset

  !P.MULTI = [0, 1, 2]
  if not keyword_set(framename) then framename = ''

  djs_plot, wave, hifvect, yr=[0.7, 1.2], /nodata, $
            xtitle = 'Wavelength', ytitle = 'Normalized Flux', $
            title = 'High Frequency Sphoto Correction to ' + framename
  for istd = 0, nstd - 1 do djs_oplot, wave, hifvect[*,istd] 
  djs_oplot, wave, hifmed, color='red', thick=3

  djs_plot, wave, residcorv, yr=[0.5, 1.5], /nodata, $
            xtitle = 'Wavelength', ytitle = 'Normalized Flux', $
            title = 'Residual Sphoto Correction to ' + framename 
  for istd = 0, nstd - 1 do djs_oplot, wave, residcorv[*,istd] 
  djs_oplot, wave, lowfmed, color='red', thick=3
  djs_oplot, wave, calibfac/avgcorv[*,0]/scalefactor, color='green', $
     thick=3
 
  xyouts, mean(wave) - 500, 1.3, 'Sigma = ' + $
    string(fsig * 100, format='(I3)') + '%'
  !P.MULTI = 0
  
  return, calibset

end
