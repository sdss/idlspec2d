
function frame_flux_calib, wave, corvector, corvivar, avgcorvset, cormed, $
         bkpts, framename, fsig = fsig, fcor = fcor, noplot = noplot, $
         blue = blue

  ;----------------
  ; Create 2d wave and flux calib vector average

  npix = n_elements(corvector[*,0])
  nstd = n_elements(corvector[0,*])
  wave2d = rebin(wave, npix, nstd)

  avgcorv = bspline_valu(wave2d, avgcorvset)  ; average from PCA
 
  ;-----------------
  ; Divide each vector by the average and fit a polynomial to the
  ; result -- this is the "low frequency" part

  residcorv = corvector
  residivar = corvivar
  divideflat, residcorv, invvar=residivar, avgcorv, $
    minval=0.001*median(avgcorv)
  xy2traceset, wave2d, residcorv, polyset, invvar=residivar, $
    ncoeff=3, lower = 2.5, upper = 2.5, maxiter = 15

  ;----------------
  ; Create an average low-frequency vector
  traceset2xy, polyset, wave2d, lowfvect
  lowfset = {func: 'legendre', xmin: polyset.xmin, xmax: polyset.xmax, $
             coeff: djs_median(polyset.coeff, 2)}
  traceset2xy, lowfset, wave, lowfmed

  ;----------------
  ; Create an average hi-frequency vector
  hifvect = residcorv
  divideflat, hifvect, lowfvect, minval=0.1
  hifmed = djs_median(hifvect, 2)

  ;--------------
  ; Recombine high & low-frequencey parts and the average 

  fcor = lowfmed * hifmed * avgcorv[*,0]

  ;--------------
  ; Measure the variance between the fluxcalib vectors derived
  ; for each of the standard stars -- this is an indicator of the
  ; spectrophotometric quality.

  meanclip, corvector - rebin(fcor, npix, nstd), fmean, fsig, $
    maxiter=3, clipsig=5

  ;---------------
  ; Restore to original scale

  meanclip, cormed, scalefactor, cormedsig, maxiter=10, clipsig=3

  splog, 'Zeropoint of corvector: ', string(cormed, format='(F6.2)')
  splog, 'Average zeropoint & stdev: ' + $
    string(scalefactor, format='(F6.2)') + ' +/- ' + $
    string(cormedsig, format='(F7.2)')

  fcor = fcor * scalefactor 

  ;----------------
  ; Original "cormed" (median of corvector) is measured where the sphoto
  ; correction vectors are relatively flat, but to optimize the match between
  ; the blue and red want to take the median near the dichroic

  frameresid = corvector * (1.0/fcor # cormed)

  if keyword_set(blue) then $
      range = where(10.0^wave gt 5000 and 10.0^wave lt 6000) $
  else range = where(10.0^wave gt 6000 and 10.0^wave lt 7000)
 
  fcor = fcor * median(frameresid[range, *])

  ;--------------
  ; Do the spline fit
  calibset = bspline_iterfit(wave, fcor, nord=4, bkpt=bkpts, $
             upper=3, lower=3, maxrej=ceil(0.10*n_elements(fcor)))

  calibfac = bspline_valu(wave, calibset)

  ;------------
  ; QA plots

  if keyword_set(noplot) then return, calibset

  !P.MULTI = [0, 1, 2]
  djs_plot, 10.0^wave, hifvect, yr=[0.7, 1.2], /nodata, $
            xtitle = 'Wavelength', ytitle = 'Normalized Flux', $
            title = 'High Frequency Sphoto Correction to Frame ' + framename
  for istd = 0, nstd - 1 do djs_oplot, 10.0^wave, hifvect[*,istd] 
  djs_oplot, 10.0^wave, hifmed, color='red', thick=3

  djs_plot, 10.0^wave, residcorv, yr=[0.5, 1.5], /nodata, $
            xtitle = 'Wavelength', ytitle = 'Normalized Flux', $
            title = 'Sphoto Correction to Frame ' + framename + $
                    ' / Average Sphoto Correction'
  for istd = 0, nstd - 1 do djs_oplot, 10.0^wave, residcorv[*,istd] 
  djs_oplot, 10.0^wave, lowfmed, color='red', thick=3
  djs_oplot, 10.0^wave, calibfac/avgcorv[*,0]/scalefactor, color='green', $
     thick=3
 
  xyouts, mean(wave) - 500, 1.3, 'Sigma = ' + $
    string(fsig * 100, format='(I3)') + '%'
  !P.MULTI = 0
  
  return, calibset

end
