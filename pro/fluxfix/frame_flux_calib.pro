
function frame_flux_calib, wave, corvector, corvivar, avgcorvset, cormed, $
         framename, fit_wiggles = fit_wiggles, fsig = fsig, fcor = fcor, $
         final = final, noplot = noplot

  ;----------------
  ; Create 2d wave and flux calib vector average

  npix = n_elements(corvector[*,0])
  nstd = n_elements(corvector[0,*])
  wave2d = rebin(wave, npix, nstd)

  ;-----------------
  ; Divide each vector by the average and fit a polynomial to the
  ; result -- this is the "low frequency" part

  residcorv = corvector
  residivar = corvivar
 
  if not keyword_set(final) then begin
    avgcorv = bspline_valu(wave2d, avgcorvset)  ; average from PCA
    divideflat, residcorv, invvar=residivar, avgcorv, $
                minval=0.001*median(avgcorv)
  endif else avgcorv = corvector * 0 + 1

  xy2traceset, wave2d, residcorv, polyset, invvar=residivar, $
    ncoeff=4, lower = 2.5, upper = 2.5, maxiter = 15

  ;----------------
  ; Create an average low-frequency vector
  traceset2xy, polyset, wave2d, lowfvect
  lowfset = {func: 'legendre', xmin: polyset.xmin, xmax: polyset.xmax, $
             coeff: djs_median(polyset.coeff, 2)}
  traceset2xy, lowfset, wave, lowfmed

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

  ;----------------
  ; Check if scale factor is ok, this time using more pixels than when
  ; "cormed" was computed 

  ;frameresid = corvector * (1.0/fcor # cormed)

  ;if min(10.0^wave) lt 5000 then $
  ;    range = where(10.0^wave gt 5000 and 10.0^wave lt 6000) $
  ;else range = where(10.0^wave gt 6000 and 10.0^wave lt 7000)
  ;if keyword_set(final) then $
  ;  range = where(10.0^wave gt 5600 and 10.0^wave lt 6900)

  ;zptadjust_factor = median(frameresid[range, *])
  ;djs_iterstat, frameresid[range, *], mean=zptadjust_factor, sigrej=2.5
  ;splog, 'Zeropoint adjustment: ', string(zptadjust_factor, format='(F7.3)')
  ;fcor = fcor * zptadjust_factor
  
  ;plot, 10.0^wave, fcor, xr=[3800, 9200], /xs, /nodata, yr=[0.6, 1.4]
  ;for i = 0, nstd - 1 do oplot, 10.0^wave, frameresid[*,i]
  ;djs_oplot, [3000, 10000], [1, 1], color='blue'
  ;djs_oplot, [3000, 10000], [zptadjust_factor, zptadjust_factor], color='red'

  ;--------------
  ; Do the spline fit

  ;if keyword_set(final) then begin
     ;bbkptfile = filepath('blue.bkpts', $
     ;  root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
     ;readcol, bbkptfile, bbkpts, silent=1
     ;rbkptfile = filepath('red.bkpts', $
     ;  root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
     ;readcol, rbkptfile, rbkpts, silent=1
     ;bkpts = alog10([bbkpts, rbkpts]) ; Convert to log-10 Angstroms
     ;ibk = where(bkpts GE min(wave) AND bkpts LE max(wave), nbk)
     ;if (nbk LT 4) then splog, 'Error selecting break points'
     ;bkpts = [min(wave), bkpts[ibk], max(wave)]
  ;endif else bkpts = avgcorvset.fullbkpt

  if keyword_set(final) then nbkpts = 100 else nbkpts = 50
  bkpts = bspline_bkpts(wave, nord=4, nbkpts=nbkpts)

  calibset = bspline_iterfit(wave, fcor, bkpt = bkpts, $
             nord=4, upper=3, lower=3, maxrej=ceil(0.10*n_elements(fcor)))

  calibfac = bspline_valu(wave, calibset)

  ;------------
  ; QA plots

  if keyword_set(noplot) then return, calibset

  !P.MULTI = [0, 1, 2]
   
  djs_plot, 10.0^wave, hifvect, yr=[0.7, 1.2], /nodata, $
            xtitle = 'Wavelength', ytitle = 'Normalized Flux', $
            title = 'High Frequency Sphoto Correction to ' + framename
  for istd = 0, nstd - 1 do djs_oplot, 10.0^wave, hifvect[*,istd] 
  djs_oplot, 10.0^wave, hifmed, color='red', thick=3

  djs_plot, 10.0^wave, residcorv, yr=[0.5, 1.5], /nodata, $
            xtitle = 'Wavelength', ytitle = 'Normalized Flux', $
            title = 'Residual Sphoto Correction to ' + framename 
  for istd = 0, nstd - 1 do djs_oplot, 10.0^wave, residcorv[*,istd] 
  djs_oplot, 10.0^wave, lowfmed, color='red', thick=3
  djs_oplot, 10.0^wave, calibfac/avgcorv[*,0]/scalefactor, color='green', $
     thick=3
 
  xyouts, mean(10.0^wave) - 500, 1.3, 'Sigma = ' + $
    string(fsig * 100, format='(I3)') + '%'
  !P.MULTI = 0
  
  return, calibset

end
