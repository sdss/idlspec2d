
pro sphoto_calib, wave, flux, invvar, mask, fibertag, fcalfile, $
    stdstarfile, blue = blue

  ;-----------------------------------------
  ; Figure out frame, fiber, & pixel dimensions (wave -> [npix, nfib*nframe])
  ; For now, we assume we are only combining 1 plugmap

  plug1 = fibertag[where(fibertag.expid eq fibertag[0].expid)]
  nfib = n_elements(plug1)
  frames = fibertag[where(fibertag.fiberid eq fibertag[0].fiberid)].expid
  nframes = n_elements(frames)
  npix = n_elements(flux[*,0])
  nspec = nfib * nframes
  if n_elements(flux[0,*]) ne nspec then begin
     splog, 'ABORT: # of frames and fibers not equal to # of spectra!'
     return
  endif

  ;-----------------------------------------
  ; Put everything on a common wavelength grid

  if wave[0,0] lt wave[npix-1,0] then endpix=[0,npix-1] else endpix=[npix-1,0]
  wavemin = max(wave[endpix[0],*])
  wavemax = min(wave[endpix[1],*])
  dloglam = 1.0d-4

  nnewpix = round((wavemax - wavemin) / dloglam)
  newwave = findgen(nnewpix) * dloglam + wavemin  ; units are log-lambda
  wave2d = rebin(newwave, nnewpix, nspec)
  newflux = fltarr(nnewpix, nspec)
  newivar = fltarr(nnewpix, nspec)
  newmask = fltarr(nnewpix, nspec)

  for ispec = 0, nspec - 1 do begin

    combine1fiber, wave[*,ispec], flux[*,ispec], invvar[*,ispec], $
      finalmask=mask[*,ispec], newloglam=newwave, newflux=fluxi, $
      newivar=invvari, andmask=maski, binsz=dloglam

    ;-------------------
    ; Mask out bad pixels and regions dominated by sky-sub residuals
    invvari = skymask(invvari, maski)
    fluxi = djs_maskinterp(fluxi, (invvari EQ 0), iaxis=0, /const)
 
    newflux[*,ispec] = fluxi
    newivar[*,ispec] = invvari
    newmask[*,ispec] = maski
  endfor

  ;-----------------------------------------
  ; If the frames are blue, spectral type the standards 

  if keyword_set(blue) then begin

    ;-----------------------------------------
    ; Rectify spectra and combine frames for each fiber  

    avgflux = fltarr(nnewpix, nfib)
    avginvar = fltarr(nnewpix, nfib)

    for ifib = 0, nfib - 1 do begin
      indx = where(fibertag.fiberid eq plug1[ifib].fiberid)
      flux1fib = rectify(newflux[*,indx], newivar[*,indx], nivar = invar1fib)

      combine1fiber, wave2d[*,indx], flux1fib, invar1fib, $
        finalmask = newmask[*,indx], newloglam=newwave, newflux=avgflux1fib, $
        newivar=avginvar1fib, maxiter=50, upper=3.0, lower=3.0, maxrej=1, $
        binsz=dloglam

      avgflux[*,ifib] = avgflux1fib
      avginvar[*,ifib] = avginvar1fib
    endfor

    ;-----------------------------------------
    ; Spectral type the averaged fibers

    stdinfo = stype_standard(newwave, avgflux, avginvar, plug1, $
              outfile = stdstarfile)
  endif 

  stdinfo = mrdfits(stdstarfile, 1)

  ;-----------------------------------------
  ; Expand the standard star info to match the full list of frames

  stdinfo_all = make_array(val = stdinfo[0], dim = nspec)
  for ifib = 0, nfib * nframes -1 do $
      stdinfo_all[ifib] = stdinfo[where(stdinfo.fiberid eq $
                                        fibertag[ifib].fiberid)]

  ;---------------------
  ; Chose high S/N spectra (otherwise spectral typing is likely to be
  ; wrong).  The velocity cut (hopefully) eliminates galaxies and
  ; other things that accidentally get targeted

  ok = where(stdinfo_all.sn gt 20 and abs(stdinfo_all.v_off) lt 450, nok)
  nok = nok / nframes
  if nok lt 3 then begin
    splog, 'WARNING:  Too few spectrophoto standards with good S/N'
    splog, 'Proceeding anyway!'
    ok = where(stdinfo_all.sn gt 10 and abs(stdinfo_all.v_off) lt 450, nok)
    nok = nok / nframes
  endif
  if nok lt 3 then begin
    splog, 'WARNING:  NO good spectrophoto standards!!!'
    splog, 'Proceeding anyway!'
    ok = where(stdinfo_all.sn gt 0, nok)
    nok = nok / nframes
  endif else begin
    splog, 'Using ' + string(nok, format='(I2)') + ' spectrophoto standards'
  endelse

  ;-----------------------------------------
  ; Compute the spectrophotometric calibration from all standards on all 
  ; frames (for best results with PCA)

  pca_flux_standard, newwave, newflux[*,ok], newivar[*,ok], $
    stdinfo_all[ok], corvector = corvector, corvivar = corvivar, $
    cormed = cormed, calibset = calibset, bkpts = bkpts 

  ;-----------------------------------------
  ; Now for each frame compute the spectrophotometric calibration

  for iframe = 0, nframes - 1 do begin
    indx = where(fibertag[ok].expid eq frames[iframe])

    ; check to be sure that fluxcalib output name matches frame ID 
    if not strmatch(fcalfile[iframe], '*' + frames[iframe] + '*') then begin
      splog, 'ABORT: Fluxcalib file names do not match frame names' 
      return
    endif
   
    thiswave = wave2d[*,indx]
    avgcorv = bspline_valu(thiswave, calibset)
    corvector_frame = corvector[*,indx]
    corvivar_frame = corvivar[*,indx]
    ;divideflat, corvector_frame, invvar=corvivar_frame, avgcorv, $
    ;    minval=0.005*avgcorv
        
    ;pcaflux = pca_solve(corvector_frame, corvivar_frame, $
    ;    niter=30, nkeep=1, nreturn = 1, acoeff=acoeff, $
    ;    maxiter=5, upper=5, lower=5, $
    ;    eigenval = eigenval, usemask = usemask, $
    ;    maxrej=ceil(0.1*nnewpix), groupsize=ceil(nnewpix/5.))

    ;fcor = pcaflux * median(acoeff)
    

     xy2traceset, thiswave[*], corvector_frame[*], polyset, $
            invvar=corvivar_frame[*], ncoeff=3, inputfunc=avgcorv[*], $
            lower = 2.5, upper = 2.5, maxiter = 15

     traceset2xy, polyset, newwave, fcor

;stop
    ;--------------
    ; Measure the variance between the fluxcalib vectors derived
    ; for each of the standard stars -- this is an indicator of the
    ; spectrophotometric quality.

    meanclip, corvector_frame/avgcorv - rebin(fcor, nnewpix, nok), $
              fmean, fsig, maxiter=3, clipsig=5

    ;---------------
    ; QA plot

    djs_plot, 10.0^newwave, corvector_frame[*,0], /nodata, $
           xtitle='\lambda [A]', $
           ytitle='Counts / (10^{-17}erg/cm^{2}/s/A)', yr = [0.5, 1.5], $
           title = ' Spectrophoto Correction: ' + frames[iframe] + '-' + $
           fibertag[0].camcolor + strtrim(fibertag[0].spectrographid, 2)

    for istd=0, nok - 1 do oplot, 10.0^newwave, corvector_frame[*,istd]/avgcorv
    djs_oplot, 10.0^newwave, fcor, color='green', thick=3
    xyouts, mean(10.0^newwave) - 500, 1.3, $
    'Standard star variation = ' + string(fsig * 100, format='(I3)') + ' %'

    ;---------------
    ; Restore to original shape + scale
    ;scalefactor = median(cormed[indx]) 
    meanclip, cormed[indx], scalefactor, cormedsig, maxiter=10, clipsig=3
    splog, 'Zeropoint of corvector: ', cormed[indx]
    splog, 'Average zeropoint & stdev: ', scalefactor, cormedsig
    fcor = fcor * scalefactor * avgcorv[*,0]
    
    ;--------------
    ; Do the spline fit
    frame_calibset = bspline_iterfit(newwave, fcor, nord=4, bkpt=bkpts, $
             upper=3, lower=3, maxrej=ceil(0.10*n_elements(fcor)))
 
    calibfac = bspline_valu(newwave, frame_calibset)
    djs_oplot, newwave, calibfac/avgcorv[*,0]/scalefactor, color='red', thick=3

    ;--------------
    ; Create header cards describing the data range and write to FITS
    ; The keyword 'SPHOTERR' holds the standard deviation of the
    ; correction vectors for the individual stars -- this is a good measure
    ; of the quality of the spectrophotometry

    hdr = ['']
    sxaddpar, hdr, 'WAVEMIN', min(10.0^newwave)
    sxaddpar, hdr, 'WAVEMAX', max(10.0^newwave)
    sxaddpar, hdr, 'SPHOTERR', fsig
    mwrfits, 0, fcalfile[iframe], hdr, /create 
    mwrfits, frame_calibset, fcalfile[iframe]  
  endfor

end
