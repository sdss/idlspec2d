
pro sphoto_calib, wave, flux, invvar, mask, plugmap, frameid, $
    fcalfile, stdstarfile, combinedir = combinedir, blue = blue

  ;-----------------------------------------
  ; Figure out frame, fiber, & pixel dimensions (wave -> [npix, nfib*nframe])
  ; For now, we assume we are only combining 1 plugmap

  plug1 = plugmap[uniq(plugmap.fiberid, sort(plugmap.fiberid))]
  nfib = n_elements(plug1)
  frames = frameid[uniq(frameid)]
  nframes = n_elements(frames)
  npix = n_elements(flux[*,0])

  ;-----------------------------------------
  ; Put everything on a common wavelength grid

  if wave[0,0] lt wave[npix-1,0] then endpix=[0,npix-1] else endpix=[npix-1,0]
  wavemin = max(wave[endpix[0],*])
  wavemax = min(wave[endpix[1],*])
  dloglam = 1.0d-4

  nnewpix = round((wavemax - wavemin) / dloglam)
  newwave = findgen(nnewpix) * dloglam + wavemin  ; units are log-lambda
  wave2d = rebin(newwave, nnewpix, nfib * nframes)
  newflux = fltarr(nnewpix, nfib * nframes)
  newivar = fltarr(nnewpix, nfib * nframes)
  newmask = fltarr(nnewpix, nfib * nframes)

  for ispec = 0, nfib * nframes - 1 do begin

    combine1fiber, wave[*,ispec], flux[*,ispec], invvar[*,ispec], $
      finalmask=mask[*,ispec], newloglam=newwave, newflux=fluxi, $
      newivar=invvari, andmask=maski, binsz=dloglam

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
      indx = where(plugmap.fiberid eq plug1[ifib].fiberid)
      flux1fib = rectify(newflux[*,indx], newivar[*,indx], nivar = invar1fib)
      ; again to get rid of high order terms ...
      flux1fib = rectify(flux1fib)

      combine1fiber, wave2d[*,indx], flux1fib, invar1fib, $
        finalmask = newmask[*,indx], newloglam=newwave, newflux=avgflux1fib, $
        newivar=avginvar1fib, maxiter=50, upper=3.0, lower=3.0, maxrej=1, $
        binsz=dloglam, andmask = mask1fib

      avgflux[*,ifib] = avgflux1fib
      ; use mask to get rid of bad pix!!
      avginvar[*,ifib] = avginvar1fib
    endfor

    ;-----------------------------------------
    ; Spectral type the averaged fibers

    stype_standard, newwave, avgflux, avginvar, plug1, stdstarfile

  endif 

  ;-----------------------------------------
  ; Now for each frame compute the spectrophotometric calibration

  for iframe = 0, nframes - 1 do begin
     indx = where(frameid eq frames[iframe])

     ; check to be sure that fluxcalib output name matches frame ID 
     if not strmatch(fcalfile[iframe], '*' + frames[iframe] + '*') then begin
       splog, 'ABORT: Fluxcalib file names do not match frame names' 
       return
     endif

     pca_flux_standard, newwave, newflux[*,indx], newivar[*,indx], $
       newmask[*,indx], stdstarfile, outname = fcalfile[iframe], blue=blue 
  endfor

end
