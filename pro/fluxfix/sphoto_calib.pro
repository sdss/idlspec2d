
pro sphoto_calib, wave, flux, invvar, mask, plugtag, fcalfile, $
    stdstarfile, stype = stype, input_calibset = input_calibset, $
    pca_calibset = pca_calibset, noplot = noplot

  ; Double check that all spcetra are from 1 camera + spectrograph
  camid = plugtag.camcolor + strtrim(plugtag.spectrographid, 2)
  if n_elements(uniq(camid)) gt 1 then begin
    splog, 'ABORT: Camera mismatch in spectrophotometric calibration!'
    return
  endif else camid = camid[0]

  ;-----------------------------------------
  ; Figure out frame, fiber, & pixel dimensions 
  ; Note, we assume we are only combining 1 plugmap

  npix = n_elements(flux[*,0])
  nspec = n_elements(flux[0,*])
  fibers = plugtag[uniq(plugtag.fiberid, sort(plugtag.fiberid))].fiberid
  nfibers = n_elements(fibers)
  frames = plugtag[uniq(plugtag.expid, sort(plugtag.expid))].expid
  nframes = n_elements(frames)

  ;----------------------------------------------------------------------------
  ; Put everything on a common wavelength grid
  ;----------------------------------------------------------------------------

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
    invvari = skymask(invvari, maski, ngrow =3)
    fluxi = djs_maskinterp(fluxi, (invvari EQ 0), iaxis=0, /const)
 
    newflux[*,ispec] = fluxi
    newivar[*,ispec] = invvari
    newmask[*,ispec] = maski
  endfor

  ;----------------------------------------------------------------------------
  ; If the frames are blue, spectral type the standards 
  ;----------------------------------------------------------------------------

  if keyword_set(stype) then begin

    ;-----------------------------------------
    ; Rectify spectra and combine frames for each fiber  

    avgflux = fltarr(nnewpix, nfibers)
    avginvar = fltarr(nnewpix, nfibers)
    plug1 = replicate(plugtag[0], nfibers)
    struct_assign, {fiberid: 0L}, plug1  ; Zero all elements

    for ifib = 0, nfibers - 1 do begin
      indx = where(plugtag.fiberid eq fibers[ifib])
      plug1[ifib] = plugtag[indx[0]]
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

  ;----------------------------------------------------------------------------
  ; Read back standard star info and expand to match the full list of spectra
  ;----------------------------------------------------------------------------

  stdinfo = mrdfits(stdstarfile, 1)
  stdinfo_all = make_array(val = stdinfo[0], dim = nspec)
  struct_assign, {fiberid: 0L}, stdinfo_all  ; Zero out all elements

  for ispec = 0, nspec -1 do begin
    fiber_match = where(stdinfo.fiberid eq plugtag[ispec].fiberid)
    if fiber_match[0] ne -1 then stdinfo_all[ispec] = stdinfo[fiber_match]
    ; Unmatched fibers will occur if the blue side had bad fibers (bad trace,
    ; etc) but the red side was ok -- these will have S/N of zero so get 
    ; eliminated in the next step
  endfor

  ;---------------------
  ; Chose high S/N spectra (otherwise spectral typing is likely to be
  ; wrong).  The velocity cut (hopefully) eliminates galaxies and
  ; other things that accidentally get targeted

  ok = where(stdinfo_all.sn gt 20 and abs(stdinfo_all.v_off) lt 450 and $
             stdinfo_all.mag[2] gt 0, nok)
  nok = nok / nframes
  if nok lt 3 then begin
    splog, 'WARNING:  Too few spectrophoto standards with good S/N'
    splog, 'Proceeding anyway!'
    ok = where(stdinfo_all.sn gt 10 and abs(stdinfo_all.v_off) lt 450 and $
               stdinfo_all.mag[2] gt 0, nok)
    nok = nok / nframes
  endif
  if nok lt 3 then begin
    splog, 'WARNING: Too few spectrophoto standards!!!'
    splog, 'Proceeding anyway!'
    ok = where(stdinfo_all.sn gt 0 and stdinfo_all.mag[2] gt 0, nok)
    nok = nok / nframes
  endif else begin
    splog, 'Using ' + string(nok, format='(I2)') + ' spectrophoto standards'
  endelse

  ;----------------------------------------------------------------------------
  ; Compute the average spectrophotometric calibration from all frames 
  ;----------------------------------------------------------------------------

  pca_flux_standard, newwave, newflux[*,ok], newivar[*,ok], $
    stdinfo_all[ok], camid, corvector = corvector, corvivar = corvivar, $
    cormed = cormed, calibset = pca_calibset, bkpts = bkpts 

  ; If desired, replace PCA average with normalized input calibset 
  ; (This option is used for smears)
  if keyword_set(input_calibset) then begin
     input_calibfac = bspline_valu(newwave, input_calibset)
     iwave = 10.0^newwave
     norm_indx = where(iwave gt 5700 and iwave lt 6300 and $
                       iwave lt max(iwave) - 200 and $
                       iwave gt min(iwave) + 200)
     input_calibfac = input_calibfac / djs_median(input_calibfac[norm_indx])
     avgcalibset = bspline_iterfit(newwave, input_calibfac, $
                   bkpt = input_calibset.fullbkpt, nord = 4)
  endif else avgcalibset = pca_calibset

  ;----------------------------------------------------------------------------
  ; Compute the spectrophoto calibration separately for each frame
  ;----------------------------------------------------------------------------

  for iframe = 0, nframes - 1 do begin
    indx = where(plugtag[ok].expid eq frames[iframe], nindx)

    ; check to be sure that fluxcalib output name matches frame ID 
    if not strmatch(fcalfile[iframe], '*' + frames[iframe] + '*') then begin
      splog, 'ABORT: Fluxcalib file names do not match frame names' 
      return
    endif

    ; compute the median S/N of the sphot fibers in the frame
    sn = newflux[*,ok[indx]] * sqrt(newivar[*,ok[indx]])
    good = where(newivar[*,ok[indx]] ne 0) 
    medsn = median(sn[good])

    ; Find the median of the low and hi-frequency componets of the 
    ; corvectors separately -- use the hi-f part only if the S/N is good
    frame_calibset = frame_flux_calib(newwave, corvector[*,indx], $
      corvivar[*,indx], avgcalibset, cormed[indx], $
      'Frame ' + camid + '-' + frames[iframe], fit_wiggles = (medsn gt 12), $
      fsig = fsig, noplot = noplot)

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
