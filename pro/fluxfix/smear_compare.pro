
function smear_compare, smearname, finalwave, sciflux, sciivar, $
         best_exp_str, plate_str, mjd_str, $
         camnames = camnames, combinedir = combinedir, $
         tsobjname = tsobjname, nord=nord, maxsep=maxsep, adderr=adderr, $
         noplot = noplot

  ;------------------------------------
  ; Read frame files into big 2-d array

  spread_frames, smearname, window=window, binsz = binsz, $
     adderr=adderr, camnames=camnames, tsobjname = tsobjname, $
     flux = flux, ivar = fluxivar, wave = wave, pixelmask = pixelmask, $
     plugtag = plugtag, camerasvec = camerasvec, expid = expid, sn2 = sn2,  $
     hdrarr = hdrarr

  splog, 'Calibrating smear frames:', expid[uniq(expid)]

  ;---------------------------------------------
  ; Now determine the spectrophotometric solution

  if not keyword_set(camnames) then camnames = ['b1', 'b2', 'r1', 'r2']
  ncam = N_elements(camnames)

  ; Loop through cameras
   for icam=0, ncam-1 do begin
     camid = camnames[icam]
     camcol = strmid(camnames[icam], 0,1)
     specnum = strmid(camnames[icam], 1,1)
     frames = expid[where(camerasvec eq camid)]

     ; File names for sphoto calibration inputs/outputs
     fcalfiles = djs_filepath('spFluxcalib' + '-' + camid + '-' + $
                 frames + '.fits', root_dir=combinedir)

     stdstarfile = djs_filepath('spStd-' + plate_str + '-' + mjd_str +  $
                   '-' + specnum + '.fits', root_dir=combinedir)

     avgcalfile = djs_filepath('spFluxcalib' + '-' + camid + '-' + $
                  best_exp_str + '.fits', root_dir=combinedir)

     ;--------------------------------
     ; Use plugmap to find standard stars

     nobadmask = reform(qgoodfiber(pixelmask[0,*]))

     isphoto = where((strtrim(plugtag.objtype) EQ 'SPECTROPHOTO_STD' OR $
               strtrim(plugtag.objtype) EQ 'REDDEN_STD') AND $
               (plugtag.spectrographid eq specnum) AND $
               (plugtag.camcolor eq camcol) AND nobadmask, nstd)

     incalib = mrdfits(avgcalfile, 1)

     sphoto_calib, wave[*,isphoto], flux[*,isphoto], fluxivar[*,isphoto], $
       pixelmask[*,isphoto], plugtag[isphoto], fcalfiles, $
       stdstarfile, input_calibset=incalib

     ;---------------------------------
     ; Apply sphoto calibration to all fibers in each frame

     for iframe = 0, n_elements(frames) - 1 do begin

       indx = where(plugtag.expid eq frames[iframe] AND $
                    plugtag.camcolor eq camcol AND $
                    plugtag.spectrographid eq specnum)

       junk = mrdfits(fcalfiles[iframe], 0, calibhdr, /silent)
       calibset = mrdfits(fcalfiles[iframe], 1)

       cwavemin = sxpar(calibhdr, 'WAVEMIN')
       cwavemax = sxpar(calibhdr, 'WAVEMAX')
       calibfac = bspline_valu(wave[*,indx], calibset)

       ; Set to bad any pixels whose wavelength is outside the known
       ; flux-calibration region.
       ibad = where(wave[*,indx] LT alog10(cwavemin) OR $
                    wave[*,indx] GT alog10(cwavemax))
       if (ibad[0] NE -1) then calibfac[ibad] = 0

       tempflux = flux[*,indx]
       tempivar = fluxivar[*,indx]
       
       divideflat, tempflux, invvar=tempivar, calibfac, $
                   minval=0.001*mean(calibfac)

       flux[*,indx] = tempflux
       fluxivar[*,indx] = tempivar
       pixelmask[*,indx] = pixelmask[*,indx] $
        OR (calibfac LE 0.001*mean(calibfac)) * pixelmask_bits('BADFLUXFACTOR')
     endfor
  endfor

  ;-----------------------
  ; Combine frames

  nfinalpix = n_elements(finalwave)
  nfiber = max(plugtag.fiberid)
  finalflux = fltarr(nfinalpix, nfiber)
  finalivar = fltarr(nfinalpix, nfiber)
  finalandmask = lonarr(nfinalpix, nfiber)
  finalplugtag = replicate(plugtag[0], nfiber)
   struct_assign, {fiberid: 0L}, finalplugtag
  binsz = finalwave[1] - finalwave[0]

  splog, 'Combining smear exposures'

  for ifiber=0, nfiber-1 do begin

    ; Find the first occurance of fiber number IFIBER+1
    indx = where(plugtag.fiberid EQ ifiber+1)

    if (indx[0] NE -1) then begin
;     splog, 'Fiber', ifiber+1, ' ', plugtag[indx[0]].objtype, $
;     plugtag[indx[0]].mag, format = '(a, i4.3, a, a, f6.2, 5f6.2)'

      finalplugtag[ifiber] = plugtag[indx[0]]
      temppixmask = pixelmask[*,indx]

      combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
          finalmask=temppixmask, andmask=bestandmask, $
          newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
          nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
          maxiter=50, upper=3.0, lower=3.0, maxrej=1

      finalflux[*,ifiber] = bestflux
      finalivar[*,ifiber] = bestivar
      finalandmask[*,ifiber] = bestandmask

    endif else begin
      splog, 'Fiber', ifiber+1, ' NO DATA'
      finalandmask[*,ifiber] = pixelmask_bits('NODATA')
    endelse
  endfor 

  ; Free memory
  wave = 0
  flux = 0
  fluxivar = 0

  ;---------------------------------------------------------------------------
  ; Compare smear & final science image
  ;---------------------------------------------------------------------------

  wave2d = finalwave # replicate(1, nfiber)

  smear_medflux = spmedian_rebin(wave2d, finalflux, finalivar, 'full', $
                  outwave = medwave, sn = smear_sn, mask = smear_mask, $
                  quality = smear_quality) 
 
  sci_medflux = spmedian_rebin(wave2d, sciflux, sciivar, 'full', $
                sn = sci_sn, mask = sci_mask, quality = sci_quality) 

  ; Fit with 4th order Legendre polynomial  
  medwave2d = float(medwave) # replicate(1, nfiber)

  xy2traceset, medwave2d, smear_medflux, smearset, invvar=smear_mask, $
    ncoeff=4, inputfunc=sci_medflux, lower=3, upper=3

  ; Save as binary table with science & smear S/N ratios
  smear_struct = {sci_sn: 0.0, smear_sn: 0.0, legendre_coeff: fltarr(4)}

  smear_hdu = make_array(dim=nfiber, value=smear_struct)
  smear_hdu.legendre_coeff = smearset.coeff
  smear_hdu.sci_sn = sci_sn
  smear_hdu.smear_sn = smear_sn
 
  ;---------------------------------------------------------------------------
  ; QA plot
  ;---------------------------------------------------------------------------

  ; Only use where science & smear have reasonable S/N 
  ok = sci_sn GT 3.0 AND smear_sn GT 1.0 AND $
       sci_quality EQ 0 AND smear_quality EQ 0
  traceset2xy, smearset, wave2d, smear_ratio

  ;!P.MULTI = [0, 1, 2]
  ;for i = 0, nfiber - 1 do begin
  ;  plot, 10.0^finalwave, smooth(finalflux[*,i], 5), xr=[3800, 9200], /xs
  ;  djs_oplot, 10.0^finalwave,  smooth(sciflux[*,i], 5), color='red'
  ;  smsci = smooth(djs_median(sciflux[*,i], width=75), 25)
  ;  djs_oplot, 10.0^finalwave,  smsci * smear_ratio[*,i], color='green'
  ;endfor
  ;!P.MULTI = 0
  
  ptsrc = where(strmatch(finalplugtag.objtype, '*GALAXY*') ne 1 and ok, nptsrc)
  qso = where(strmatch(finalplugtag.objtype, '*QSO*') eq 1 and ok, nqso)

  linwave=10.0^finalwave
  if not keyword_set(noplot) then begin
    djs_plot, linwave, smear_ratio[*,0], xr=[3800, 9200], /xs, $
      yr=[0, 2.5], /ys, xtitle = 'Wavelength', ytitle = 'Smear / Science', $
      title = 'Smear Correction Vectors for Point Sources', /nodata
    for iobj = 0, nptsrc - 1 do $
      djs_oplot, linwave, smear_ratio[*,ptsrc[iobj]], nsum=10, color='green'
    for iobj = 0, nqso - 1 do $
      djs_oplot, linwave, smear_ratio[*,qso[iobj]], nsum=10, color='blue'
    djs_oplot, [2000, 10000], [1, 1], thick= 4
    legend, ['Quasars', 'Other Point Sources'], $
            color=djs_icolor(['blue', 'green']), psym=[0,0]
  endif

  splog, 'Median smear correction for point sources: ' +  $
          string(djs_median(smear_ratio[*,ptsrc]), format = '(F6.2)')

  return, smear_hdu
  
end
