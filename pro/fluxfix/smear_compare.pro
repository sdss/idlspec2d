
function smear_compare, smearname, finalwave, sciflux, sciivar, $
         best_exp_str, plate_str, mjd_str, $
         camnames = camnames, combinedir = combinedir, $
         tsobjname = tsobjname, nord=nord, maxsep=maxsep, adderr=adderr

  ;------------------------------------
  ; Read frame files into big 2-d array

  spread_frames, smearname, window=window, binsz = binsz, $
     adderr=adderr, camnames=camnames, tsobjname = tsobjname, $
     flux = flux, ivar = fluxivar, wave = wave, pixelmask = pixelmask, $
     fibertag = fibertag, $
     camerasvec = camerasvec, expid = expid, sn2 = sn2,  $
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

     isphoto = where((strtrim(fibertag.objtype) EQ 'SPECTROPHOTO_STD' OR $
               strtrim(fibertag.objtype) EQ 'REDDEN_STD') AND $
               (fibertag.spectrographid eq specnum) AND $
               (fibertag.camcolor eq camcol), nstd)

     incalib = mrdfits(avgcalfile, 1)

     sphoto_calib, wave[*,isphoto], flux[*,isphoto], fluxivar[*,isphoto], $
       pixelmask[*,isphoto], fibertag[isphoto], fcalfiles, $
       stdstarfile, input_calibset=incalib

     ;---------------------------------
     ; Apply sphoto calibration to all fibers in each frame

     for iframe = 0, n_elements(frames) - 1 do begin

       indx = where(fibertag.expid eq frames[iframe] AND $
                    fibertag.camcolor eq camcol AND $
                    fibertag.spectrographid eq specnum)

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
  nfiber = max(fibertag.fiberid)
  finalflux = fltarr(nfinalpix, nfiber)
  finalivar = fltarr(nfinalpix, nfiber)
  finalandmask = lonarr(nfinalpix, nfiber)
  finalfibertag = replicate(fibertag[0], nfiber)
   struct_assign, {fiberid: 0L}, finalfibertag
  binsz = finalwave[1] - finalwave[0]

  splog, 'Combining smear exposures'

  for ifiber=0, nfiber-1 do begin

    ; Find the first occurance of fiber number IFIBER+1
    indx = (where(fibertag.fiberid EQ ifiber+1))[0]

    if (indx NE -1) then begin
;     splog, 'Fiber', ifiber+1, ' ', fibertag[indx].objtype, $
;     fibertag[indx].mag, format = '(a, i4.3, a, a, f6.2, 5f6.2)'

      finalfibertag[ifiber] = fibertag[indx]
      ; Identify all objects within 2 arcsec of this position, and

      adist = djs_diff_angle(fibertag.ra, fibertag.dec, $
          fibertag[indx].ra, fibertag[indx].dec, units='degrees')
      indx = where(adist LT 2./3600. AND strtrim(fibertag.objtype,2) NE 'NA')
    endif

    if (indx[0] NE -1) then begin
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

  ; Create structure to hold poly coeff & fill with zeros  
  medwave2d = medwave # replicate(1, nfiber)
  fitimg = smear_medflux*0.0 + 1.0
  xy2traceset, medwave2d, fitimg, smearset, ncoeff=3, /silent
  smearset.coeff[*] = 0.0

  ; Only do fitting if science & smear have reasonable S/N 
  ok = where(sci_sn GT 3.0 AND smear_sn GT 1.0 $
             AND sci_quality EQ 0 AND smear_quality EQ 0, nok)

  xy2traceset, medwave2d[*,ok], smear_medflux[*,ok], polyset, $
    invvar=smear_mask[*,ok], ncoeff=3, inputfunc=sci_medflux[*,ok], $
    lower = 3, upper = 3
  smearset.coeff[*,ok] = polyset.coeff

  ;---------------------------------------------------------------------------
  ; QA plot
  ;---------------------------------------------------------------------------

  traceset2xy, smearset, wave2d, smear_ratio

  ;!P.MULTI = [0, 1, 2]
  ;for i = 0, nfiber - 1 do begin
  ;  plot, 10.0^finalwave, smooth(finalflux[*,i], 5), xr=[3800, 9200], /xs
  ;  djs_oplot, 10.0^finalwave,  smooth(sciflux[*,i], 5), color='red'
  ;  smsci = smooth(djs_median(sciflux[*,i], width=75), 25)
  ;  djs_oplot, 10.0^finalwave,  smsci * smear_ratio[*,i], color='green'
  ;endfor
  ;!P.MULTI = 0
  
  ptsrc = where(strmatch(finalfibertag.objtype, '*GALAXY*') ne 1 and $
                smearset.coeff[0,*] ne 0.0, nptsrc)

  djs_plot, 10.0^wave2d, smear_ratio[*,ptsrc], xr=[3800, 9200], /xs, $
        yr=[0, 2.5], /ys, xtitle = 'Wavelength', $
        ytitle = 'Smear / Science', $
        title = 'Smear Correction Vectors for Point Sources', /nodata
  for iobj = 0, nptsrc - 1 do $
    djs_oplot, 10.0^finalwave, smear_ratio[*,ptsrc[iobj]], nsum=10
  djs_oplot, [2000, 10000], [1, 1], color='red', thick= 4

  lowsn = where(smearset.coeff[0,*] eq 0.0, nlowsn)
  splog, 'Number of fibers with S/N too low to measure smear correction: ', $
          nlowsn
  
  splog, 'Median smear correction for point sources: ' +  $
          string(djs_median(smear_ratio[*,ptsrc]), format = '(F6.2)')

  return, smearset
  
end
