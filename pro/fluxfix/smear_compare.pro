;+
; NAME:
;   smear_compare
;
; PURPOSE:
;   Ratio calibrated smear and science spectra and fit with a 3rd order 
;   polynomial.
;
; CALLING SEQUENCE:
;   smearset = smear_compare(smearname, sciloglam, sciflux, sciivar, $
;     best_exp_str, plate_str, mjd_str, [camnames=, combinedir=, $
;     tsobjname=, nord=nord, adderr=, noplot=])
;
; INPUTS:
;   smearname - list of smear exposure names [nsmear * nexp]
;   sciloglam - wavelength (log10 (A)) of calibrated science image [npix]
;   sciflux   - flux of fully calibrated science image [npix,nfiber]
;   sciivar   - inverse variance of calibrated science image [npix,nfiber]
;   best_exp_str - string containing the exposure ID of the best science image.
;                  The fluxcalib vector from this image takes the place of the
;                  PCA average vector for the smear.
;   plate_str - plate id -- used to reconstruct standard star FITS filename
;   mjd_str - MJD -- used to reconstruct standard star FITS filename
;
; OUTPUT:
;   A structure containing the median smear and science S/N and the 
;   Legendere coeff of the smear/sci flux ratio for each fiber is returned
;
; KEYWORDS:
;   camnames   - defaults to ['b1', 'b2', 'r1', 'r2']
;   combinedir - optional output directory
;   tsobjname  - full path name of tsobjfile 
;   nord       - order of spline fit when combining exposures (defaults to 3)
;   adderr     - additional error to add to formal error 
;   noplot     - toggle some of the plotting
;
; COMMENTS:
;   The smear images are calibrated in a manner nearly identical to the 
;   calibration of the science images.  The smear and science image are
;   coarsely binned with locally deviant points rejected.  The ratio 
;   (smear/sci) is fit with a 3rd order Legendre polynomial.  The structure
;   which is returned contains the Legendre coeff as well as the median S/N
;   of the binned smear and science frames.  The smear/sci ratio is probably
;   meaningless for object with a S/N < 3 in the science image or < 1 in the
;   smear.
;   
; BUGS:
;   The smear/science flux ratio is only meaningful for high S/N point 
;   sources!!!!
;   
; EXAMPLES:
;   To turn the coeff back into vectors do:
;  
;   smear_hdu = mrdfits('spPlate-0503-51999.fits', 8) 
;   smearset = {func: 'legendre', xmin: alog10(3800.0), xmax: alog10(9200.0), $
;               coeff: dblarr(3,640)}
;   smearset.coeff = smear_hdu.legendre_coeff
;   readspec, 503, mjd=51999, wave=wave
;   traceset2xy, smearset, alog10(wave), sratio
;   ok = where(smear_hdu.smear_sn gt 1 and smear_hdu.sci_sn gt 3)
;   plot, wave[*,ok], sratio[*,ok], psym=3
;
; PROCEDURES CALLED:
;   bspline_valu
;   combine1fiber
;   divideflat
;   djs_filepath
;   djs_icolor
;   djs_oplot
;   djs_plot
;   legend
;   sphoto_calib
;   splog
;   spmedian_rebin
;   spread_frames
;   traceset2xy
;   xy2traceset
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   12-Aug-2003  Created by C. Tremonti, Steward Observatory
;-
;------------------------------------------------------------------------------

function smear_compare, smearname, sciloglam, sciflux, sciivar, $
         best_exp_str, plate_str, mjd_str, noplot = noplot, $
         camnames = camnames, combinedir = combinedir, $
         tsobjname = tsobjname, nord=nord, adderr=adderr 

  ;------------------------------------
  ; Read frame files into big 2-d array

  spread_frames, smearname, window=window, binsz = binsz, $
     adderr=adderr, camnames=camnames, tsobjname = tsobjname, $
     flux = flux, ivar = fluxivar, loglam = loglam, pixelmask = pixelmask, $
     plugtag = plugtag, camerasvec = camerasvec, expid = expid

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

     sphoto_calib, loglam[*,isphoto], flux[*,isphoto], fluxivar[*,isphoto], $
       pixelmask[*,isphoto], plugtag[isphoto], fcalfiles, $
       stdstarfile, input_calibset=incalib, noplot=noplot

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
       calibfac = bspline_valu(loglam[*,indx], calibset)

       ; Set to bad any pixels whose wavelength is outside the known
       ; flux-calibration region.
       ibad = where(loglam[*,indx] LT alog10(cwavemin) OR $
                    loglam[*,indx] GT alog10(cwavemax))
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

  nfinalpix = n_elements(sciloglam)
  nfiber = max(plugtag.fiberid)
  finalflux = fltarr(nfinalpix, nfiber)
  finalivar = fltarr(nfinalpix, nfiber)
  finalandmask = lonarr(nfinalpix, nfiber)
  finalplugtag = replicate(plugtag[0], nfiber)
   struct_assign, {fiberid: 0L}, finalplugtag
  binsz = sciloglam[1] - sciloglam[0]

  splog, 'Combining smear exposures'

  for ifiber=0, nfiber-1 do begin

    ; Find the first occurrence of fiber number IFIBER+1
    indx = where(plugtag.fiberid EQ ifiber+1)

    if (indx[0] NE -1) then begin
;     splog, 'Fiber', ifiber+1, ' ', plugtag[indx[0]].objtype, $
;     plugtag[indx[0]].mag, format = '(a, i4.3, a, a, f6.2, 5f6.2)'

      finalplugtag[ifiber] = plugtag[indx[0]]
      temppixmask = pixelmask[*,indx]

      combine1fiber, loglam[*,indx], flux[*,indx], fluxivar[*,indx], $
          finalmask=temppixmask, andmask=bestandmask, $
          newloglam=sciloglam, newflux=bestflux, newivar=bestivar, $
          nord=nord, binsz=binsz, maxiter=50, upper=3.0, lower=3.0, maxrej=1

      finalflux[*,ifiber] = bestflux
      finalivar[*,ifiber] = bestivar
      finalandmask[*,ifiber] = bestandmask

    endif else begin
      splog, 'Fiber', ifiber+1, ' NO DATA'
      finalandmask[*,ifiber] = pixelmask_bits('NODATA')
    endelse
  endfor 

  ; Free memory
  loglam = 0
  flux = 0
  fluxivar = 0

  ;---------------------------------------------------------------------------
  ; Compare smear & final science image
  ;---------------------------------------------------------------------------

  loglam2d = sciloglam # replicate(1, nfiber)

  smear_medflux = spmedian_rebin(loglam2d, finalflux, finalivar, 'full', $
                  outwave = medwave, sn = smear_sn, mask = smear_mask, $
                  quality = smear_quality) 
 
  sci_medflux = spmedian_rebin(loglam2d, sciflux, sciivar, 'full', $
                sn = sci_sn, mask = sci_mask, quality = sci_quality) 

  ; Fit with 4th order Legendre polynomial  
  medwave2d = float(medwave) # replicate(1, nfiber)

  xy2traceset, medwave2d, smear_medflux, smearset, invvar=smear_mask, $
    ncoeff=3, inputfunc=sci_medflux, lower=3, upper=3

  ; Save as binary table with science & smear S/N ratios
  smear_struct = {sci_sn: 0.0, smear_sn: 0.0, legendre_coeff: fltarr(3)}

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
  traceset2xy, smearset, loglam2d, smear_ratio

  ;!P.MULTI = [0, 1, 2]
  ;for i = 0, nfiber - 1 do begin
  ;  plot, 10.0^sciloglam, smooth(finalflux[*,i], 5), xr=[3800, 9200], /xs
  ;  djs_oplot, 10.0^sciloglam,  smooth(sciflux[*,i], 5), color='red'
  ;  smsci = smooth(djs_median(sciflux[*,i], width=75), 25)
  ;  djs_oplot, 10.0^sciloglam,  smsci * smear_ratio[*,i], color='green'
  ;endfor
  ;!P.MULTI = 0
  
  ptsrc = where(strmatch(finalplugtag.objtype, '*GALAXY*') ne 1 and ok, nptsrc)
  qso = where(strmatch(finalplugtag.objtype, '*QSO*') eq 1 and ok, nqso)

  wave=10.0^sciloglam
  if not keyword_set(noplot) then begin
    djs_plot, wave, smear_ratio[*,0], xr=[3800, 9200], /xs, $
      yr=[0, 2.5], /ys, xtitle = 'Wavelength', ytitle = 'Smear / Science', $
      title = 'Smear Correction Vectors for Point Sources', /nodata
    for iobj = 0, nptsrc - 1 do $
      djs_oplot, wave, smear_ratio[*,ptsrc[iobj]], nsum=10, color='green'
    for iobj = 0, nqso - 1 do $
      djs_oplot, wave, smear_ratio[*,qso[iobj]], nsum=10, color='blue'
    djs_oplot, [2000, 10000], [1, 1], thick= 4
    legend, ['Quasars', 'Other Point Sources'], $
            color=djs_icolor(['blue', 'green']), psym=[0,0]
  endif

  splog, 'Median smear correction for point sources: ' +  $
          string(djs_median(smear_ratio[*,ptsrc]), format = '(F6.2)')

  return, smear_hdu
  
end
