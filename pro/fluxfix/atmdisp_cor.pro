;+
; NAME:
;   atmdisp_cor
; 
; PURPOSE:
;   Atmospheric dispersion causes wavelength dependent lightloss from 
;   fiber spectra.  The mean atmospheric dispersion experienced by the 
;   standard stars on each plate is automatically corrected for as part of
;   the flux calibration.  However, any shift in the centering of the 
;   spectroscopic targets in the fibers across the plate will result in 
;   differences in atmospheric dispersion as a function of plate position.
;   These centering errors -- caused by small errors in the plate scale or
;   rotation -- can produce flux calibration errors of up to 20%.  This
;   procedure is designed to correct for these.
;
;   Synthetic magnitudes are computed from the spectra on 1 half-plate and
;   compared to the photo fiber magnitudes.  The residuals are fit as a 
;   function of plate x/y position with a 3rd order polynomial. The (g-r) 
;   offsets of this fit are then mapped into an atmospheric dispersion 
;   correction.  This mapping makes use of models of the atmospheric dispersion
;   of point sources observed with 2" seeing through 3" fibers at a variety
;   of airmasses.  The absolute light loss is mapped using the r-band mag
;   residuals.  A vector of wavelength dependent corrections is returned.
;
; CALLING SEQUENCE:
;
;   atm_model = atmdisp_cor(loglam, flux, plugtag, hdr, title = title, $
;                           surfgr_sig = surfgr_sig)
;
; INPUTS:
;   loglam  - wavelength array of input spectra in log10(Angstroms) [npix]
;   flux    - flux array from 1 half plate [npix, nfiber] (nfiber~320)
;   plugtag - plugmap [nfiber] -- used for mags, plate x/y, sky & standard ID
;   hdr     - image header -- used for aquiring airmass info (maybe also
;             seeing and guiding in the future)
;
; OUTPUT:
;   A model of the atmospheric dispersion correction for each fiber is
;   returned [npix, nfiber].  (To correct: flux_cor = flux / atm_model)
;   Plots are generated showing the (g-r) and r mag offsets as a function 
;   of plate x/y before and after the correction.  Histograms are produced
;   as well.
; 
; KEYWORDS:
;   title     -  title of output plots (usually plate/mjd/specid)
;   surfgr_sig - sigma of 2d fit to the (g-r) offsets -- this is a good
;                indicator of how big the correction is.
;
; COMMENTS:
;   This function requires an external source of info about the effects
;   of atmospheric dispersion as a function of wavelength.  This info is
;   stored in structures in IDLSPEC2D_DIR/etc/.  They are named 
;   'atmdisp_vec_am' + airmass + 'see2.0.fit' where airmass is a value from
;   1.1 - 1.5.
;
; BUGS:
;   Galaxies and stars respond differently to atmospheric dispersion beause
;   of thier different spatial extent.  The mean correction which is 
;   derived here is really that appropriate for the average galaxy -- 
;   so corrections for point sources are likely to be underestimated a bit.
;   To do this right the spatial profile of each object should be taken into
;   account.
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_icolor()
;   djs_iterstat()
;   djs_oplot
;   djs_sfit_iter()
;   filter_thru()
;   gaussfit()
;   legend
;   linterp
;   mrdfits
;   plothist
;   traceset2xy
;   
; INTERNAL SUPPORT ROUTINES:
;   sphoto_xy
;
; REVISION HISTORY:
;   Created 12-Aug-2003 by C. Tremonti, Steward Observatory
;-
;-----------------------------------------------------------------------------

;-----------------------------------------------------------------------------
; Subroutine: sphoto_xy
;-----------------------------------------------------------------------------

; Make plot of spectrophtometry vs plate x or y

pro sphoto_xy, platepos, sphotoff, xx = xx, yy=yy, ytag = ytag, $
    first = first, ok = ok, std = std, rej = rej

  if keyword_set(xx) then xtitle = 'Plate X-position'
  if keyword_set(yy) then xtitle = 'Plate Y-position'
  if keyword_set(first) then ymargin = [2,5] else ymargin = [4,2]

  plot, platepos, sphotoff, psym=3, xtitle = xtitle, $
        ytitle = ytag + ' [Spectro - Photo]', $
        charsize = 1.5, ymargin = ymargin, xr=[-350,350], $
        /xs, yr = [-0.2, 0.4], /nodata

  if keyword_set(ok) then $
    djs_oplot, platepos[ok], sphotoff[ok], psym=6, symsize=0.2, color='blue'

  if keyword_set(std) then $
    djs_oplot, platepos[std], sphotoff[std], psym=6, symsize=0.3, color='red', $
      thick=4

  if keyword_set(rej) then $
    djs_oplot, platepos[rej], sphotoff[rej], psym=2, symsize=0.4, color='green'

  oplot, [-500, 500], [0, 0], thick=3

end

;-----------------------------------------------------------------------------

;-----------------------------------------------------------------------------
; Main pro: atmdisp_cor
;-----------------------------------------------------------------------------

function atmdisp_cor, loglam, flux, plugtag, hdr, title = title, $
  surfgr_sig = surfgr_sig

  npix = n_elements(flux[*,0])
  nfib = n_elements(flux[0,*])
  std = where(strmatch(plugtag.objtype, '*R*STD*') eq 1)
  wave = 10.0^rebin(loglam, npix, nfib)

  ;-----------------------------------------------------------------------------
  ; Read atmospheric dispersion vecotor matched in airmass 
  ; Ignore seeing since not all objects are point sources (seeing mostly 
  ; changes the amplitude of the correction not the shape, so this is ok)
  ;-----------------------------------------------------------------------------

  airmass = sxpar(hdr, 'AIRMASS')
  airmass_fid = ['1.1', '1.2', '1.3', '1.4', '1.5']
  junk = min(abs(airmass - airmass_fid), aindx)

  atmdisp_file = filepath('atmdisp_vec_am' + airmass_fid[aindx] + $
                 'see2.0.fit', $
                 root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')

  dispset = mrdfits(atmdisp_file, 1)
  ndisp = n_elements(dispset.coeff[0,*])
  traceset2xy, dispset, wave[*,0:ndisp-1], dispvec
  mag = -2.5 * alog10(filter_thru(dispvec, wave = wave[*,0], /toair)) 
  groffsynth = mag[*,1] - mag[*,2]  

  ;---------------------------------------------------------------------------
  ; Compute (g-r) offsets of spectro and photo
  ;---------------------------------------------------------------------------

  flam2fnu = (wave*wave / 2.99792e18)
  filt_flux = filter_thru(flux * flam2fnu, wave = wave[*,0], /toair) 
  mag = filt_flux * 0.0
  nonneg = where(filt_flux gt 0)
  mag[nonneg] = -2.5 * alog10(filt_flux[nonneg]) - 48.6 + 2.5*17.0
  smag = transpose(mag[*,[1,2,3]]) 
  pmag = plugtag.mag[[1,2,3]]

  groff = (smag[0,*] - smag[1,*])  - (pmag[0,*] - pmag[1,*])
  rioff = (smag[1,*] - smag[2,*])  - (pmag[1,*] - pmag[2,*])
  xpos = plugtag.xfocal
  ypos = plugtag.yfocal

  ;-----------------------------------------------------------------------------
  ; Fit 2d surface to (g-r) offsets
  ;-----------------------------------------------------------------------------

  ok = where(abs(groff - 0.03) lt 0.5 and pmag[1,*] gt 0 and smag[1,*] gt 0 $
             and strmatch(plugtag.objtype, '*SKY*') ne 1, nok)

  ; Do iterative fit
  mask = groff * 0
  mask[ok] = 1
  acoeff = djs_sfit_iter(groff, xpos, ypos, 3, 3, yfit=zfit, mask=mask, $
           maxrej = 25, maxdev=0.4, lower=2.5, upper=3.0, freeiter=10, $
           maxiter=19, outmask = outmask)

  djs_iterstat, zfit[ok], sigrej=5, sigma=surfgr_sig

  ;-----------------------------------------------------------------------------
  ; Plot (g-r) offsets as a function of x & y before and after
  ;-----------------------------------------------------------------------------

  !P.MULTI = [0, 2, 3]
  rej = where(outmask ne 1)

  ; Before x plot 
  sphoto_xy, xpos, groff, /xx,  ytag = '(g-r)', /first, ok = ok, std = std
  xyouts, 0, 0.32, 'Original', align=0.5, charsize = 1.5, $
          color = djs_icolor('red'), charthick=3
  xyouts, 0.5, 0.97, title, /norm, align=0.5, charsize = 1.5
  legend, ['Standard', 'Rejected'], color=djs_icolor(['red', 'green']), $
          psym=[6, 2], symsize = [0.5, 0.5], charsize = 0.7, thick=3

  ; After x plot   
  sphoto_xy, xpos, groff - zfit, /xx,  ytag = '(g-r)', /first, ok = ok, $
    std = std, rej = rej
  xyouts, 0, 0.32, 'Corrected', align=0.5, charsize = 1.5, $
          color = djs_icolor('red'), charthick=3

  ; Before y plot   
  sphoto_xy, ypos, groff, /yy, ytag = '(g-r)', ok=ok, std = std

  ; After y plot   
  sphoto_xy, ypos, groff - zfit, /yy, ytag = '(g-r)', ok=ok, std = std, rej=rej

  ;-----------------------------------------------------------------------------
  ; Translate (g-r) offsets into correction vectors
  ;-----------------------------------------------------------------------------

  zfit = zfit - median(zfit[std])
  model = flux * 0.0
  ncoeff = n_elements(dispset.coeff[*,0])

  ; Interpolate values of legendre coefficients
  for ii = 0, nfib - 1 do begin
    dispseti = {func: 'legendre', xmin: 3500.0, xmax: 9400.0, $
                coeff: fltarr(ncoeff)}
 
    for jj=0, ncoeff - 1 do begin
      linterp, groffsynth, reform(dispset.coeff[jj,*],ndisp), zfit[ii], coeffi
      dispseti.coeff[jj] = coeffi
    endfor
    traceset2xy, dispseti, wave[*,ii], modeli
    model[*,ii] = modeli
  endfor

  ;-----------------------------------------------------------------------------
  ; Plot (g-r) histogram before/after correction
  ;-----------------------------------------------------------------------------

  ymax = 50
  ok = where(abs(rioff) lt 1 and abs(groff) lt 1 and $
             pmag[1,*] gt 0 and smag[1,*] gt 0)

  plothist, rioff[ok], bin=0.01, xr=[-0.4, 0.4], yr = [0, ymax], charsize=1.5, $
      xtitle = 'color [Spectro - Photo]', ytitle = 'Number of Spectra'
  plothist, rioff[ok], rixhist, riyhist, bin=0.01, color=djs_icolor('red'), $
      /over, thick=4
  plothist, groff[ok], grxhist, gryhist, bin=0.01, color=djs_icolor('green'), $
      /over, thick=4
  oplot, [0, 0], [0, 1e4], thick=4
  grfit = gaussfit(grxhist, gryhist, nterms=3, grcoef)
  rifit = gaussfit(rixhist, riyhist, nterms=3, ricoef)
  xyouts, 0.10, ymax*0.90, '(g-r) offset = ' + $
          string(grcoef[1], format='(F6.3)'), charsize = 0.6
  xyouts, 0.10, ymax*0.80, '(g-r) sigma = ' + $
          string(grcoef[2], format='(F6.3)'), charsize = 0.6
  xyouts, 0.10, ymax*0.70, '(r-i) offset = ' + $
          string(ricoef[1], format='(F6.3)'), charsize = 0.6
  xyouts, 0.10, ymax*0.60, '(r-i) sigma = ' + $
          string(ricoef[2], format='(F6.3)'), charsize = 0.6

  legend, ['(g-r)', '(r-i)'], color=djs_icolor(['green', 'red']), $
          psym=[0,0], thick=4, charsize = 0.6

  ;----------------
  ; correct the flux

  cflux = flux / model
  flam2fnu = (wave*wave / 2.99792e18) 
  filt_flux = filter_thru(cflux * flam2fnu, wave = wave[*,0], /toair) 
  mag = filt_flux * 0.0
  nonneg = where(filt_flux gt 0)
  mag[nonneg] = -2.5 * alog10(filt_flux[nonneg]) - 48.6 + 2.5*17.0
  nsmag = transpose(mag[*,[1,2,3]]) 
      
  ngroff = (nsmag[0,*] - nsmag[1,*])  - (pmag[0,*] - pmag[1,*])
  nrioff = (nsmag[1,*] - nsmag[2,*])  - (pmag[1,*] - pmag[2,*])
  nroff = nsmag[1,*] - pmag[1,*]
  ok = where(abs(nrioff) lt 1 and abs(ngroff) lt 1 and $
             pmag[1,*] gt 0.0 and nsmag[1,*] gt 0.0)

  plothist, nrioff[ok], bin=0.01, xr=[-0.4, 0.4], yr = [0, ymax], $
    charsize=1.5, xtitle = 'color [Spectro - Photo]', $
    ytitle = 'Number of Spectra' 
  plothist, nrioff[ok], rixhist, riyhist, bin=0.01, $
    color=djs_icolor('red'), /over, thick=4
  plothist, ngroff[ok], grxhist, gryhist, bin=0.01, $
    color=djs_icolor('green'), /over, thick=4
  oplot, [0, 0], [0, 1e4], thick=4

  grfit = gaussfit(grxhist, gryhist, nterms=3, grcoef)
  rifit = gaussfit(rixhist, riyhist, nterms=3, ricoef)
  xyouts, 0.10, ymax*0.90, '(g-r) offset = ' + $
          string(grcoef[1], format='(F6.3)'), charsize = 0.6
  xyouts, 0.10, ymax*0.80, '(g-r) sigma = ' + $
          string(grcoef[2], format='(F6.3)'), charsize = 0.6
  xyouts, 0.10, ymax*0.70, '(r-i) offset = ' + $
          string(ricoef[1], format='(F6.3)'), charsize = 0.6
  xyouts, 0.10, ymax*0.60, '(r-i) sigma = ' + $
          string(ricoef[2], format='(F6.3)'), charsize = 0.6

;-----------------------------------------------------------------------------
; Now correct the total fluxes using the r-band mags
;-----------------------------------------------------------------------------

  ;--------------------
  ; Do iterative 2d fit

  ok = where(abs(nroff) lt 1 and pmag[1,*] gt 0 and nsmag[1,*] gt 0 $
             and strmatch(plugtag.objtype, '*SKY*') ne 1, nok)
  mask = nroff * 0
  mask[ok] = 1
  acoeff = djs_sfit_iter(nroff, xpos, ypos, 3, 3, yfit=zfit, mask=mask, $
           maxrej = 25, maxdev=0.4, lower=2.5, upper=3.0, freeiter=10, $
           maxiter=19, outmask = outmask)

  ;--------------------
  ; Compute the zeropoint correction and add it to the model

  ;zfit = zfit - median(zfit[std])
  fluxcor = 10.0^(-0.4*zfit)
  model = model * transpose(rebin(fluxcor, nfib, npix))

  ;------------------------------
  ; Plot r mag offset versus plate x/y before & after the correction

  rej = where(outmask ne 1)

  ; Before x plot 
  sphoto_xy, xpos, nroff, /xx,  ytag = 'r mag', /first, ok=ok, std=std
  xyouts, 0, 0.32, 'Original', align=0.5, charsize = 1.5, $
          color = djs_icolor('red'), charthick=3
  xyouts, 0.5, 0.97, title, /norm, align=0.5, charsize = 1.5
  legend, ['Standard', 'Rejected'], color=djs_icolor(['red', 'green']), $
          psym=[6, 2], symsize = [0.5, 0.5], charsize = 0.7, thick=3

  ; After x plot   
  sphoto_xy, xpos, nroff - zfit, /xx,  ytag = 'r mag', /first, ok=ok, $
    std=std, rej=rej
  xyouts, 0, 0.32, 'Corrected', align=0.5, charsize = 1.5, $
        color = djs_icolor('red'), charthick=3

  ; Before y plot   
  sphoto_xy, ypos, nroff, /yy,  ytag = 'r mag', ok=ok, std=std

  ; After y plot   
  sphoto_xy, ypos, nroff - zfit, /yy,  ytag = 'r mag', ok=ok, std=std, rej=rej

  ;----------------------------------
  ; Plot R-band Histogram before & after

  ymax = 50
  ; before
  plothist, nroff[ok], bin=0.01, xr=[-0.4, 0.4], yr=[0,ymax], charsize=1.5, $
      xtitle = 'r mag [Spectro - Photo]', ytitle = 'Number of Spectra'
  plothist, nroff[ok], rxhist, ryhist, bin=0.01, color=djs_icolor('red'), $
      /over, thick=4
  oplot, [0, 0], [0, 1e4], thick=4

  rfit = gaussfit(rxhist, ryhist, nterms=3, rcoef)
  xyouts, 0.10, ymax*0.90, 'r offset = ' + $
        string(rcoef[1], format='(F6.3)'), charsize = 0.6
  xyouts, 0.10, ymax*0.80, 'r sigma = ' + $
        string(rcoef[2], format='(F6.3)'), charsize = 0.6

  ;----------------
  ; correct the flux

  cflux = flux / model
  flam2fnu = (wave*wave / 2.99792e18) 
  filt_flux = filter_thru(cflux * flam2fnu, wave = wave[*,0], /toair) 
  mag = filt_flux * 0.0
  nonneg = where(filt_flux gt 0)
  mag[nonneg] = -2.5 * alog10(filt_flux[nonneg]) - 48.6 + 2.5*17.0
  nsmag = transpose(mag[*,[1,2,3]]) 
  nroff = nsmag[1,*] -  pmag[1,*]

  ; after
  plothist, nroff[ok], bin=0.01, xr=[-0.4, 0.4], yr = [0, ymax], charsize=1.5, $
      xtitle = 'r mag [Spectro - Photo]', ytitle = 'Number of Spectra'
  plothist, nroff[ok], rxhist, ryhist, bin=0.01, color=djs_icolor('red'), $
      /over, thick=4
  oplot, [0, 0], [0, 1e4], thick=4

  rfit = gaussfit(rxhist, ryhist, nterms=3, rcoef)
  xyouts, 0.10, ymax*0.90, 'r offset = ' + $
  	string(rcoef[1], format='(F6.3)'), charsize = 0.6
  xyouts, 0.10, ymax*0.80, 'r sigma = ' + $
	string(rcoef[2], format='(F6.3)'), charsize = 0.6

  !P.MULTI = 0
  return, model

end
