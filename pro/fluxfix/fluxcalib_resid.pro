;+
; NAME:
;  fluxcalib_resid
;
; PURPOSE:
;   Solve for flux-calibration vectors on a single camera + spectrograph.
;   This is accomplished by matching the normalized standard star spectra 
;   to Kurucz models and then ratioing the spectrum (in counts) to the best 
;   fitting model (in units of ergs/s/cm^2/A) scaled to match the photo fiber 
;   mag.  The final flux correction vector is made from the average of the 
;   vectors for each of the standard stars after rejecting outliers.
;
; CALLING SEQUENCE:
;    flux_standard, loglam, stdflux, stdivar, stdmask, stdstarfile, $
;    outname = , combinedir = , calibset = ,  waveminmax = , corvector = , $
;    goodv = , fcor = ,  fsig = , noplot = 
;
;
; INPUTS:
;   loglam       - wavelength in log10(Angstroms) of the input spectra [npix]
;   stdflux      - array of standard star spectra [npix, nstar]
;   stdivar      - inverse variance of standard star spectra [npix, nstar]
;   stdmask      - and mask of the standard star spectra [npix, nstar]
;   stdstarfile  - name of FITS file containing info about the best model fit 
;                  produced by stype_standard (e.g. spStd-0506-52022-1.fits)
;
; OPTIONAL INPUTS:
;   outname      - name of FITS file to contain the final correction vector
;                  (e.g. fluxcalib-0506-52022-b1.fits) -- also the plot title
;   combinedir   - directory to write the FITS file to
;   waveminmax   - wavelength range to spline fit (in Angstroms) [wlow, whigh]
;   noplot       - toggle plotting
;
; OUTPUT:  FITS file containing the coefficients of a spline fit to the ratio 
;          of the flux in counts to the flux in real units (1e-17 erg/s/cm^2/s).
;          The FITS header keyword 'SPHOTERR' stores a measure of the quality 
;          of the spectrophotometry.  This number is also passed in the keyword
;          "fsig".
;
;          Diagnostic plots are also produced.  The flux correction vector 
;          produced from each high S/N standard is plotted in black; the mean 
;          of all the fluxcor vectors in green; and the bspline of the mean 
;          vector in red.
;
; OPTIONAL OUTPUTS:
;   calibset     - spline coefficients of the fit
;   corvector    - flux calibrations derived from the individual
;                  standard stars [npix, nstds]
;   goodv        - array to indicate if corvector is used in fcor [nstds]
;   fcor         - final flux correction vector [npix] (not spline fit) 
;   fsig         - standard deviation of the flux calibration derived from
;                  individual stars about the final answer (scalar)
;
; COMMENTS:  The file IDLSPEC2D_DIR/etc/kurucz_stds_interp.fit is needed.  This
;            file is described in more detail in "stype_standard".
;
; EXAMPLES:
;
; BUGS: The kurucz models we are using have not been fully tested.  Do they 
;       yield reliable broad band fluxes??
;
; PROCEDURES CALLED:
;   bspline_iterfit()
;   bspline_valu()
;   divideflat
;   djs_filepath()
;   djs_maskinterp()
;   djs_median()
;   djs_oplot
;   djs_plot
;   djs_reject()
;   ext_odonnell()
;   ext_ccm()
;   fibermask_bits()
;   fileandpath()
;   glactc()
;   mrdfits()
;   mwrfits
;   readcol
;   sky mask
;   splog
;   sxaddpar
;   sxpar()
;   traceset2xy
;
;
; INTERNAL SUPPORT ROUTINES
;   kfluxratio()
;
; REVISION HISTORY:
;   25-Nov-2002  Added empirical correction to Kurucz models
;   22-Sep-2002  Revised to use Kurucz models by C. Tremonti
;   08-Sep-2000  Formerly newfluxcalib Written by D. Schlegel & S. Burles
;-
;------------------------------------------------------------------------------
; Compute flux correction vector from ratio of observed star flux to 
; reddened model flux

function kfluxratio, wave, objflux, objivar, kflux, fiber, $
         fluxvivar = fluxvivar

   ;------------
   ; Get extinction from SFD maps
   
   A_v = 3.1 * fiber.e_bv_sfd
   a_odonnell = ext_odonnell(wave, 3.1)
   red_kflux = kflux * exp(-1 * a_odonnell * A_v / 1.086) 

   ;-------------
   ; Get zeropoint from phot fiber mag 
   
   ;fluxfnu = red_kflux * wave^2 / 2.99792e18
   ;fthru=filter_thru(fluxfnu, waveimg=wave, $
   ;                  filter_prefix = 'sdss_jun2001', /toair)
   scalefactor = 10.0^(-0.4 * (fiber.mag[2] - fiber.red_model_mag[2]))

   red_kflux = red_kflux * scalefactor / 1e-17

   ;-----------
   ; Divide star flux in counts by model flux in erg/s/cm^2/A

   fluxvect = objflux
   fluxvivar = objivar
   divideflat, fluxvect, invvar=fluxvivar,  red_kflux, minval = 0.1
   
   ;-----------
   ; Smooth the flux vector to reduce noise (do we need to smooth invar?)

   fluxvect = djs_median(fluxvect, width = 100, boundary = 'reflect')
   fluxvect = smooth(fluxvect, 30, /NAN)

   return, fluxvect
end

;------------------------------------------------------------------------------
; Average the flux correction vectors together weighting by the S/N
; of the spectrum (since this influences the quality of the spectral typing). 
; First average the high frequency information, iteratively rejecting 
; outliers at each pixel, then combine the low frequency info, rejecting 
; whole fibers where necessary

function avg_fluxvect, wave, corvector, sn, cvectok = cvectok, noplot = noplot

   npix = n_elements(corvector[*,0])
   nok = n_elements(corvector[0,*])

   ;-----------------
   ; Fit out low order terms

   tracewave = where(wave gt 4400, ntw)
   wave2d = rebin(wave, npix, nok)
   xy2traceset, wave2d[tracewave, *], corvector[tracewave, *], polyset, $
      ncoeff=4, lower = 3, upper = 3
   traceset2xy, polyset, wave2d, corvect_lowf 

   corvect_hif = corvector / corvect_lowf

   ;-----------------
   ; Call iterative rejection on hi-frequency correction vectors

   iiter = 0
   maxiter = 3 
   outmask_hif = fltarr(npix, nok) + 1
   weight_hif = transpose(rebin(sn, nok, npix))

   fcor_hif = djs_median(corvect_hif, 2)

   while (NOT keyword_set(qdone) AND iiter LE maxiter) do begin
     qdone = djs_reject(corvect_hif, rebin(fcor_hif, npix,  nok), $
             outmask=outmask_hif, maxrej = 3 < (nok - 2), /sticky, $
             upper=3, lower=3, groupdim = 1)

     ; Toss out fibers where more than 1/3 of the pixels get rejected
     badfiber = where(total(outmask_hif, 1) lt npix * 2.0/3.0)
     if badfiber[0] ne -1 then outmask_hif[*,badfiber] = 0

     fcor_hif = total(corvect_hif * outmask_hif * weight_hif, 2) / $
                total(weight_hif * outmask_hif, 2)
     fcor_hif = djs_median(fcor_hif, width = 50, boundary = 'reflect')

     iiter = iiter + 1
   endwhile

   ok = where(total(outmask_hif, 1) gt 0)
   nok = n_elements(ok)
   cvectok = ok

   ;-----------------
   ; Now find the median of the low-frequency correction vectors
   
   lowfset = {func: 'legendre', xmin: polyset.xmin, xmax: polyset.xmax, $
             coeff: djs_median(polyset.coeff, 2)}
   traceset2xy, lowfset, wave, fcor_lowf

   ;---------------- 
   ; Put low and high pieces back together
 
   fcor = fcor_hif * fcor_lowf

   ;---------------
   ; Clip out any extreme outliers (usually @ blue edge)

   bad = where(fcor lt 0.4 or fcor gt 1.6)
   if bad[0] ne -1 then fcor[bad] = 1

   ;---------------------
   ; QA plot

   if not keyword_set(noplot) then begin
     plot, indgen(5), xr = [3800, 9200], /xs, yr = [0.8, 1.2], $
       xtitle = 'Wavelength', ytitle = 'Data / Model', $
       title = 'Rectified Spectrophotometry Residuals'
     for i = 0, nok - 1 do oplot, wave, corvect_hif[*,i]
     djs_oplot, wave, fcor_hif, color='green', thick=3
   endif

   return, fcor
end

;------------------------------------------------------------------------------
pro fluxcalib_resid, loglam, stdflux, stdivar, stdmask, stdstarfile, $
    outname = outname, combinedir = combinedir, calibset = calibset, $
    waveminmax = waveminmax, corvector = corvector, goodv = goodv, $
    fcor = fcor, fsig = fsig, noplot = noplot, title_tag = title_tag

   if not keyword_set(title_tag) then title_tag = ''

   ;----------
   ; Read info about f-star standards 
   stdstars = mrdfits(stdstarfile, 1)

   nstd = n_elements(stdstars)
   if (nstd EQ 0) then $
     splog, 'No SPECTROPHOTO or REDDEN stars for flux calibration'

   ;-------------------
   ; Mask out bad pixels and regions dominated by sky-sub residuals

   stdivar = skymask(stdivar, stdmask, ngrow=3)
   stdflux = djs_maskinterp(stdflux, stdivar EQ 0, iaxis=0, /const)

   ;--------------
   ; Read in Kurucz model files

   kurucz_restore, kwave, kflux, kindx = kindx

   ;-------------
   ; Compute the offset of the models from the data in pixels 

   dloglam = 1e-4
   model_offset = (loglam[0] - alog10(kwave[0])) / dloglam
   v_pix = alog10(stdstars.v_off / 3e5 + 1) / dloglam
   pixshift = round(model_offset - v_pix) 

   ;-----------------
   ; Compute the flux correction vector from the ratio of model/data

   npix = (size(stdflux,/dimens))[0]
   corvector = fltarr(npix, nstd)
   corvivar = fltarr(npix, nstd)
   wave = 10.0^loglam 

   for iobj=0, nstd-1 do begin
     model_index = (where(kindx.model eq stdstars[iobj].model))[0]
     pix_index = lindgen(n_elements(kwave)) + pixshift[iobj]

     corvector[*,iobj] = kfluxratio(wave, stdflux[*,iobj], stdivar[*,iobj], $
                         kflux[pix_index, model_index], stdstars[iobj], $
                         fluxvivar = fluxvivar)
     corvivar[*, iobj] = fluxvivar

   endfor

   ;---------------------
   ; Chose high S/N spectra (otherwise spectral typing is likely to be 
   ; wrong).  The velocity cut (hopefully) eliminates galaxies and 
   ; other things that accidentally get targeted
   
   ok = where(stdstars.sn gt 20 and abs(stdstars.v_off) lt 450, nok)

   if nok lt 3 then begin
     splog, 'WARNING:  Too few spectrophoto standards with good S/N'
     splog, 'Proceeding anyway!'
     ok = where(stdstars.sn gt 10 and abs(stdstars.v_off) lt 450, nok)
   endif
   if nok lt 3 then begin
     splog, 'WARNING:  NO good spectrophoto standards!!!'
     splog, 'Proceeding anyway!'
     ok = where(stdstars.sn gt 0, nok)
   endif else begin
     splog, 'Using ' + string(nok, format = '(I2)') + ' spectrophoto standards'
   endelse  

   fcor = fltarr(npix)  
   corsig = fltarr(npix)

   ;----------------------
   ; Trim to range of good pixels
 
   if keyword_set(waveminmax) then begin 
     goodpix = where(wave GT waveminmax[0] AND wave LT waveminmax[1]) 
     if goodpix[0] EQ -1 then $
       splog, 'Bad wavelength range for spectrophotometric calibration'
   endif else goodpix = lindgen(npix)

   npix = n_elements(goodpix)
   wave = wave[goodpix,*]
   corvector = corvector[goodpix,*]
   corvivar = corvivar[goodpix,*]

   ;----------------------
   ; Normalize the flux correction vectors and set to the mean
  
   ;normwave = where(wave gt 5491 and wave lt 6873) ; r-band
   normwave = where(wave gt 5300 and wave lt 5700) ; center of guiding???

   corrmed = fltarr(nstd)
   for i = 0, nstd - 1 do corrmed[i] = median(corvector[normwave, i])
   meanclip, corrmed[ok], corrmean
   corvector = corvector / (corrmed ## replicate(1, npix)) * corrmean
 
   ;---------------
   ; Now find the average vector, with iterative rejection 

   fcor = avg_fluxvect(wave, corvector[*,ok], stdstars[ok].sn, $
                       cvectok = cvectok)

   ;--------------
   ; Recompute the mean without rejected fibers

   ok = ok[cvectok]
   nok = n_elements(ok)
   goodv = intarr(nstd)
   goodv[ok] = 1
   meanclip, corrmed[ok], newcorrmean

   corvector = corvector / corrmean * newcorrmean
   fcor = fcor / corrmean * newcorrmean

   ;----------
   ; Now measure the variance between the fluxcalib vectors derived 
   ; for each of the standard stars -- this is an indicator of the 
   ; spectrophotometric quality.

   meanclip, corvector[*,ok] - rebin(fcor, npix,  nok), fmean, fsig, $
             maxiter=3, clipsig=5

   ;----------
   ; Select break points for spline

   logmin = min(loglam[goodpix])
   logmax = max(loglam[goodpix])

   bbkptfile = filepath('blue.bkpts', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc') 
   readcol, bbkptfile, bbkpts, silent=1
   rbkptfile = filepath('red.bkpts', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, rbkptfile, rbkpts, silent=1
   bkpts = alog10([bbkpts, rbkpts]) ; Convert to log-10 Angstroms

   ibk = where(bkpts GE logmin AND bkpts LE logmax, nbk)
   if (nbk LT 4) then splog, 'Error selecting break points'
   bkpts = [logmin, bkpts[ibk], logmax]

   ;----------
   ; Do the spline fit

   calibset = bspline_iterfit(loglam[goodpix], fcor, nord=4, bkpt=bkpts, $
              upper=3, lower=3, maxrej=ceil(0.05*n_elements(fcor)))

   calibvector = bspline_valu(loglam[goodpix],calibset)

   if keyword_set(outname) then begin
     ;----------
     ; Create header cards describing the data range and write to FITS
     ; The keyword 'SPHOTERR' holds the standard deviation of the
     ; correction vectors for the individual stars -- this is a good measure
     ; of the quality of the spectrophotometry

     hdr = ['']
     wavemin = 10.^min(loglam)
     wavemax = 10.^max(loglam)
     sxaddpar, hdr, 'WAVEMIN', wavemin
     sxaddpar, hdr, 'WAVEMAX', wavemax
     sxaddpar, hdr, 'SPHOTERR', fsig  

     fluxcalibfile = djs_filepath(outname, root_dir=combinedir)
     mwrfits, calibset, fluxcalibfile, hdr, /create
   endif else outname = ''

   ;----------
   ; QA plot

   if keyword_set(noplot) then return

   djs_plot, [min(wave)-100,max(wave)+100], [0.6,1.4], /nodata, $
             /xstyle, /ystyle, xtitle='Wavelength', $
             ytitle='Data / Model', $
             title = 'Spectrophotometry Residuals: ' + title_tag

   for i=0, nok - 1 do oplot, wave, corvector[*,ok[i]]
   djs_oplot, wave, fcor, color='green', thick=3
   djs_oplot, wave, calibvector, color='red', thick=3
   djs_oplot, 10^bkpts, bspline_valu(bkpts,calibset), psym=4, color='red'
 
   xyouts, min(wave) + 500, 1.3, $
     'Standard star variation = ' + string(fsig * 100, format = '(I3)') + ' %' 

   return

end
;------------------------------------------------------------------------------
