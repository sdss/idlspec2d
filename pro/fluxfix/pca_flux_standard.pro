;+
; NAME:
;  flux_standard
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
;    outname = , calibset = ,  waveminmax = , corvector = , $
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
   ; (assume reddening does not change the SED too much)

   ebv_to_A_r = 2.6980 ; derived for F-star SED
   model_rmag = fiber.model_mag[2] + (ebv_to_A_r * fiber.e_bv_sfd)
   scalefactor = 10.0^(-0.4 * (fiber.mag[2] - model_rmag))

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
pro pca_flux_standard, loglam, stdflux, stdivar, stdmask, stdstarfile, $
    outname = outname, calibset = calibset, blue = blue, $
    waveminmax = waveminmax, corvector = corvector, goodv = goodv, $
    fcor = fcor, fsig = fsig, noplot = noplot

   cspeed = 2.99792458e5

   ;----------
   ; Read info about f-star standards 
   stdstars = mrdfits(stdstarfile, 1)

   nstd = n_elements(stdstars)
   if (nstd EQ 0) then $
     splog, 'No SPECTROPHOTO or REDDEN stars for flux calibration'

   ;-------------------
   ; Mask out bad pixels and regions dominated by sky-sub residuals

   stdivar = skymask(stdivar, stdmask)
   stdflux = djs_maskinterp(stdflux, stdivar EQ 0, iaxis=0, /const)

   ;--------------
   ; Read in Kurucz model files

   kurucz_restore, kwave, kflux, kindx = kindx

   ;----------------------
   ; Trim to range of good pixels
 
   wave = 10.0^loglam 
   if not keyword_set(waveminmax) then waveminmax = [min(wave), max(wave)]
   goodpix = where(wave GE waveminmax[0] AND wave LE waveminmax[1]) 
   if goodpix[0] EQ -1 then $
       splog, 'Bad wavelength range for spectrophotometric calibration'

   npix = n_elements(goodpix)
   wave = wave[goodpix]

   ;-----------------
   ; Compute the flux correction vector from the ratio of model/data

   corvector = fltarr(npix, nstd)
   corvivar = fltarr(npix, nstd)

   for iobj=0, nstd-1 do begin
     model_index = (where(kindx.model eq stdstars[iobj].model))[0]

     linterp, kwave*(1 + stdstars[iobj].v_off/cspeed), kflux[*,model_index], $
              wave, kfluxi

     corvector[*,iobj] = kfluxratio(wave, stdflux[goodpix,iobj], $
                         stdivar[goodpix,iobj], kfluxi, stdstars[iobj], $
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

   ;---------------
   ; Now find the average vector, with iterative rejection 

   pcaflux = pca_solve(corvector[*,ok], corvivar[*,ok], $
             niter=30, nkeep=1, usemask=usemask, eigenval=eigenval, $
             acoeff=acoeff, maxiter=5, upper=5, lower=5, $
             maxrej=ceil(0.01*npix), groupsize=ceil(npix/5.))

   ; Demand that the first eigenspectrum is positive-valued.
   ; (The routine PCA_SOLVE() can return a negative-valued spectrum even
   ; if all the input spectra are positive-valued.)
   if (median(pcaflux[*,0]) LT 0) then begin
      pcaflux[*,0] = -pcaflux[*,0]
      acoeff[*,0] = -acoeff[*,0]
   endif

   fcor = pcaflux / n_elements(ok)

   ;--------------
   ; Measure the variance between the fluxcalib vectors derived 
   ; for each of the standard stars -- this is an indicator of the 
   ; spectrophotometric quality.

   meanclip, corvector[*,ok] - rebin(fcor, npix,  nok), fmean, fsig, $
             maxiter=3, clipsig=5
   fsig = fsig / median(fcor)

   ;--------------
   ; Select break points for spline

   logmin = min(loglam[goodpix])
   logmax = max(loglam[goodpix])

   if (min(wave) lt 5000) then bkptfile = filepath('blue.bkpts', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc') $
   else bkptfile = filepath('red.bkpts', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, bkptfile, bkpts, silent=1
   bkpts = alog10(bkpts) ; Convert to log-10 Angstroms

   ibk = where(bkpts GE logmin AND bkpts LE logmax, nbk)
   if (nbk LT 4) then splog, 'Error selecting break points'
   ;bkpts = [logmin, bkpts[ibk[1:nbk-2]], logmax]
   bkpts = [logmin, bkpts[ibk], logmax]

   ;--------------
   ; Do the spline fit

   calibset = bspline_iterfit(loglam[goodpix], fcor, nord=4, bkpt=bkpts, $
              upper=3, lower=3, maxrej=ceil(0.05*n_elements(fcor)))

   calibvector = bspline_valu(loglam[goodpix],calibset)

   ;--------------
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

   if keyword_set(outname) then $
     mwrfits, calibset, outname, hdr, /create $
   else outname = ''

   ;----------
   ; QA plot

   if keyword_set(noplot) then return

   djs_plot, [min(wave)-100,max(wave)+100], [0,1.1*max(corvector)], /nodata, $
             /xstyle, /ystyle, xtitle='\lambda [A]', $
             ytitle='Counts / (10^{-17}erg/cm^{2}/s/A)', $
             title = ' Spectrophoto Correction: ' + outname

   for i=0, nok - 1 do oplot, wave, corvector[*,ok[i]]
   djs_oplot, wave, fcor, color='green', thick=3
   djs_oplot, wave, calibvector, color='red', thick=3
   djs_oplot, 10^bkpts, bspline_valu(bkpts,calibset), psym=4, color='red'
 
   xyouts, (max(wave) - min(wave))/2 + min(wave) - 500, 0.9*max(calibvector), $
     'Standard star variation = ' + string(fsig * 100, format = '(I3)') + ' %' 

   return

end
;------------------------------------------------------------------------------
