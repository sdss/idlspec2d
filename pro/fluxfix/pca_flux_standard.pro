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

   fluxfnu = red_kflux * wave^2 / 2.99792e18
   fthru=filter_thru(fluxfnu, waveimg=wave, $
                     filter_prefix = 'sdss_jun2001', /toair)
   model_rmag = -2.5*alog10(fthru[2]) - 48.6
   scalefactor = 10.0^(-0.4 * (fiber.mag[2] - model_rmag))
   red_kflux = red_kflux * scalefactor / 1e-17

   ;-----------
   ; Divide star flux in counts by model flux in erg/s/cm^2/A

   fluxvect = objflux
   fluxvivar = objivar
   divideflat, fluxvect, invvar=fluxvivar,  red_kflux, minval = 0.1
    
   ;-----------
   ; Smooth the flux vector to reduce noise (do we need to smooth invar?)

   fluxvect = djs_median(fluxvect, width = 75, boundary = 'reflect')
   fluxvect = smooth(fluxvect, 25, /NAN)

   return, fluxvect
end

;------------------------------------------------------------------------------

pro pca_flux_standard, loglam, stdflux, stdivar, stdinfo, $
    corvector = corvector, corvivar = corvivar, cormed = cormed, $
    fcor = fcor, fsig = fsig, bkpts = bkpts, calibset = calibset, $
    noplot = noplot

   ;--------------
   ; Read in Kurucz model files

   kurucz_restore, kwave, kflux, kindx = kindx
 
   ;-----------------
   ; Compute the flux correction vector from the ratio of model/data

   cspeed = 2.99792458e5
   npix = n_elements(stdflux[*,0])
   nstd = n_elements(stdflux[0,*])

   corvector = fltarr(npix, nstd)
   corvivar = fltarr(npix, nstd)
   cormed = fltarr(nstd)

   wave = 10.0^loglam 
   medwave = djs_median(wave)
   ; Wavelength range over which to compute normalization
   norm_indx = where(wave gt medwave - 500 and wave lt medwave + 500)

   for istd=0, nstd-1 do begin
     model_index = (where(kindx.model eq stdinfo[istd].model))[0]

     linterp, kwave*(1 + stdinfo[istd].v_off/cspeed), kflux[*,model_index], $
              wave, kfluxi

     corvectori = kfluxratio(wave, stdflux[*,istd], $
                  stdivar[*,istd], kfluxi, stdinfo[istd], $
                  fluxvivar = corvivari)

     ; Normalize by median
     cormed[istd] = djs_median(corvectori[norm_indx])
     corvector[*,istd] = corvectori / cormed[istd]
     corvivar[*, istd] = corvivari * cormed[istd]^2 
   endfor

   ;---------------
   ; Now find the average of the vectors with iterative rejection

   pcaflux = pca_solve(corvector, corvivar, $
             niter=30, nkeep=1, usemask=usemask, eigenval=eigenval, $
             acoeff=acoeff, maxiter=5, upper=5, lower=5, $
             maxrej=ceil(0.01*npix), groupsize=ceil(npix/5.))

   ; Demand that the first eigenspectrum is positive-valued.
   ; (The routine PCA_SOLVE() can return a negative-valued spectrum even
   ; if all the input spectra are positive-valued.)
   if (median(pcaflux[*,0]) LT 0) then begin
      pcaflux[*,0] = -pcaflux[*,0]
      acoeff[0,*] = -acoeff[0,*]
   endif

   fcor = pcaflux * median(acoeff) 

   ;--------------
   ; Measure the variance between the fluxcalib vectors derived 
   ; for each of the standard stars -- this is an indicator of the 
   ; spectrophotometric quality.

   meanclip, corvector - rebin(fcor, npix,  nstd), fmean, fsig, $
             maxiter=3, clipsig=5

   ;--------------
   ; Select break points for spline

   minwave = min(wave)
   maxwave = max(wave)
   if (minwave lt 5000) then bkptfile = filepath('blue.bkpts', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc') $
   else bkptfile = filepath('red.bkpts', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, bkptfile, bkpts, silent=1

   ibk = where(bkpts GE minwave AND bkpts LE maxwave, nbk)
   if (nbk LT 4) then splog, 'Error selecting break points'
   bkpts = [minwave, bkpts[ibk], maxwave]
   bkpts = alog10(bkpts) ; Convert to log-10 Angstroms

   ;--------------
   ; Do the spline fit

   calibset = bspline_iterfit(loglam, fcor, nord=4, bkpt=bkpts, $
              upper=3, lower=3, maxrej=ceil(0.05*n_elements(fcor)))

   calibvector = bspline_valu(loglam, calibset)

   ;----------
   ; QA plot

   if keyword_set(noplot) then return

   djs_plot, [minwave-100,maxwave+100], [0,1.1*max(corvector)], /nodata, $
             /xstyle, /ystyle, xtitle='\lambda [A]', $
             ytitle='Counts / (10^{-17}erg/cm^{2}/s/A)', $
             title = 'Average Spectrophoto Correction'  

   for istd=0, nstd - 1 do oplot, wave, corvector[*,istd] 
   djs_oplot, wave, fcor, color='green', thick=3
   djs_oplot, wave, calibvector, color='red', thick=3
   djs_oplot, 10^bkpts, bspline_valu(bkpts,calibset), psym=4, color='red'
 
   xyouts, mean(wave) - 500, [0.9*max(corvector)], $
     'Standard star variation = ' + string(fsig * 100, format='(I3)') + ' %' 

   return

end
;------------------------------------------------------------------------------
