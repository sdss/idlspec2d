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

pro pca_flux_standard, loglam, stdflux, stdivar, stdinfo, camid, $
    corvector = corvector, corvivar = corvivar, cormed = cormed, $
    fcor = fcor, fsig = fsig, bkpts = bkpts, calibset = calibset, $
    noplot = noplot 

   ; Compute ratio of data to model for each standard
   corvector = spdata2model_ratio(loglam, stdflux, stdivar, stdmask, stdinfo, $
               corvivar = corvivar)

   nstd = n_elements(corvector[0,*])
   npix = n_elements(corvector[*,0])
   cormed = fltarr(nstd)

   ; Normalize in the dichroic region but avoiding the exact edges
   wave = 10.0^loglam 
   norm_indx = where(wave gt 5700 and wave lt 6300 and $
                     wave lt max(wave) - 200 and wave gt min(wave) + 200)

   for istd=0, nstd-1 do begin
     cormed[istd] = djs_median(corvector[norm_indx,istd])
     corvector[*,istd] = corvector[*,istd] / cormed[istd]
     corvivar[*, istd] = corvivar[*,istd] * cormed[istd]^2 
   endfor

   ;---------------
   ; Now find the average of the vectors with iterative rejection

   pcaflux = pca_solve(corvector, corvivar, $
             niter=30, nkeep=1, usemask=usemask, eigenval=eigenval, $
             acoeff=acoeff, maxiter=5, upper=5, lower=5, $
             maxrej=ceil(0.01*npix), groupsize=ceil(npix/5.), /quiet)

   ; Demand that the first eigenspectrum is positive-valued.
   ; (The routine PCA_SOLVE() can return a negative-valued spectrum even
   ; if all the input spectra are positive-valued.)
   if (median(pcaflux[*,0]) LT 0) then begin
      pcaflux[*,0] = -pcaflux[*,0]
      acoeff[0,*] = -acoeff[0,*]
   endif

   ;fcor = pcaflux * median(acoeff) 
   fcor = pcaflux / djs_median(pcaflux[norm_indx])

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
             title = 'Average Spectrophoto Correction for ' + camid + ' Frames' 

   for istd=0, nstd - 1 do oplot, wave, corvector[*,istd] 
   djs_oplot, wave, fcor, color='green', thick=3
   djs_oplot, wave, calibvector, color='red', thick=3
   djs_oplot, 10^bkpts, bspline_valu(bkpts,calibset), psym=4, color='red'
 
   ;xyouts, mean(wave) - 500, [0.9*max(corvector)], $
   ;  'Standard star variation = ' + string(fsig * 100, format='(I3)') + ' %' 

   return

end
;------------------------------------------------------------------------------
