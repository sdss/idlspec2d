;+
; NAME:
;  pca_flux_standard
;
; PURPOSE:
;   Solve for the mean flux-calibration vector of one camera + spectrograph
;   using data from all available exposures. This is accomplished by matching 
;   the normalized standard star spectra to Kurucz models and then ratioing 
;   the spectrum (in counts) to the best fitting model (in units of 
;   ergs/s/cm^2/A) scaled to match the photo fiber mag.  Because the individual
;   exposures were taken under different conditions, the zeropoint of the
;   flux calibration vectors from different exposures will be different.  To
;   eliminate this problem all of the vectors are normalized between
;   5700-6300 A. (The red/blue cameras sample different parts of this
;   wavelength region.)  The average flux correction vector is the first 
;   PCA eigenspectrum.
;
; CALLING SEQUENCE:
;   pca_flux_standard, loglam, stdflux, stdivar, stdmask, stdinfo, [camid=, $
;     corvector=, corvivar=, cormed=, fcor=, fsig=, noplot=] 
;
; INPUTS:
;   loglam   - wavelength in log10(Angstroms) of the input spectra [npix]
;   stdflux  - array of standard star spectra [npix, nstar*nexp]
;   stdivar  - inverse variance of standard star spectra [npix, nstar*nexp]
;   stdmask  - or mask of the standard star spectra [npix, nstar*exp]
;   stdinfo  - structure containing info about the best model fit produced by 
;              "stype_standard"
;
; OPTIONAL INPUTS:
;   camid        - string giving the camera id for plot titles -- i.e. 'b1'
;   noplot       - toggle plotting
;
; OUTPUT:  
;   A structure containing the coefficients of a spline fit the average
;   flux correction vector (normalized at 5700-6300 A).
;
;   Diagnostic plots are also produced.  The flux correction vector produced
;   from each high S/N standard is plotted in black; the mean of all the 
;   fluxcor vectors in green; and the bspline of the mean vector in red.
;
; OPTIONAL OUTPUTS:
;   corvector    - flux calibrations derived from the individual standard
;                  stars [npix, nstar*nexp]
;   corvivar     - inverse variance of flux calib vectors [npix, nstar*nexp]
;   fcor         - final flux correction vector [npix] (not spline fit) 
;   fsig         - standard deviation of the flux calibration derived from
;                  individual stars about the final answer (scalar)
;
; COMMENTS:  
;
; EXAMPLES:
;
; BUGS: 
;   The Kurucz models we are using have not been fully tested.  Do they 
;   yield reliable broad band fluxes at all T_eff & [Fe/H]??
;
;   The mean flux calib vector produced by the PCA sometimes shows abrupt
;   jumps in flux -- this is somewhat mitigated by the spline fit.
;
; PROCEDURES CALLED:
;   bspline_bkpts()
;   bspline_iterfit()
;   bspline_valu()
;   djs_iterstat()
;   djs_median()
;   pca_solve
;   skymask()
;   spdata2model_ratio
;
; INTERNAL SUPPORT ROUTINES
;
; REVISION HISTORY:
;   12-Aug-2003  Split off corvector creation into "spdata2model"
;   22-Sep-2002  Revised to use Kurucz models by C. Tremonti, JHU
;   08-Sep-2000  Formerly newfluxcalib Written by D. Schlegel & S. Burles
;-
;------------------------------------------------------------------------------

function pca_flux_standard, loglam, stdflux, stdivar, stdmask, stdinfo, $
    camid, corvector = corvector, corvivar = corvivar, cormed = cormed, $
    fcor = fcor, fsig = fsig, noplot = noplot 

   ; Compute ratio of data to model for each standard
   corvector = spdata2model_ratio(loglam, stdflux, stdivar, stdmask, stdinfo, $
               corvivar = corvivar, cormed=cormed, /norm)

   nstd = n_elements(corvector[0,*])
   npix = n_elements(corvector[*,0])

   ; chose same wavelengths for normalization as used for corvectors
   wave = 10.0^loglam 
   norm_indx = where(wave gt 5700 and wave lt 6300 and $
                     wave lt max(wave) - 200 and wave gt min(wave) + 200)

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

   djs_iterstat, corvector - rebin(fcor, npix,  nstd), mean=fmean, $
     sigma=fsig, maxiter=3, sigrej=5

   ;--------------
   ; Do the spline fit

   padbkpts = bspline_bkpts(loglam, nord=4, nbkpts=50, bkpt=bkpts, /silent)

   calibset = bspline_iterfit(loglam, fcor, nord=4, bkpt=bkpts, $
              upper=3, lower=3, maxrej=ceil(0.05*n_elements(fcor)))

   calibvector = bspline_valu(loglam, calibset)

   ;----------
   ; QA plot

   if keyword_set(noplot) then return, calibset
   if not keyword_set(camid) then camid = ''

   djs_plot, [min(wave)-100,max(wave)+100], [0,1.1*max(corvector)], /nodata, $
             /xstyle, /ystyle, xtitle='\lambda [A]', $
             ytitle='Counts / (10^{-17}erg/cm^{2}/s/A)', $
             title = 'Average Spectrophoto Correction for ' + camid + ' Frames' 

   for istd=0, nstd - 1 do oplot, wave, corvector[*,istd] 
   djs_oplot, wave, fcor, color='green', thick=3
   djs_oplot, wave, calibvector, color='red', thick=3
   djs_oplot, 10^bkpts, bspline_valu(bkpts,calibset), psym=4, color='red'
 
   ;xyouts, mean(wave) - 500, [0.9*max(corvector)], $
   ;  'Standard star variation = ' + string(fsig * 100, format='(I3)') + ' %' 

   return, calibset

end
;------------------------------------------------------------------------------
