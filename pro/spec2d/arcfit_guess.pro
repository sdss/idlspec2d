;+
; NAME:
;   arcfit_guess
;
; PURPOSE:
;   Determine initial wavelength solution by comparing spectrum to arc spectrum
;
; CALLING SEQUENCE:
;   wset = arcfit_guess( spec, loglam, intensity, color=color, $
;    [ func=func, bestcorr=bestcorr, acoeff=, dcoeff=, nsteps= ] )
;
; INPUTS:
;   spec       - 1-D spectrum
;   loglam     - Log-lambda of arc lines
;   intensity  - Intensity of arc lines
;
; REQUIRED KEYWORDS:
;   color      - 'red' or 'blue'
;
; OPTIONAL KEYWORDS:
;   func       - Name of fitting function; default to 'legendre'
;   acoeff     - central values of coefficents to explore
;   dcoeff     - range (-.5*dcoeff to .5*dcoeff) of values to explore
;   nsteps     - array of steps to use
;
; OUTPUTS:
;   wset       - traceset (pix -> lambda)
;
; OPTIONAL OUTPUTS:
;   bestcorr   - Correlation coefficient with simulated arc spectrum
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; INTERNAL SUPPORT PROCEDURES:
;   tset_struc()
;   arcfit_iter()
;
; PROCEDURES CALLED:
;   traceset2xy()
;   xy2traceset
;
; REVISION HISTORY:
;   18-Nov-1999  Written by D. Schlegel, Princeton.
;                Excised code from FITARCIMAGE.
;   01-Dec-2000  added acoeff, dcoeff, nsteps keywords
;-
;------------------------------------------------------------------------------

; Define traceset structure
function tset_struc, func, ncoeff, ntrace

   if (NOT keyword_set(ntrace)) then ntrace = 1

   tset = $      
    { func    :    func               , $
      xmin    :    0.0d               , $
      xmax    :    0.0d               , $
      coeff   :    dblarr(ncoeff, ntrace) $
    }

   return, tset
end

;------------------------------------------------------------------------------
; NSTEPS = total number of steps to explore in the range
;   [COEFF-0.5*DCOEFF, COEFF+0.5*DCOEFF]
;   except for the first coefficient, for which NSTEPS is unused.
; NSMOOTH = Sigma in pixels for Gaussian-smoothing kernal
; DLAG    = Step in lags for cross-correlation

function arcfit_iter, spec, loglam, intensity, $
 aset, dcoeff, nsteps, nsmooth=nsmooth, dlag=dlag, bestcorr=bestcorr

   npix = N_elements(spec)
   nline = N_elements(loglam)
   pixarray = 2.0d0 * dindgen(npix) / (npix-1) - 1.0d0

   ncoeff = N_elements(aset.coeff)

   nsz = 4*fix(nsmooth) + 1 ; kernal size (odd)
   gausskern = exp( -( ((nsz-1)/2 - findgen(nsz))^2 ) / nsmooth^2 )
   gausskern = gausskern / total(gausskern)

   npad = fix((npix+1)/2) ; Padding on left and right of spectra before
                          ; cross-correlating

   ; Pad, subtract baseline (with median filter), and smooth input spectrum
   ; Also, take the square-root of the intensity.
   speccorr = spec
   medval1 = median(spec[0:50<npix-1]) ; Median value of left of spectrum
   medval2 = median(spec[npix-50>0:npix-1]) ; Median value of right of spectrum
   speccorr = [fltarr(npad)+medval1, speccorr, fltarr(npad)+medval2] ; Pad
   speccorr = speccorr - median(speccorr, 101<npix) > 1 ; Median-filter
   speccorr = convol(speccorr, gausskern, /center, /edge_truncate) ; Smooth
   speccorr = sqrt(speccorr > 1) ; Weight by the square-root of the intensity

   bestcorr = -2.0
   bestcoeff = 0
   coeff_lo = aset.coeff - dcoeff * (nsteps-1)/(2.*nsteps)
   coeff_lo[0] = aset.coeff[0] ; Set this coefficient to its central value,
                               ; e.g. that passed in ASET.  We will be doing
                               ; the cross-correlation left and right of this.
   tempset = aset


   ; Set minimum number of lags to check equal to 10 pixels
   lagmax = fix( dcoeff[0] / (2. * abs(aset.coeff[1]) / npix) + 1 ) > 10
   nlag = fix(lagmax/dlag)
   lags = (indgen(nlag) - fix(nlag/2)) * dlag
   splog, 'Searching lags from ', min(lags), ' to ', max(lags), $
    ' spaced ', dlag, ' pixels', format='(a,i5,a,i5,a,i3,a)'

   ; Loop over all coefficients except the first one
   nsteptot = 1
   for ic=1, ncoeff-1 do nsteptot = nsteptot * nsteps[ic]
   for istep=0, nsteptot-1 do begin
      ; Set the coefficients COEFFI for this step number
      ntemp = nsteptot
      itemp = istep
      for ic=1, ncoeff-1 do begin
         ntemp = ntemp / nsteps[ic]
         j = fix(itemp/ntemp)
         tempset.coeff[ic] = coeff_lo[ic] + (dcoeff[ic]/nsteps[ic]) * j
         itemp = itemp - j * ntemp
      endfor

      ; Construct the simulated arc spectrum
      traceset2xy, tempset, xtemp, loglambda
      model = fltarr(npix+2*npad)

      if (loglambda[1] GT loglambda[0]) then begin ; Ascending wavelengths
         for iline=0, nline-1 do begin
            qless = (loglam[iline] LT loglambda)
            iloc = (where(qless))[0] - 1
            if (iloc GE 1 AND iloc LE npix-2) then begin
               dx = loglam[iline] - loglambda[iloc]
               dpix = loglambda[iloc+1] - loglambda[iloc]
               model[npad+iloc:npad+iloc+1] = model[npad+iloc:npad+iloc+1] $
                + intensity[iline] * [1-dx/dpix, dx/dpix]
            endif
         endfor
      endif else begin ; Descending wavelengths
         for iline=0, nline-1 do begin
            qless = (loglam[iline] GT loglambda)
            iloc = (where(qless))[0] - 1
            if (iloc GE 1 AND iloc LE npix-2) then begin
               dx = loglambda[iloc] - loglam[iline]
               dpix = loglambda[iloc] - loglambda[iloc+1]
               model[npad+iloc:npad+iloc+1] = model[npad+iloc:npad+iloc+1] $
                + intensity[iline] * [1-dx/dpix, dx/dpix]
            endif
         endfor
      endelse

      ; Smooth the model spectrum.
      ; Also, take the square-root of the intensity.
      model = convol(model, gausskern, /center, /edge_truncate)
      model = sqrt(model > 1) - 1 ; Weight by the square-root of the intensity

      ; Test for the cross-correlation being ill-defined...
      if (stddev(speccorr) EQ 0 OR stddev(model) EQ 0) then begin
         corrval = -1
         icorr = 0
      endif else begin
         corrval = max( c_correlate(speccorr, model, lags), icorr)
      endelse

      if (corrval GT bestcorr) then begin
         bestcorr = corrval
         bestlambda = loglambda
         lagbest = lags[icorr]
; DEBUG PLOTS
; splot,speccorr,xr=[npad,npix+npad]
; soplot,shift(model,-lagbest)*mean(speccorr)/mean(model),color='red'
; print,bestcorr,lagbest,tempset.coeff
      endif

      ; Schlegel counter of step number...
      print, format='("Step ",i5," of ",i5,a1,$)', $
       istep, nsteptot, string(13b)

   endfor

   ; Convert to a trace set with XMIN=0, XMAX=NPIX-1
   xy2traceset, dindgen(npix)-lagbest, bestlambda, wset, ncoeff=ncoeff, $
    xmin=0, xmax=npix-1, maxiter=1

   return, wset
end

;------------------------------------------------------------------------------

function arcfit_guess, spec, loglam, intensity, color=color, func=func, $
 bestcorr=bestcorr, acoeff=acoeff, dcoeff=dcoeff, nsteps=nsteps

   if (NOT keyword_set(func)) then func = 'legendre'

   npix = N_elements(spec)

   ;---------------------------------------------------------------------------
   ; INITIAL WAVELENGTH SOLUTION
   ;---------------------------------------------------------------------------

   ; Give arcfit_iter initial starting point for wavelength solutions.

   if not (keyword_set(acoeff) and keyword_set(dcoeff)) then begin 
      if (color EQ 'blue') then begin
;      acoeff = [3.6846, -0.1060, -0.0042, 0.00012] ; Blue-1 (01)
;      acoeff = [3.7014, -0.1028, -0.0040, 0.00020] ; Blue-2 (03)
         acoeff = [3.6930, -0.1044, -0.0041, 0.00016]
         dcoeff = [0.0500,  0.0080,  0.0003, 0.00010]
      endif else if (color EQ 'red') then begin
;      acoeff = [ 3.8640, 0.1022, -0.0044, -0.00024] ; Red-1 (01)
;      acoeff = [ 3.8740, 0.0994, -0.0043, -0.00020] ; Red-2 (02)
;      acoeff = [ 3.8808, 0.0980, -0.0044, -0.00023] ; Another Red-2 (02)
         acoeff = [ 3.8700, 0.1008, -0.0044, -0.00022]
         dcoeff = [ 0.0500, 0.0080,  0.0003,  0.00010]
      endif
   endif

   nacoeff = N_elements(acoeff)
   aset = tset_struc(func, nacoeff)
   aset.xmin = 0
   aset.xmax = npix-1
   aset.coeff = acoeff

   ; First search only varies the first 3 coefficients
   splog, 'Searching about coefficients = ', aset.coeff, format='(a,99f9.5)'
   if not keyword_set(nsteps) then nsteps = [1, 20, 5, 1]
   wset = arcfit_iter(spec, loglam, intensity, $
    aset, dcoeff, nsteps, nsmooth=6.0, dlag=2, bestcorr=bestcorr)

   ; Second search varies 4 coefficients, but narrows search window
   ; for the first two coefficients
   splog, 'Searching about coefficients = ', wset.coeff, format='(a,99f9.5)'
   nsteps = [1, 5, 5, 5]
   dcoeff[0] = dcoeff[0] / 5. ; Narrow search on coefficient #0
   dcoeff[1] = dcoeff[1] / 6. ; Narrow search on coefficient #1
   wset = arcfit_iter(spec, loglam, intensity, $
    wset, dcoeff, nsteps, nsmooth=4.0, dlag=1, bestcorr=bestcorr)

   splog, 'Initial wavelength fit = ', wset.coeff, format='(a,99f9.5)'
   splog, 'Initial correlation = ', bestcorr

   return, wset
end
;------------------------------------------------------------------------------
