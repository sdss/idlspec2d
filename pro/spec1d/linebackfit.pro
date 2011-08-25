;+
; NAME:
;   linebackfit
;
; PURPOSE:
;   Fit emission lines and background terms to a spectrum.
;
; CALLING SEQUENCE:
;   result = linebackfit(lambda, loglam, flux, invvar=, $
;    linename=, zindex=, windex=, findex=, fvalue=, zlimits=, siglimits=, $
;    background=, contrange=, zguess=, sigguess=, backguess=, $
;    yfit=, bfit=, bterms=, /silent )
;
; INPUTS:
;   lambda     - Rest-frame vacuum wavelength of lines to fit in Ang [NLINE]
;   loglam     - Wavelengths for spectrum in vacuum log-10 Angstroms [NPIX]
;   flux       - Flux for spectrum [NPIX]
;   invvar     - Inverse variance for spectrum [NPIX]
;
; OPTIONAL INPUTS:
;   linename   - String name(s) of line(s) to be copied to output structure.
;   zindex     - Lines with the same ZINDEX are constrained to have the
;                same redshift; default to a unique list of indices [NLINE].
;   windex     - Lines with the same WINDEX are constrained to have the
;                same width [NLINE]; default to a unique list of indices.
;   findex     - Lines with the same FINDEX are constrained to have the
;                same flux ratio as input in FVALUE; defult to a unique
;                list of indices [NLINE].
;   fvalue     - If FINDEX is specified, then constrain the lines to have
;                the same flux ratios as in FVALUE [NLINE]; default to
;                all values of 1.  These values are also used as the
;                initial guesses for the line strengths.
;   zlimits    - Optional limits in velocity, either as a 2-elements array
;                for the low/high limit of every line, or as a [2,NLINE] array
;                with different limits for each line [km/s]
;   siglimits  - Optional limits in sigma, either as a 2-elements array
;                for the low/high limit of every line, or as a [2,NLINE] array
;                with different limits for each line [km/s]
;   background - Background vector(s) to fit simultaneously with the lines,
;                where the scaling of each vector is fit [NPIX,NBACK].
;                The redshift of these background vectors is not fit, but
;                rather we maintain the one-to-one correspondence of each
;                BACKGROUND pixel to each FLUX pixel.  The initial guess
;                for the scaling of each background vector is unity.
;   contrange  - Continuum-level fitting range.  If not set or set to zero,
;                then use the background level at the line center.
;                Or set as 2-element array for fitting ranges on either side
;                of the line center in km/s.  A mean (unweighted) continuum
;                is evaluated in the blueward and redward bands, and then
;                the average of those two bands.
;   zguess     - Initial guess for redshifts of all lines (scalar or vector
;                with one entry per line); default to 0.
;   sigguess   - Initial guess for sigmas of all lines in log-10 Angstroms
;                (scalar or vector with one entry per line); default to 1.5d-4
;                (105 km/sec).
;   backguess  - Initial guess for background term coefficients [NBACK];
;                default to 1 for all coefficients
;   silent     - If set, then make no print statements with SPLOG
;
; OUTPUTS:
;   result     - Output structure with result of line fits [NLINE].
;                The elements are as follows:
;                LINENAME        - String name of line copied from input
;                                  parameter by the same name, or else
;                                  constructed from the rounded-down wavelength
;                                  of each line in air (not vacuum)
;                LINEWAVE        - Rest-frame wavelength in vacuum [Ang]
;                                  (copy of LAMBDA)
;                LINEZ           - Redshift [dimensionless]
;                LINEZ_ERR       - Error in above
;                LINESIGMA       - Sigma of gaussian [km/s]
;                LINESIGMA_ERR   - Error in above
;                LINEAREA        - Area of line [flux-units * Ang]
;                LINEAREA_ERR    - Error in above
;                LINEEW          - Line equivalent width [Ang]; or 0 if the
;                                  background level is <= 0.
;                LINEEW_ERR      - Error in above; or -1L if the background
;                                  level is <= 0.  We approximate this error
;                                  as sqrt((LINEAREA_ERR/LINEAREA)^2
;                                         +(LINECONTLEVEL_ERR/LINECONTLEVEL)^2) * LINEEW
;                LINECONTLEVEL   - Continuum level at line center [flux-units];
;                                  if the line center is outside the wavelength
;                                  range, then return the nearest value (either
;                                  the first or last value)
;                LINECONTLEVEL_ERR - Error in above, or -1L if the line center
;                                  is outside the wavelength range.
;                LINENPIXLEFT    - Number of pixels from the line center to
;                                  -3 sigma that have INVVAR > 0.
;                LINENPIXRIGHT   - Number of pixels from the line center to
;                                  +3 sigma that have INVVAR > 0.
;                LINEDOF         - LINENPIXLEFT+LINENPIXRIGHT minus the number
;                                  of terms fit for that line, which could be
;                                  fractional (if one parameter is fixed between
;                                  N lines, then we say only 1/N-th of that
;                                  parameter is fit in each of those lines).
;                                  This can be zero or negative.
;                LINECHI2       -  Chi^2 for all points within +/- 3 sigma of
;                                  the line center; -1L if no such points.
;
; OPTIONAL OUTPUTS:
;   yfit       - Fit spectrum including lines and background terms [NPIX].
;   bfit       - Fit spectrum including only background terms [NPIX].
;   bterms     - Coefficients for background terms [NBACK].
;
; COMMENTS:
;   If a line was dropped from the fit (for example, no points to fit),
;   then set the LINEAREA to 0 and the LINEAREA_ERR to -1L.
;
;   Also, if LINENPIX=0 for a line, then remove that line from the fit.
;
; EXAMPLES:
;
; BUGS:
;   Do not use lines with no points to fit in the computation of
;   degrees of freedom for other lines.
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   create_linestruct()
;   mpfitfun
;
; INTERNAL SUPPORT ROUTINES:
;   onegauss()
;   manygauss()
;
; REVISION HISTORY:
;   05-Feb-2002  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function onegauss, xval, pp, dp

   term1 = exp( - (xval - pp[1])^2 / (2. * pp[2]^2) )
   yval = pp[0] * term1 / (sqrt(2.*!pi) * pp[2])
   if (arg_present(dp)) then begin
      message, 'Derivatives not supported yet!!!???'
; These derivatives are for a different parameterization...
;      dp0 = term1
;      dp1 = yval * ((xval - pp[1]) / pp[2]^2)
;      dp2 = yval * ((xval - pp[1])^2 / pp[2]^3)
;      dp = [dp0,dp1,dp2]
   endif

   return, yval
end
;------------------------------------------------------------------------------
function manygauss, xindx, pp, nline=nline, nback=nback, loglam=loglam, $
 background=background

   yval = 0.d
   for iline=0L, nline-1L do $
    yval = yval + onegauss(loglam[xindx], pp[iline*3:iline*3+2])
   for iback=0L, nback-1L do $
    yval = yval + background[xindx,iback] * pp[nline*3+iback]

   return, yval
end
;------------------------------------------------------------------------------
function linebackfit, lambda, loglam, flux, invvar=invvar, linename=linename, $
 zindex=zindex, windex=windex, findex=findex, fvalue=fvalue, $
 zlimits=zlimits1, siglimits=siglimits1, $
 background=background, contrange=contrange, $
 zguess=zguess1, sigguess=sigguess1, $
 backguess=backguess1, yfit=yfit, bfit=bfit, bterms=bterms, silent=silent

   cspeed = 2.99792458e5

   nline = n_elements(lambda)
   ndim = size(background, /n_dimen)
   dims = size(background, /dimens)
   if (ndim EQ 1) then nback = 1 $
    else if (ndim EQ 2) then nback = dims[1] $
    else nback = 0
   if (n_elements(background) EQ 0) then background = 0

   if (keyword_set(invvar)) then $
    if (n_elements(loglam) NE n_elements(invvar)) then $
     message, 'Number of elements in LOGLAM and INVVAR must agree'
   if (keyword_set(background)) then $
    if (n_elements(loglam) NE dims[0]) then $
     message, 'Number of elements in LOGLAM and BACKGROUND[*,0] must agree'

   if (NOT keyword_set(zindex)) then zindex = lindgen(nline)
   if (NOT keyword_set(windex)) then windex = lindgen(nline)
   if (NOT keyword_set(findex)) then findex = lindgen(nline)
   if (NOT keyword_set(fvalue)) then fvalue = fltarr(nline) + 1.
   if (keyword_set(zguess1)) then begin
      if (n_elements(zguess1) EQ 1) then begin
         zguess = replicate(zguess1[0], nline)
      endif else if (n_elements(zguess1) EQ nline) then begin
         zguess = zguess1
      endif else begin
         message, 'Wrong number of elements for ZGUESS'
      endelse
   endif else begin
      zguess = replicate(0.d0, nline)
   endelse
   if (keyword_set(sigguess1)) then begin
      if (n_elements(sigguess1) EQ 1) then begin
         sigguess = replicate(sigguess1[0], nline)
      endif else if (n_elements(sigguess1) EQ nline) then begin
         sigguess = sigguess1
      endif else begin
         message, 'Wrong number of elements for SIGGUESS'
      endelse
   endif else begin
      sigguess = replicate(1.5d-4, nline)
   endelse
   if (nback GT 0) then begin
      if (keyword_set(backguess1)) then begin
         if (n_elements(backguess1) EQ 1) then begin
            backguess = replicate(backguess1[0], nback)
         endif else if (n_elements(backguess1) EQ nback) then begin
            backguess = backguess1
         endif else begin
            message, 'Wrong number of elements for BACKGUESS'
         endelse
      endif else begin
         backguess = replicate(1., nback)
      endelse
   endif

   npix = n_elements(flux)

   ;----------
   ; Initialize the output structure, and set default return values

   linestruct = create_linestruct(nline)

   linestruct.linewave = lambda
   if (keyword_set(linename)) then begin
      linestruct.linename = linename
   endif else begin
      airwave = lambda
      vactoair, airwave
      linestruct.linename = strtrim(string(long(airwave),format='(i)'),2)
   endelse

   ;----------
   ; Initialize the structures to be passed to the fitting routine

   parinfo = replicate({value:0.D, fixed:0, limited:[0,0], tied:'', $
    limits:[0.D,0]}, nline*3+nback)
   functargs = {nline: nline, nback: nback, loglam: loglam, $
    background: background}

   ;----------
   ; Set the initial guesses

   for iline=0L, nline-1L do begin
      parinfo[0+iline*3].value = fvalue[iline]
      parinfo[1+iline*3].value = alog10( lambda[iline] * (1. + zguess[iline]) )
      parinfo[2+iline*3].value = sigguess[iline]
   endfor

   ; Set the initial guess for the scaling of each background vector
   parinfo[nline*3+lindgen(nback)].value = backguess

   ;----------
   ; Make a list of the number of fitting terms per line.
   ; If a parameter is constrained between N lines, then we say each
   ; of those lines is only fitting 1/N-th of that parameter

   nfitterms = fltarr(nline)

   ;----------
   ; Apply constraints to peak flux values (not integrated flux in a line)

   allindex = findex[ uniq(findex,sort(findex)) ]
   for iall=0L, n_elements(allindex)-1L do begin
      ii = where(findex EQ allindex[iall], ct)
      nfitterms[ii] = nfitterms[ii] + 1.0/ct
      if (ct GT 1) then begin
         for jj=1L, ct-1L do begin
            fratio = fvalue[ii[jj]] / fvalue[ii[0]]
            parinfo[0+ii[jj]*3].tied = $
             string(fratio, 0+ii[0]*3, $
             format='(e12.5," * P(",i,")")')
         endfor
      endif
   endfor

   ;----------
   ; Apply constraints to couple redshifts

   allindex = zindex[ uniq(zindex,sort(zindex)) ]
   for iall=0L, n_elements(allindex)-1L do begin
      ii = where(zindex EQ allindex[iall], ct)
      nfitterms[ii] = nfitterms[ii] + 1.0/ct
      if (ct GT 1) then begin
        for jj=1L, ct-1L do begin
            lamshift = alog10(lambda[ii[jj]] / lambda[ii[0]])
            parinfo[1+ii[jj]*3].tied = $
             string(lamshift, 1+ii[0]*3, $
             format='(e12.5," + P(",i,")")')
         endfor
      endif
   endfor

   ;----------
   ; Apply constraints to couple widths

   allindex = windex[ uniq(windex,sort(windex)) ]
   for iall=0L, n_elements(allindex)-1L do begin
      ii = where(windex EQ allindex[iall], ct)
      nfitterms[ii] = nfitterms[ii] + 1.0/ct
      if (ct GT 1) then begin
        for jj=1L, ct-1L do begin
            parinfo[2+ii[jj]*3].tied = $
             string(2+ii[0]*3, format='("P(",i,")")')
         endfor
      endif
   endfor

   ;----------
   ; Constrain min/max redshift relative to initial guess

   if (keyword_set(zlimits1)) then begin
      if (size(zlimits1,/n_dimen) EQ 1) then $
       zlimits = rebin(zlimits1, 2, nline) $
      else $
       zlimits = zlimits1
      parinfo[lindgen(nline)*3+1].limits[0] = $
       alog10(lambda * (1. + reform(zlimits[0,*])))
      parinfo[lindgen(nline)*3+1].limits[1] = $
       alog10(lambda * (1. + reform(zlimits[1,*])))
      parinfo[lindgen(nline)*3+1].limited = 1
   endif

   ;----------
   ; Constrain min/max sigmas relative to initial guess

   if (keyword_set(siglimits1)) then begin
      if (size(siglimits1,/n_dimen) EQ 1) then $
       siglimits = rebin(siglimits1, 2, nline) $
      else $
       siglimits = siglimits1
      parinfo[lindgen(nline)*3+2].limits[0] = $
       reform(siglimits[0,*]) / (alog(10.)*cspeed)
      parinfo[lindgen(nline)*3+2].limits[1] = $
       reform(siglimits[1,*]) / (alog(10.)*cspeed)
      parinfo[lindgen(nline)*3+2].limited = 1
   endif

   ;----------
   ; Do the fit!

   yfit = fltarr(npix)

   igood = where(invvar GT 0, ngood)
   status = 0
   if (ngood GT 0) then begin
      lfit = mpfitfun('manygauss', igood, flux[igood], $
       1./sqrt(invvar[igood]), parinfo=parinfo, $
       covar=covar, perror=perror, yfit=yfit1, functargs=functargs, $
       nfev=nfev, niter=niter, status=status, /quiet)

      if (NOT keyword_set(silent)) then begin
         splog, 'MPFIT nfev=', nfev, ' niter=', niter, ' status=', status
         if (status EQ 5) then $
          splog, 'Warning: Maximum number of iterations reached: ', niter
      endif
      yfit[igood] = yfit1
   endif
   if (ngood EQ 0 OR status EQ 0) then begin
      if (NOT keyword_set(silent)) then $
       splog, 'Too few points to fit ', ngood
      nparam = 3 * nline + nback
      lfit = fltarr(nparam)
      perror = fltarr(nparam)
      covar = fltarr(nparam,nparam)
   endif

   ;----------
   ; For parameters that are fixed, assign them the same errors as
   ; those parameters to which they are tied.  For the flux area,
   ; scale those errors according to the ratio of the flux levels.

   allindex = findex[ uniq(findex,sort(findex)) ]
   for iall=0L, n_elements(allindex)-1L do begin
      ii = where(findex EQ allindex[iall], ct)
      if (ct GT 1) then perror[ii*3+0] = perror[ii[0]*3+0] $
       * fvalue[ii] / fvalue[ii[0]]
   endfor

   allindex = zindex[ uniq(zindex,sort(zindex)) ]
   for iall=0L, n_elements(allindex)-1L do begin
      ii = where(zindex EQ allindex[iall], ct)
      if (ct GT 1) then perror[ii*3+1] = perror[ii[0]*3+1]
   endfor

   allindex = windex[ uniq(windex,sort(windex)) ]
   for iall=0L, n_elements(allindex)-1L do begin
      ii = where(windex EQ allindex[iall], ct)
      if (ct GT 1) then perror[ii*3+2] = perror[ii[0]*3+2]
   endfor

   ;----------
   ; Construct the line-measure outputs (and their errors)

   linestruct.linearea      = lfit[lindgen(nline)*3+0] $
    * alog(10.) * 10.^lfit[lindgen(nline)*3+1]
   linestruct.linez         =  (10.^lfit[lindgen(nline)*3+1] / lambda - 1) $
    * (lfit[lindgen(nline)*3+1] GT 0) ; Set to to zero if the parameter is zero
   linestruct.linesigma     = abs(lfit[lindgen(nline)*3+2]) * alog(10.) * cspeed
   linestruct.linearea_err  = perror[lindgen(nline)*3+0] $
    * alog(10.) * 10.^lfit[lindgen(nline)*3+1]
   linestruct.linez_err     = perror[lindgen(nline)*3+1] $
    * alog(10.) * (linestruct.linez + 1)
   linestruct.linesigma_err = perror[lindgen(nline)*3+2] * alog(10.) * cspeed

   ;----------
   ; If a line was dropped from the fit (for example, no points to fit),
   ; then set the LINEAREA to 0 and the LINEAREA_ERR to -1L.

   ibad = where(perror[lindgen(nline)*3+0] LE 0)
   if (ibad[0] NE -1) then begin
      linestruct[ibad].linearea = 0
      linestruct[ibad].linearea_err = -1L
   endif

   ;----------
   ; Find the background levels only

   if (nback EQ 0) then begin
      bterms = 0
      bfit = fltarr(npix)
      berr = fltarr(npix)
   endif else begin
      bterms = lfit[nline*3:nline*3+nback-1]
      bcovar = covar[nline*3:nline*3+nback-1,nline*3:nline*3+nback-1]

;      The following two methods for evaluating bfit are equivalent.
;      bfit = manygauss(lindgen(npix), bterms, nline=0, nback=nback, $
;       loglam=loglam, background=background)
      bfit = bterms ## background

      berr = fltarr(npix)
      for ipix=0L, npix-1L do $
       berr[ipix] = sqrt( (transpose(background[ipix,*]) $
        # bcovar # transpose(background[ipix,*])) )

      ; If the returned uncertainty for a template is ill-defined in MPFIT,
      ; then any non-zero pixels for that template should propogate to
      ; an ill-defined template fit
      for iback=0L, nback-1L do begin
         if (perror[nline*3+iback] LE 0) then $
          berr *= (background[*,iback] EQ 0)
      endfor
   endelse

   ;----------
   ; For each line, determine the background level at the line center
   ; and the number of pixels and chi^2 of each line fit.

   logmin = min(loglam)
   logmax = max(loglam)
   for iline=0L, nline-1L do begin
      junk = min(abs(loglam - lfit[iline*3+1]), ipix)
      if (keyword_set(contrange)) then begin
         if (n_elements(contrange) NE 2) then $
          message, 'CONTRANGE must be 2-element array!'
         if (contrange[0] LT 0 OR contrange[1] LT 0 $
          OR contrange[0] GT contrange[1]) then $
          message, 'CONTRANGE has invalid values!'
         loglam1 = lfit[iline*3+1] - alog10(1+contrange[1]/cspeed)
         loglam2 = lfit[iline*3+1] - alog10(1+contrange[0]/cspeed)
         loglam3 = lfit[iline*3+1] - alog10(1-contrange[0]/cspeed)
         loglam4 = lfit[iline*3+1] - alog10(1-contrange[1]/cspeed)
         junk = min(abs(loglam - loglam1), ipix1)
         junk = min(abs(loglam - loglam2), ipix2)
         junk = min(abs(loglam - loglam3), ipix3)
         junk = min(abs(loglam - loglam4), ipix4)
         cont_left = mean(bfit[ipix1:ipix2])
         cont_right = mean(bfit[ipix3:ipix4])
         if (loglam2 LT logmin) then begin
            ; Case where the left continuum band is blueward of the spectrum.
            linestruct[iline].linecontlevel = cont_right
            linestruct[iline].linecontlevel_err = -1L
         endif else if (loglam3 GT logmax) then begin
            ; Case where the right continuum band is redward of the spectrum.
            linestruct[iline].linecontlevel = cont_left
            linestruct[iline].linecontlevel_err = -1L
         endif else begin
            ; Evaluate the mean of the left and right continuum bands.
            linestruct[iline].linecontlevel = 0.5 * (cont_left + cont_right)
            linestruct[iline].linecontlevel_err = berr[ipix] ; ???
         endelse
      endif else begin
         if (lfit[iline*3+1] LT logmin) then begin
            ; Case where the line center is blueward of the entire spectrum.
            linestruct[iline].linecontlevel = bfit[0]
            linestruct[iline].linecontlevel_err = -1L
         endif else if (lfit[iline*3+1] GT logmax) then begin
            ; Case where the line center is redward of the entire spectrum.
            linestruct[iline].linecontlevel = bfit[npix-1]
            linestruct[iline].linecontlevel_err = -1L
         endif else begin
            ; Select the nearest pixel for evaluating the background
            ; level at this line center.
            linestruct[iline].linecontlevel = bfit[ipix]
            linestruct[iline].linecontlevel_err = berr[ipix]
         endelse
      endelse
   endfor

   ;----------
   ; Find the pixels that are within +/- 3 sigma of the line center
   ; Note that if the line is very (unphysically) narrow, it is
   ; possible to have no lines within this domain.  Reject those fits.

   for iline=0L, nline-1L do begin
      qleft = loglam GE lfit[iline*3+1] - 3 * lfit[iline*3+2] $
       AND loglam LE lfit[iline*3+1]
      qright = loglam LE lfit[iline*3+1] + 3 * lfit[iline*3+2] $
       AND loglam GT lfit[iline*3+1]
      qgood = invvar GT 0
      linestruct[iline].linenpixleft = total(qleft AND qgood)
      linestruct[iline].linenpixright = total(qright AND qgood)
      linenpix = linestruct[iline].linenpixleft $
       + linestruct[iline].linenpixright

      if (linenpix GT 0) then begin
         linestruct[iline].linedof = linenpix - nfitterms[iline]

         if (linestruct[iline].linedof GT 0) then begin
            indx = where(qgood AND (qleft OR qright))
            linestruct[iline].linechi2 = $
             total( (flux[indx] - yfit[indx])^2 * invvar[indx] )
         endif
      endif

      ; Special-case rejection -- set AREA=0,AREA_ERR=-2
      ; if there are no data points within the line-fitting region.
      if (linestruct[iline].linenpixleft + linestruct[iline].linenpixright $
       EQ 0) then begin
         linestruct[iline].linearea = 0
         linestruct[iline].linearea_err = -2L

         ; Set these line-fit coefficients equal to zero for when
         ; we re-evaluate YFIT.
         lfit[iline*3+0] = 0
      endif
   endfor

   ;----------
   ; Construct the line equivalent widths

   for iline=0L, nline-1L do begin
      if (linestruct[iline].linenpixleft + linestruct[iline].linenpixright $
       LE 0) then begin
         linestruct[iline].lineew_err = -2L
      endif else if (linestruct[iline].linecontlevel LE 0) then begin
         linestruct[iline].lineew_err = -1L
      endif else begin
         linestruct[iline].lineew = linestruct[iline].linearea $
          / linestruct[iline].linecontlevel

         ; The following is a crude approximation of the EW error,
         ; assuming that the line area + continuum level are uncorrelated.
         if (linestruct[iline].linearea_err LE 0) then begin
            linestruct[iline].lineew_err = -3L
         endif else if (linestruct[iline].linecontlevel_err LE 0) then begin
          linestruct[iline].lineew_err = -3L
         endif else begin
            linestruct[iline].lineew_err = linestruct[iline].lineew * $
             (linestruct[iline].linearea_err / linestruct[iline].linearea $
             - linestruct[iline].linecontlevel_err $
             / linestruct[iline].linecontlevel )
            if (linestruct[iline].lineew_err LT 0) then $
             linestruct[iline].lineew_err = -4L
         endelse
      endelse
   endfor

   ;----------
   ; Re-evaluate such that we get the functional fit at the rejected
   ; wavelengths too.  This also re-evaluates the fit for lines that
   ; have been rejected and removed due to no valid data points within
   ; the line-fitting region.

   if (arg_present(yfit)) then $
    yfit = manygauss(lindgen(npix), lfit, _EXTRA=functargs)

   return, linestruct
end
;------------------------------------------------------------------------------
