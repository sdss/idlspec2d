;+
; NAME:
;   fitarcimage
;
; PURPOSE:
;   Determine wavelength calibration from arclines
;
; CALLING SEQUENCE:
;   fitarcimage, arc, arcinvvar, color, linelist, xnew, ycen, wset, $
;    func=func, aset=aset, ncoeff=ncoeff, lambda=lambda, $
;    thresh=thresh, row=row, $
;    xdif_lfit=xdif_lfit, xdif_tset=xdif_tset
;
; INPUTS:
;   arc        - Extracted arc spectra with dimensions [NY,NFIBER]
;   arcinvvar  - Inverse variance of ARC
;   color      - 'red' or 'blue'; not required if ANS is set
;   linelist   - Array with wavelengths of possible arc lines
;
; OPTIONAL KEYWORDS:
;   func       - Name of fitting function; default to 'legendre'
;   aset       - Trace set for initial wavelength solution in row number ROW.
;   ncoeff     - Number of coefficients in fits.  This may be different than
;                the number of coefficients in the initial guess ASET.
;                Default to 5.
;   thresh     - Threshhold counts for significant lines;
;                default to 200 if COLOR='blue' or 500 if COLOR='red'
;   row        - Row to use in initial guess of wavelength solution;
;                default to (NFIBER-30)/2
;
; OUTPUTS:
;   xnew       - pixel position of lines [nfiber, nlambda]
;   ycen       - fiber number [nfiber, nlambda]
;   wset       - traceset (pix -> lambda)
;
; OPTIONAL OUTPUTS:
;   lambda     - returns alog10(wavelength) of good lines
;   xdif_lfit  - fit residual for individual arclines
;   xdif_tset  - fit residual of traceset
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Check for case where initial ROW is an un-illuminated fiber.
;   Not making sure that only the same lines are fit for each fiber.
;      (Different lines can be rejected in xy2traceset.)
;   THRESH is unused.
;   TRACESET2PIX maybe returns the transpose of what is natural?
;   Check QA stuff at end.
;   When constructing MX, exclude any fibers with bad (e.g., saturated) pixels.
;
; INTERNAL PROCEDURES:
;   tset_struc()
;   fullfit()
;   fitmx()
;
; PROCEDURES CALLED:
;   djs_median
;   djsig()
;   fit_tset
;   fitwithmx()
;   trace_crude()
;   trace_fweight()
;   traceset2pix()
;   xy2traceset
;
; REVISION HISTORY:
;   15-Oct-1999  Written by S. Burles, D. Finkbeiner, & D. Schlegel, APO.
;   09-Nov-1999  Major modifications by D. Schlegel, Ringberg.
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
; [COEFF-0.5*DCOEFF, COEFF+0.5*DCOEFF]
; except for the first coefficient, for which NSTEPS is unused.

function fullfit, spec, linelist, aset, dcoeff, nsteps, bestcorr=bestcorr

   npix = N_elements(spec)
   logline = (alog10(linelist[*,0]))[*]
   intensity = (linelist[*,1])[*]
   nline = N_elements(logline)
   pixarray = 2.0d0 * dindgen(npix) / (npix-1) - 1.0d0

   ncoeff = N_elements(aset.coeff)

   nsmooth = 5
   nsz = 4*fix(nsmooth) + 1 ; kernal size (odd)
   gausskern = exp( -( ((nsz-1)/2 - findgen(nsz))^2 ) / nsmooth^2 )
   gausskern = gausskern / total(gausskern)

   ; Clip, pad, and smooth input spectrum
; NOTE THAT THIS MEDIAN-FILTERING CLIPS DATA WITHIN 50 PIXELS OF THE EDGE !???
   speccorr = spec
   speccorr = speccorr - median(speccorr,101) > 1
   speccorr = [speccorr, fltarr(npix)+1]
   speccorr = convol(speccorr, gausskern, /center, /edge_truncate)
   speccorr = sqrt(speccorr > 1) ; Weight by the square-root of the intensity

   bestcorr = -1.0
   bestcoeff = 0
   coeff_lo = aset.coeff - dcoeff * (nsteps-1)/(2.*nsteps)
   coeff_lo[0] = aset.coeff[0] ; Set this coefficient to its central value,
                               ; e.g. that passed in ASET.  We will be doing
                               ; the cross-correlation left and right of this.
   tempset = aset

   ; Loop over all coefficients except the first one

   ; Set minimum number of lags to check equal to 10 pixels
   nlag = fix( dcoeff[0] / (2. * abs(aset.coeff[1]) / npix) + 1 ) > 10
print, 'nlag', nlag
   lags = indgen(nlag) - fix(nlag/2)
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
      model = fltarr(2*npix)

      if (loglambda[1] GT loglambda[0]) then begin ; Ascending wavelengths
         for iline=0, nline-1 do begin
            qless = (logline[iline] LT loglambda)
            iloc = (where(qless))[0]
            if (iloc GE 1 AND iloc LE npix-2) then begin
               dx = logline[iline] - loglambda[iloc]
               dpix = loglambda[iloc+1] - loglambda[iloc]
               model[iloc:iloc+1] = model[iloc:iloc+1] $
                + intensity[iline] * [1-dx, dx] / dpix
            endif
         endfor
      endif else begin ; Descending wavelengths
         for iline=0, nline-1 do begin
            qless = (logline[iline] GT loglambda)
            iloc = (where(qless))[0]
            if (iloc GE 1 AND iloc LE npix-2) then begin
               dx = loglambda[iloc] - logline[iline]
               dpix = loglambda[iloc] - loglambda[iloc+1]
               model[iloc:iloc+1] = model[iloc:iloc+1] $
                + intensity[iline] * [1-dx, dx] / dpix
            endif
         endfor
      endelse

      model = convol(model, gausskern, /center, /edge_truncate)
      model = sqrt(model > 1) - 1 ; Weight by the square-root of the intensity

      corrval = max( c_correlate(speccorr, model, lags), icorr)
      if (corrval GT bestcorr) then begin
         bestcorr = corrval
         bestlambda = loglambda
         lagbest = lags[icorr]
; PLOT ???
splot,speccorr,xr=[0,2048]
soplot,shift(model,-lagbest)*mean(speccorr)/mean(model),color='red'
print,bestcorr,lagbest,tempset.coeff
      endif

   endfor
print, 'Correlation = ', bestcorr

   ; Convert to a trace set with XMIN=0, XMAX=NPIX-1
   xy2traceset, dindgen(npix)-lagbest, bestlambda, wset, ncoeff=ncoeff, $
    xmin=0, xmax=npix-1, maxiter=0

; TEST ???
;traceset2xy,wset,xtemp,ytemp
;print, traceset2pix(wset,alog10([5769.6,5790.663,5852.48]))

   return, wset
end

;------------------------------------------------------------------------------
; LAMBDA = log10-wavelength

function fitmx, wset, lambda, xpos, nord=nord

   if (NOT keyword_set(nord)) then nord = 4
   dims = size(xpos, /dim)
   nfiber = dims[0]
   nlambda = dims[1]

   ; Evaluate the trace set to get the fit pixel position at each arc lambda
   pix1 = traceset2pix(wset, lambda)

   x = findgen(nfiber) / float(nfiber) ; From [0,1)
   xnew = xpos

   ; Loop through each arc line...

   for i=0, nlambda-1 do begin
      mx = transpose(pix1[i,*]) ; Pixel position where wavelength should fall
      dif = xpos[*,i] - mx ; Measured position minus MX

      ; Fit to DIF as a function of fiber number
      junk = poly_fit(x, dif, nord, yfit)
      res1 = yfit - dif

      ; Find which points are most deviant (4-sigma cut; 2 iterations)
      qgood = abs(res1) LT 4.*stdev(res1)
      qgood = abs(res1) LT 4.*stdev(res1*qgood)

      ; Re-fit to DIF as a function of fiber number (rejecting outliers)
      junk = polyfitw(x, dif, qgood, nord, yfit)

      ; Return the position where the wavelength should fall (MX) plus
      ; a fitted-function to the deviations from those positions.
      xnew[*,i] = mx + yfit
   endfor

   return, xnew
end
;------------------------------------------------------------------------------

pro fitarcimage, arc, arcinvvar, color, linelist, xnew, ycen, wset, $
 func=func, aset=aset, ncoeff=ncoeff, lambda=lambda, thresh=thresh, row=row, $
 xdif_lfit=xdif_lfit, xdif_tset=xdif_tset

   ;---------------------------------------------------------------------------
   ; Preliminary stuff
   ;---------------------------------------------------------------------------

   if (NOT keyword_set(aset)) then begin
      if (color NE 'blue' AND color NE 'red') then $
       message, "SIDE must be set to 'blue' or 'red' if ANS not specified"
   endif
   if (NOT keyword_set(func)) then func = 'legendre'
   if (NOT keyword_set(ans)) then ans = 0
   if (func EQ 'legendre') then function_name = 'flegendre'
   if (func EQ 'chebyshev') then function_name = 'fchebyshev'

   t_begin = systime(1)

   ndim = size(arc, /n_dim)
   if (ndim NE 2) then $
    message, 'Expecting 2-D arc image'
   dims = size(arc, /dim)
   npix = dims[0]
   nfiber = dims[1]
   if (total(dims NE size(arcinvvar, /dim))) then $
    message, 'ARC and ARCINVVAR must have same dimensions'

   if (NOT keyword_set(row)) then row = (nfiber-30)/2

   if (NOT keyword_set(thresh)) then begin
      if (color EQ 'blue') then thresh = 200
      if (color EQ 'red') then thresh = 500
   endif

   if (NOT keyword_set(ncoeff)) then ncoeff = 5

   ;---------------------------------------------------------------------------
   ; INITIAL WAVELENGTH SOLUTION
   ;---------------------------------------------------------------------------

   ; Find the initial wavelength solution if ANS is not passed
   ; One might want to change nave and nmed for first pass???

   if (NOT keyword_set(aset)) then begin

      ; Extract one spectrum from the 5 spectra around fiber number ROW
      ; by taking the median value at each wavelength

      spec = djs_median(arc[*,row-2:row+2], 2) ; ???

      ; Give fullfit initial starting point for wavelength solutions.

      if (color EQ 'blue') then begin
;         acoeff = [3.6846, -0.1060, -0.0042, 0.00012] ; Blue-1 (01)
;         acoeff = [3.7014, -0.1028, -0.0040, 0.00020] ; Blue-2 (03)
         acoeff = [3.6930, -0.1044, -0.0041, 0.00016]
         dcoeff = [0.0200,  0.0040,  0.0003, 0.00010]
         nsteps = [1, 10, 5, 5]
      endif else if (color EQ 'red') then begin
         ; For red-1 (04) or red-2 (02)
         acoeff = [ 3.8640, 0.1022, -0.0044, -0.00024]
         dcoeff = [ 0.0100, 0.0020,  0.0003,  0.00010]
         nsteps = [1, 5, 5, 5]
      endif

      nacoeff = N_elements(acoeff)
      aset = tset_struc(func, nacoeff)
      aset.xmin = 0
      aset.xmax = npix-1
      aset.coeff = acoeff
      wset = fullfit(spec, linelist, aset, dcoeff, nsteps, $
       bestcorr=bestcorr)
stop

      if (color EQ 'blue' AND bestcorr LT 0.70) then $
       print, 'Initial wavelength solution looks suspicious'
      if (color EQ 'red' AND bestcorr LT 0.70) then $
       print, 'Initial wavelength solution looks suspicious'

   endif

   ; Select which lines from the line list to trace and use in the final fits

   igood = where(linelist[*,2] GT 0.0)
   if (igood[0] EQ -1) then $
    message, 'No unblended lines in linelist'
   uselambda = alog10(linelist[igood,0])

; THIS DOES NOT WORK !!! ???
   xstart = traceset2pix(wset, uselambda)
   itrim = where(xstart GT 1 AND xstart LT npix-2)
   if (itrim[0] EQ -1) then $
    message, 'No arc lines in wavelength range'
   xstart = xstart[itrim]
   uselambda = uselambda[itrim]
;window,0
;plot,spec,xrange=[1500,1900]
;for i=0,N_elements(xstart)-1 do $ 
; djs_oplot, [xstart[i],xstart[i]],[0,1e4], color='red'
;print, traceset2pix(wset,alog10([5769.6,5790.663,5852.48]))

   ;---------------------------------------------------------------------------
   ; Trace
   ;---------------------------------------------------------------------------

   ; Allow for a shift of up to 2 pixels in the initial centers,
   ; but only 0.3 pixels while tracing

   xcen = trace_crude(arc, yset=ycen, nave=1, nmed=1, xstart=xstart, $
    ystart=row, maxshifte=0.3d, maxshift0=2.0d)

   ; Iterate the flux-weighted centers
   xnew = trace_fweight(arc, xcen, ycen)

   ; Make use of the errors??? - Seems to just mess things up???
   ; Well... the reason for that is satured lines, which return infinite errors
;   xnew = trace_fweight(arc, xcen, ycen, invvar=arcinvvar, xerr=xerr)

   ;---------------------------------------------------------------------------
   ; Reject bad (i.e., saturated) lines
   ;---------------------------------------------------------------------------

   ; Reject any arc line with more than 10% of the fibers are bad.
   ; Bad fibers are any with an infinite error (ARCINVVAR=0) within 1 pixel
   ; of the central wavelength.  Note that saturated lines should then
   ; show up as bad.

   nmatch = (size(xnew))[2] ; Number of lamp lines traced
   qbad = bytarr(nmatch)

   for i=0, nmatch-1 do begin
      xpix = round(xnew[*,i]) ; Nearest X position (wavelength) in all traces
      mivar = fltarr(nfiber) + 1
      for ix=-1, 1 do begin
         mivar = mivar * arcinvvar[ ((xpix+ix)>0)<(npix-1), * ]
      endfor
      junk = where(mivar EQ 0, nbad)
      fracbad = float(nbad) / nfiber
      qbad[i] = fracbad GT 0.10
      if (qbad[i]) then $
       print, 'Discarding trace', i, ',   fraction bad', fracbad
   endfor

   igood = where(qbad EQ 0, ngood)
   print, 'Number of good arc lines: ', ngood
   if (ngood EQ 0) then $
    message, 'No good arc lines'

   ; Trim linelist

   xnew = xnew[*,igood]
   ycen = ycen[*,igood]
   lambda = uselambda[igood]

   ;---------------------------------------------------------------------------
   ; Do the first traceset fit
   ;---------------------------------------------------------------------------

; ??? Let maxdev be a parameter; should be about 3.0d-5 = 20 km/s
maxdev = 3.0d-5

   nlambda = N_elements(lambda)
   xy2traceset, transpose(double(xnew)), lambda # (dblarr(nfiber)+1), wset, $
    func=func, ncoeff=ncoeff, maxdev=maxdev, maxiter=nlambda, /singlerej, $
    xmask=xmask, xmin=0, xmax=npix-1

   print,'Pass 1 complete'

   ;---------------------------------------------------------------------------
   ; Do the second traceset fit
   ;---------------------------------------------------------------------------

   ; Keep only "good" lines.
   ; The following logic means that an arc line is rejected if any
   ; bundle has more than 3 bad centers.

   if (nfiber NE 320) then $
    message, 'Not 320 fibers -- Cannot figure out bundle test'
   testg = reform(xmask, nlambda, 20, 16)   
   gind = where(total(total(testg EQ 0,2) GT 3,2) EQ 0)
   if (gind[0] EQ -1) then $
    message, 'No good arcs common to all fiber bundles'

   nlambda = N_elements(gind) ; number of good lines
   xnew = xnew[*,gind]
   lambda = lambda[gind]
   ycen = ycen[*,gind]

   xy2traceset, transpose(xnew), lambda # (dblarr(nfiber)+1), wset, $
    func=func, ncoeff=ncoeff, maxdev=maxdev, maxiter=nlambda, /singlerej, $
    xmask=xmask, xmin=0, xmax=npix-1
   print, 'Pass 2 complete'

   ;---------------------------------------------------------------------------
   ; Do the third traceset fit
   ;---------------------------------------------------------------------------

   ; Fit arc lines subtracting out scatter term
   xmeasured = xnew
   xnew = fitmx(wset, lambda, xmeasured)

   ; In this final fit, do no rejection

   xy2traceset, transpose(xnew), lambda # (dblarr(nfiber)+1), wset, $
    func=func, ncoeff=ncoeff, maxdev=0, maxiter=nlambda, $
    xmin=0, xmax=npix-1

   print, 'Pass 3 complete'

   ;---------------------------------------------------------------------------
   ; Quality Assurance
   ;---------------------------------------------------------------------------

   ; pixel positions derived from the traceset

   tset_pix = transpose( traceset2pix(wset, lambda) )

   xdif_tset = (xmeasured-tset_pix)  ; difference between measured line 
                                     ;  positions and fit positions
   xdif_lfit = (xmeasured-xnew)      ; dif between measured line positions
                                     ;  and best fit for each line

   print
   print, 'Arcline fit summary'
   print, 'All sigma values are in millipixels'
   print, format='(71("-"))'
   for k=0, nlambda-1 do $
      print,'Arcline',k,':  lambda =',10.^lambda[k], $
            '    sig_lfit =', djsig(1e3*xdif_lfit[k,*]), $
            '    sig_tset =', djsig(1e3*xdif_tset[k,*]), $
            format='(A,I3,A,F8.2,A,F7.2,A,F7.2)'

   highones = where(lambda GT alog10(8000.), highct)
   if (nlambda LT 6) then $
     message, 'WARNING: only '+string(nlambda)+ ' good arclines found',/cont

   print, 'Found ', nlambda, ' good arc lines'
   if (highct GT 0) then print, '----', highct, ' are above 8000 A'
   print,'> FITARCIMAGE: ',systime(1)-t_begin, ' seconds elapsed', $
    format='(A,F8.2,A)'

   return
end
;------------------------------------------------------------------------------
