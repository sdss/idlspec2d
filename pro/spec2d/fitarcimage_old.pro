;+
; NAME:
;   fitarcimage
;
; PURPOSE:
;   Determine wavelength calibration from arclines
;
; CALLING SEQUENCE:
;   fitarcimage, arc, arcivar, xnew, ycen, wset, $
;    [ color=color, lampfile=lampfile, fibermask=fibermask, $
;    func=func, aset=aset, ncoeff=ncoeff, lambda=lambda, $
;    thresh=thresh, row=row, nmed=nmed, $
;    xdif_tset=xdif_tset, bestcorr=bestcorr ]
;
; INPUTS:
;   arc        - Extracted arc spectra with dimensions [NY,NFIBER]
;   arcivar    - Inverse variance of ARC
;
; OPTIONAL KEYWORDS:
;   color      - 'red' or 'blue'; not required if ANS is set
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in the IDL path.
;   fibermask  - Mask of 0 for bad fibers and 1 for good fibers [NFIBER]
;   func       - Name of fitting function; default to 'legendre'
;   aset       - Trace set for initial wavelength solution in row number ROW.
;   ncoeff     - Number of coefficients in fits.  This may be different than
;                the number of coefficients in the initial guess ASET.
;                Default to 5.
;   thresh     - Threshhold counts for significant lines;
;                default to 200 if COLOR='blue' or 500 if COLOR='red'
;   row        - Row to use in initial guess of wavelength solution;
;                default to (NFIBER-30)/2
;   nmed       - Number of rows around ROW to median filter for initial
;                wavelengths solution; default to 5
;
; OUTPUTS:
;   aset       - (Modified)
;   xnew       - pixel position of lines [nfiber, nlambda]
;   ycen       - fiber number [nfiber, nlambda]
;   wset       - traceset (pix -> lambda)
;
; OPTIONAL OUTPUTS:
;   lampfile   - Modified from input to include full path name of file
;   lambda     - Returns wavelengths of good lamp lines [Angstroms]
;   fibermask  - (Modified)
;   xdif_tset  - Fit residual of lamp lines to fit positions [pixels]
;   bestcorr   - Correlation coefficient with simulated arc spectrum
;
; COMMENTS:
;   Return from routine after computing BESTCORR if XCEN, YCEN and WSET
;   are not to be returned.
;
;   THIS IS REVISION 1.27 OF FITARCIMAGE, SINCE THE CURRENT VERSION
;   APPEARS TO BE BUGGY!!!
;
; EXAMPLES:
;
; BUGS:
;   Not making sure that only the same lines are fit for each fiber.
;      (Different lines can be rejected in xy2traceset.)
;   THRESH is unused.
;   TRACESET2PIX maybe returns the transpose of what is natural?
;   Check QA stuff at end.
;   FIBERMASK not yet modified if an arc is atrociously bad.
;
;
; PROCEDURES CALLED:
;   arcfit_guess()
;   djs_median
;   djsig()
;   trace_crude()
;   trace_fweight()
;   traceset2pix()
;   traceset2xy()
;   xy2traceset
;
; REVISION HISTORY:
;   15-Oct-1999  Written by S. Burles, D. Finkbeiner, & D. Schlegel, APO.
;   09-Nov-1999  Major modifications by D. Schlegel, Ringberg.
;-
;------------------------------------------------------------------------------

pro fitarcimage_old, arc, arcivar, xnew, ycen, wset, $
 color=color, lampfile=lampfile, fibermask=fibermask, $
 func=func, aset=aset, ncoeff=ncoeff, lambda=lambda, thresh=thresh, $
 row=row, nmed=nmed, xdif_tset=xdif_tset, bestcorr=bestcorr

   ;---------------------------------------------------------------------------
   ; Preliminary stuff
   ;---------------------------------------------------------------------------

   if (NOT keyword_set(aset)) then begin
      if (color NE 'blue' AND color NE 'red') then $
       message, "SIDE must be set to 'blue' or 'red' if ANS not specified"
   endif
   if (NOT keyword_set(func)) then func = 'legendre'
   if (NOT keyword_set(ans)) then ans = 0

   t_begin = systime(1)

   ndim = size(arc, /n_dim)
   if (ndim NE 2) then $
    message, 'Expecting 2-D arc image'
   dims = size(arc, /dim)
   npix = dims[0]
   nfiber = dims[1]
   if (total(dims NE size(arcivar, /dim))) then $
    message, 'ARC and ARCIVAR must have same dimensions'
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber) + 1

   if (NOT keyword_set(row)) then row = (nfiber-30)/2
   if (NOT keyword_set(nmed)) then nmed = 5

   if (NOT keyword_set(thresh)) then begin
      if (color EQ 'blue') then thresh = 200
      if (color EQ 'red') then thresh = 500
   endif

   if (NOT keyword_set(ncoeff)) then ncoeff = 5

   ;---------------------------------------------------------------------------
   ; Read LAMPLIST file for wavelength calibration
   ;---------------------------------------------------------------------------
   ; Read this file into a structure

   if (keyword_set(lampfile)) then begin
      lampfilename = (findfile(lampfile, count=ct))[0]
      if (ct EQ 0) then message, 'No LAMPFILE found '+lampfile
   endif else begin
      lampdefault = getenv('IDLSPEC2D_DIR') + '/etc/lamphgcdne.dat'
      lampfilename = (findfile(lampdefault, count=ct))[0]
      if (NOT keyword_set(lampfilename)) then $
       message, 'No LAMPFILE found '+lampdefault
   endelse

   splog, 'Reading lamp file ', lampfilename
   readcol, lampfilename, lampwave, lampinten, lampquality, format='D,F,A'
   lamps = {lambda: 0.0d0, loglam: 0.0d0, intensity: 0.0d0, good: 0.0d0}
   lamps = replicate(lamps, N_elements(lampwave))
   lamps.lambda = lampwave
   lamps.loglam = alog10(lampwave)
   lamps.intensity = lampinten
   lamps.good = strupcase(lampquality) EQ 'GOOD' AND lampinten GT 0

   ;---------------------------------------------------------------------------
   ; INITIAL WAVELENGTH SOLUTION
   ;---------------------------------------------------------------------------

   ; Find the initial wavelength solution if ANS is not passed
   ; One might want to change nave and nmed for first pass???

   if (NOT keyword_set(aset)) then begin
      ; Extract one spectrum from the NMED spectra around fiber number ROW
      ; by taking the median value at each wavelength.
      ; Find the NMED fibers nearest to ROW that are not masked.

      ii = where(fibermask, ngfiber)
      if (ngfiber EQ 0) then $
       message, 'No unmasked fibers according to FIBERMASK'
      ii = ii[ sort(abs(ii-row)) ]
      ii = ii[0:(nmed<ngfiber)] ; At most NGFIBER good fibers

      spec = djs_median(arc[*,ii], 2)

      aset = arcfit_guess( spec, lamps.loglam, lamps.intensity, color=color, $
       bestcorr=bestcorr )

      splog, 'Best correlation = ', bestcorr
   endif

   ; Return from routine if XCEN, YCEN and WSET are not to be returned
   if (N_params() LE 2) then return

   ; Trim lamp list to only those within the wavelength range
   ; and denoted as good in the LAMPS structure.

   xstart = traceset2pix(aset, lamps.loglam)
   qtrim = xstart GT 1 AND xstart LT npix-2 AND lamps.good
   itrim = where(qtrim, ct)
   if (ct EQ 0) then $
    message, 'No arc lines in wavelength range'
   xstart = xstart[itrim]
   lamps = lamps[itrim]

   ;---------------------------------------------------------------------------
   ; Trace
   ;---------------------------------------------------------------------------

   ; Allow for a shift of up to 2 pixels in the initial centers,
   ; but only 0.5 pixels while tracing

   splog, 'Tracing', N_elements(lamps), ' arc lines'
   xcen = trace_crude(arc, yset=ycen, nave=1, nmed=1, xstart=xstart, $
    ystart=row, maxshifte=0.5d, maxshift0=2.0d)

   ; Iterate the flux-weighted centers
   splog, 'Iterating with gaussian-weighted centers'
   xnew = trace_gweight(arc, xcen, ycen, sigma = 1.0, invvar=arcivar, xerr=xerr)
   stop

   ; Make use of the errors??? - Seems to just mess things up???
   ; Well... the reason for that is satured lines, which return infinite errors
;   xnew = trace_fweight(arc, xcen, ycen, invvar=arcivar, xerr=xerr)

   ;---------------------------------------------------------------------------
   ; Reject bad (i.e., saturated) lines
   ;---------------------------------------------------------------------------

   ; Reject any arc line with more than 10% of the "good" fibers have bad arcs.
   ; Bad fibers are any with an infinite error (ARCIVAR=0) within 1 pixel
   ; of the central wavelength.  Note that saturated lines should then
   ; show up as bad.

   nmatch = N_elements(xstart) ; Number of lamp lines traced
   igfiber = where(fibermask, ngfiber) ; Number of good fibers
   qgood = bytarr(nmatch)

   for i=0, nmatch-1 do begin
      xpix = round(xnew[*,i]) ; Nearest X position (wavelength) in all traces
      mivar = fltarr(ngfiber) + 1
      for ix=-1, 1 do begin
         mivar = mivar * arcivar[ (((xpix+ix)>0)<(npix-1))[igfiber], igfiber ]
      endfor
      junk = where(mivar EQ 0, nbad)
      fracbad = float(nbad) / ngfiber
      qgood[i] = fracbad LE 0.10
      if (qgood[i] EQ 0) then $
       splog, 'Discarding trace', i, ',   fraction bad', fracbad
   endfor

   ; Trim linelist
   igood = where(qgood, ngood)
   splog, 'Number of good arc lines: ', ngood
   if (ngood EQ 0) then $
    message, 'No good arc lines'
   xnew = xnew[*,igood]
   ycen = ycen[*,igood]
   lamps = lamps[igood]

   ;---------------------------------------------------------------------------
   ; Do the first traceset fit
   ;---------------------------------------------------------------------------

; ??? Let maxdev be a parameter; should be about 3.0d-5 = 20 km/s
maxdev = 3.0d-5

   nlamp = N_elements(lamps)
   xy2traceset, transpose(double(xnew)), lamps.loglam # (dblarr(nfiber)+1), $
    wset, func=func, ncoeff=ncoeff, maxdev=maxdev, maxiter=nlamp, /singlerej, $
    xmask=xmask, xmin=0, xmax=npix-1

   print, 'Pass 1 complete'

   ;---------------------------------------------------------------------------
   ; Do the second traceset fit
   ;---------------------------------------------------------------------------

   ; Keep only "good" lines.
   ; The following logic means that an arc line is rejected if any
   ; bundle has more than 3 bad centers.

   if (nfiber NE 320) then $
    message, 'Not 320 fibers -- Cannot figure out bundle test'
   testg = reform(xmask, nlamp, 20, 16)   
   gind = where(total(total(testg EQ 0,2) GT 3, 2) EQ 0, nlamp)
   if (nlamp EQ 0) then $
    message, 'No good arcs common to all fiber bundles'
   xnew = xnew[*,gind]
   ycen = ycen[*,gind]
   lamps = lamps[gind]

   xy2traceset, transpose(xnew), lamps.loglam # (dblarr(nfiber)+1), wset, $
    func=func, ncoeff=ncoeff, maxdev=maxdev, maxiter=nlamp, /singlerej, $
    xmask=xmask, xmin=0, xmax=npix-1
   print, 'Pass 2 complete'

   ;---------------------------------------------------------------------------
   ; Do the third traceset fit
   ;---------------------------------------------------------------------------

wsave = wset

   ;------
   ; Fit and replace all the coefficients numbered 1 through NCOEFF-1
   ; with their fit value.  Fit a Chebyshev with 8 coefficients, and
   ; a split in the baseline at the central fibers.
   for ic=1, ncoeff-1 do begin
      xy2traceset, dindgen(nfiber), transpose(wset.coeff[ic,*]), tmpset, $
       func='chebyshev', ncoeff=8, maxsig=2.0, xmask=xmask, yfit=yfit, $
       /halfintwo
      wset.coeff[ic,*] = yfit
; jj=where(xmask EQ 0)
; plot,transpose(wset.coeff[ic,*]),/yno
; djs_oplot,yfit,color='green'
; djs_oplot,jj,(transpose(wset.coeff[ic,*]))[jj],color='red',ps=1
   endfor

   ;------
   ; Re-fit the 0th coefficient for each fiber, e.g. the 0-pt shift term
   ; THIS CODE SHOULD BE SPAWNED OFF TO A ROUTINE FUNC_EVAL() ???

   ifit = [0,1] ; List of the coefficients to fit for each fiber
   ifix = indgen(ncoeff-2) + 2 ; List of coefficients to keep fixed
   ncfit = N_elements(ifit)
   tset = wset
   tset.coeff[ifit,*] = 0
   dims = size(tset.coeff, /dim)
   ncoeff = dims[0]
   nTrace = dims[1]
   xmid = 0.5 * (tset.xmin + tset.xmax)
   xrange = tset.xmax - tset.xmin
   xpos = transpose(xnew)
   ypos = fltarr(nlamp, ntrace)
   for i=0, ntrace-1 do begin
      xvec = 2.0 * xpos[*,i]/(npix-1) - 1.0
      if (tset.func EQ 'legendre') then legarr = flegendre(xvec, ncoeff)
      if (tset.func EQ 'chebyshev') then legarr = fchebyshev(xvec, ncoeff)
      ypos[*,i] = legarr # tset.coeff[*,i]
   endfor

   ;-----
   ; Fit the first NCFIT coefficients to the residuals, with same rejection
   ; as before

   xy2traceset, transpose(xnew), lamps.loglam # (dblarr(nfiber)+1) - ypos, $
    tset2, func=func, ncoeff=ncfit, maxdev=maxdev, maxiter=nlamp, /singlerej, $
    xmask=xmask, xmin=0, xmax=npix-1

   ;-----
   ; Now replace WSET with the re-fit coefficients

   wset.coeff[ifit,*] = tset2.coeff
   wset.coeff[ifix,*] = tset.coeff[ifix,*]

   ;-----
   ; Comment out this FITMEANX code...

   xmeasured = xnew

;   ; Fit arc lines subtracting out scatter term
;   xnew = fitmeanx(wset, lamps.loglam, xmeasured)
;
;   ; In this final fit, do no rejection
;
;   xy2traceset, transpose(xnew), lamps.loglam # (dblarr(nfiber)+1), wset, $
;    func=func, ncoeff=ncoeff, maxiter=0, $
;    xmin=0, xmax=npix-1

   print, 'Pass 3 complete'

   ;---------------------------------------------------------------------------
   ; Quality Assurance
   ;---------------------------------------------------------------------------

   ; pixel positions derived from the traceset

   tset_pix = transpose( traceset2pix(wset, lamps.loglam) )

   xdif_tset = (xmeasured-tset_pix)  ; difference between measured line
                                     ; positions and fit positions

   splog, 'Time ',systime(1)-t_begin, ' seconds elapsed', $
    format='(A,F6.0,A)'

   lampfile = lampfilename ; Return the name of the lamp file used.
   lambda = lamps.lambda

   ; Replace the "measured" arc line positions (XNEW) with the fit positions
   ; Do this so that the sky-line fitting will use those fit positions for
   ; the arc lines
   xnew = tset_pix

   return
end
;------------------------------------------------------------------------------
