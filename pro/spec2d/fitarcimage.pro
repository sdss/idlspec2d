;+
; NAME:
;   fitarcimage
;
; PURPOSE:
;   Determine wavelength calibration from arclines
;
; CALLING SEQUENCE:
;   fitarcimage, arc, arcivar, xcen, ycen, wset, [wfirst=, $
;    color=color, lampfile=lampfile, fibermask=fibermask, $
;    func=func, aset=aset, ncoeff=ncoeff, lambda=lambda, $
;    thresh=thresh, row=row, nmed=nmed, /gauss, $
;    xdif_tset=xdif_tset, bestcorr=bestcorr ]
;
; INPUTS:
;   arc        - Extracted arc spectra with dimensions [NY,NFIBER]
;   arcivar    - Inverse variance of ARC
;
; OPTIONAL KEYWORDS:
;   color      - 'red' or 'blue'; not required if ANS is set
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in $IDLSPEC2D_DIR/etc.
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
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
;   maxdev     - max deviation in log lambda to allow (default 1.0e-5=7 km/s)
;   gauss      - Use gaussian profile fitting for final centroid fit
;
; OUTPUTS:
;   aset       - (Modified)
;   xcen       - pixel position of lines [nfiber, nlambda]
;   ycen       - fiber number [nfiber, nlambda]
;   wset       - traceset (pix -> lambda)
;
; OPTIONAL OUTPUTS:
;   lampfile   - Modified from input to include full path name of file
;   lambda     - Returns wavelengths of good lamp lines [Angstroms]
;   fibermask  - (Modified)
;   xdif_tset  - Fit residual of lamp lines to fit positions [pixels]
;   bestcorr   - Correlation coefficient with simulated arc spectrum
;   wfirst     - traceset from first iteration on arc fits
;
; COMMENTS:
;   Return from routine after computing BESTCORR if XCEN, YCEN and WSET
;   are not to be returned.
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
;   fibermask_bits()
;   trace_crude()
;   trace_fweight()
;   traceset2pix()
;   traceset2xy()
;   xy2traceset
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/lamphgcdne.dat
;
; REVISION HISTORY:
;   15-Oct-1999  Written by S. Burles, D. Finkbeiner, & D. Schlegel, APO.
;   09-Nov-1999  Major modifications by D. Schlegel, Ringberg.
;   20-Jan-2000  Gone back to very simple procedure: replacement (S. Burles)
;-
;------------------------------------------------------------------------------

pro fitarcimage, arc, arcivar, xcen, ycen, wset, wfirst=wfirst, $
 color=color, lampfile=lampfile, fibermask=fibermask, $
 func=func, aset=aset, ncoeff=ncoeff, lambda=lambda, thresh=thresh, $
 row=row, nmed=nmed, xdif_tset=xdif_tset, bestcorr=bestcorr, $
 gauss=gauss, maxdev=maxdev

   ;---------------------------------------------------------------------------
   ; Preliminary stuff
   ;---------------------------------------------------------------------------

   if (NOT keyword_set(aset)) then begin
      if (color NE 'blue' AND color NE 'red') then $
       message, "SIDE must be set to 'blue' or 'red' if ANS not specified"
   endif
   if (NOT keyword_set(func)) then func = 'legendre'
   if (NOT keyword_set(ans)) then ans = 0
   if NOT keyword_set(maxdev) then maxdev = 1.0d-5

   t_begin = systime(1)

   ndim = size(arc, /n_dim)
   if (ndim NE 2) then $
    message, 'Expecting 2-D arc image'
   dims = size(arc, /dim)
   npix = dims[0]
   nfiber = dims[1]
   if (total(dims NE size(arcivar, /dim))) then $
    message, 'ARC and ARCIVAR must have same dimensions'
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber)

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
      lampdefault = filepath('lamphgcdne.dat', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
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

      ii = where(fibermask EQ 0, ngfiber)
      if (ngfiber EQ 0) then $
       message, 'No unmasked fibers according to FIBERMASK'
      ii = ii[ sort(abs(ii-row)) ]
      ii = ii[0:(nmed<ngfiber)] ; At most NGFIBER good fibers

      spec = djs_median(arc[*,ii], 2)

      aset = arcfit_guess( spec, lamps.loglam, lamps.intensity, color=color, $
       bestcorr=bestcorr )

      splog, 'Best correlation = ', bestcorr
      splog, 'ASET = ', aset.coeff, format='(a,99f9.5)'
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
   ; Trace arc lines on the 2D image
   ;---------------------------------------------------------------------------

   ; Allow for a shift of up to 2 pixels in the initial centers,
   ; but only 0.5 pixels while tracing

   splog, 'Tracing', N_elements(lamps), ' arc lines'
   xcen = trace_crude(arc, yset=ycen, nave=1, nmed=1, xstart=xstart, $
    ystart=row, maxshifte=0.5d, maxshift0=2.0d)

   ; Now fit traceset to trace_crude, this will "interpolate" or bad fibers
   ; as well as extend trace off CCD if need be.

   xy2traceset, ycen, xcen, crudeset, yfit=xcrudefit

   ; Iterate the flux-weighted centers
   ; In the last iteration, use the formal errors in the arc image

   splog, 'Iterating flux-weighted centers'
   x1 = trace_fweight(arc, xcrudefit, ycen)
   x1 = trace_fweight(arc, x1, ycen)
   xmid = trace_fweight(arc, x1, ycen)

   if (keyword_set(gauss)) then begin
     x2  = trace_gweight(arc, xmid, ycen, sigma=1.0, invvar=arcivar, xerr=xerr) 
     x2  = trace_gweight(arc, x2, ycen, sigma=1.0, invvar=arcivar, xerr=xerr) 
     xcen = trace_gweight(arc, x2, ycen, sigma=1.0, invvar=arcivar, xerr=xerr) 
   endif else begin  
     x2  = trace_fweight(arc, xmid, ycen, radius=2.0, invvar=arcivar, xerr=xerr)
     x2  = trace_fweight(arc, x2, ycen, radius=2.0, invvar=arcivar, xerr=xerr)
     xcen = trace_fweight(arc, x2, ycen, radius=2.0, invvar=arcivar, xerr=xerr)
   endelse

   ;  xdiff is used to determine whether a trace has converged
   xdiff = xcen - xmid

   ;---------------------------------------------------------------------------
   ; Reject bad (i.e., saturated) lines
   ;---------------------------------------------------------------------------

   ; Reject any arc line with more than 10% of the "good" fibers have bad arcs.
   ; Bad fibers are any with an infinite error (ARCIVAR=0) within 1 pixel
   ; of the central wavelength.  Note that saturated lines should then
   ; show up as bad.

   nmatch = N_elements(xstart) ; Number of lamp lines traced
   igfiber = where(fibermask EQ 0, ngfiber) ; Number of good fibers
   qgood = bytarr(nmatch)

   for i=0, nmatch-1 do begin
      xpix = round(xcen[*,i]) ; Nearest X position (wavelength) in all traces
      mivar = fltarr(ngfiber) + 1
      for ix=-1, 1 do begin  ; Loop +/- 1 pix in the wavelength direction
         mivar = mivar * arcivar[ (((xpix+ix)>0)<(npix-1))[igfiber], igfiber ]
      endfor
      junk = where(mivar EQ 0, nbad)
      fracbad = float(nbad) / ngfiber
      qgood[i] = fracbad LE 0.10
      if (qgood[i] EQ 0) then $
       splog, 'Discarding trace', i, ',   fraction bad', fracbad

      djs_iterstat, xdiff[*,i], sigma=sigma
      if sigma GT 0.2 then begin
       qgood[i] = 0
       splog, 'Discarding trace', i, ',   Did not converge', sigma
      endif
      
   endfor

   ; Trim linelist

   igood = where(qgood, ngood)
   splog, 'Number of good arc lines: ', ngood
   if (ngood LT 6) then begin
      splog, 'WARNING: Reject arc. Number of arclines  = ', ngood
      wset = 0
      return
   endif else begin
      splog, 'Number of arclines  = ', ngood
   endelse

   xcen = xcen[*,igood]
   ycen = ycen[*,igood]
   lamps = lamps[igood]
   xerr = xerr[*,igood]

   ;---------------------------------------------------------------------
   ; Mask all bad centers

   ; Mask all centers with xerr = 999.0 (from trace_fweight)

   xmask = (xerr LT 990)

   ;---------------------------------------------------------------------------
   ; Do the first traceset fit
   ;---------------------------------------------------------------------------

   nlamp = N_elements(lamps)
   if (nlamp EQ 1) then ytmp = transpose(lamps.loglam * (dblarr(nfiber)+1)) $
    else ytmp = lamps.loglam # (dblarr(nfiber)+1)
   xy2traceset, transpose(double(xcen)), ytmp, $
     wfirst, invvar=transpose(xmask), func=func, ncoeff=ncoeff, $
     maxdev=maxdev, maxiter=nlamp, maxrej=1, /sticky, $
     outmask=xfitmask, xmin=0, xmax=npix-1, yfit=yfit

   print, 'Pass 1 complete'

   xmask = xmask AND transpose(xfitmask)

   ;------------------------------------------------------------------------
   ; Select good lines with the <=3 bad per bundle test
   ;------------------------------------------------------------------------

   ; The following logic means that an arc line is rejected if any bundle
   ; has more than 3 bad centers (not including bad fibers in FIBERMASK).

   badfibers = where(fibermask GT 0)
   testg = transpose(xmask)
   if (badfibers[0] NE -1) then testg[*,badfibers] = 1

   if (nfiber NE 320) then $
    message, 'Not 320 fibers -- Cannot figure out bundle test'
   testg = reform(testg, nlamp, 20, 16)
   gind = where(total(total(testg EQ 0,2) GT 3, 2) EQ 0, nlamp)
   bind = where(total(total(testg EQ 0,2) GT 3, 2) NE 0)
   if (nlamp LT 6) then begin
      splog, 'WARNING: Reject arc. Number of arclines after bundle test = ', nlamp
      wset = 0
      return
   endif else begin
      splog, 'Number of arclines after bundle test = ', nlamp
   endelse

   if bind[0] NE -1 then begin
       splog, ' # Lambda                   Bundles           ', /noname
       for jj=0,n_elements(bind)-1 do $
       splog, bind[jj], lamps[bind[jj]].lambda, (total(testg EQ 0,2))[bind[jj],*], $
           format='(i2,f8.2,16i3)', /noname
   endif

   xcen = xcen[*,gind]
   xmask = xmask[*,gind]
   ycen = ycen[*,gind]
   lamps = lamps[gind]

   ;----------
   ; Look for large gaps before the first arc line, between any consecutive
   ; arc lines, are after the last arc line.

   ; Get the starting and ending wavelengths on the image, but LOGWMIN
   ; and LOGWMAX will be switched if wavelengths are descending; that's
   ; why I sort that list.

   traceset2xy, wfirst, xtmp, ytmp
   logwmin = median(ytmp[0,*])
   logwmax = median(ytmp[npix-1,*])
   logwlist = [logwmin, alog10(lamps.lambda), logwmax]
   logwlist = logwlist[sort(logwlist)]
   logwdiff = logwlist[1:nlamp-1] - logwlist[0:nlamp-2]

   for i=0, N_elements(logwdiff)-1 do begin
      if (logwdiff[i] GT 0.10) then begin
         ; Abort for gaps as large as [4000,5035] or [8000,10071] Angstroms.
         splog, 'WARNING: Reject arc. Big wavelength gap from ', 10^logwlist[i], $
          ' to ', 10^logwlist[i+1], ' Ang (LAMPS NOT WARM?)'
         wset = 0
         return
      endif else if (logwdiff[i] GT 0.04) then begin
         ; Warning for gaps as large as [4000,4386] or [8000,8771] Angstroms.
         splog, 'WARNING: Big wavelength gap from ', 10^logwlist[i], $
          ' to ', 10^logwlist[i+1], ' Ang (LAMPS NOT WARM?)'
      endif
   endfor

   ;---------------------------------------------------------------------------
   ;  Now look to replace pixels masked 
   ;---------------------------------------------------------------------------

   ;--------------------------------------------------------------------
   ;  If more than half the lines are masked, than don't use shift

   xy2traceset, ycen, xcen, meanset, yfit=meanfit, invvar=xmask

   meandiff = xcen - meanfit
   badcount = total(xmask EQ 0, 2)
   maxbad = nlamp/4

   fixthese = where(badcount GT 0 AND badcount LT maxbad, nfix)
   splog, 'Fixing centroids in '+string(nfix)+' fibers'
   splog, fixthese

   replacethese = where(badcount GE maxbad, nreplace) 
   splog, 'Replacing all centroids in '+string(nreplace)+' fibers'
   splog, replacethese

   if (fixthese[0] NE -1) then begin
     xy2traceset, transpose(xcen[fixthese,*]), transpose(meandiff[fixthese,*]),$
      diffset, ncoeff=2, invvar=transpose(xmask[fixthese,*]), yfit=diffit, $
      maxdev=0.1

     ; Just replace masked centroids here, with linear fit above

     badmask = where(xmask[fixthese,*] EQ 0)
     xcentemp = xcen[fixthese,*]
     xcentemp[badmask] = $
             (meanfit[fixthese,*])[badmask] + (transpose(diffit))[badmask]
     xcen[fixthese,*] = xcentemp

   endif

   if (replacethese[0] NE -1) then begin
     fibermask[replacethese] = fibermask[replacethese] OR $
            fibermask_bits('BADARC')

     ; Replace the entire fiber with fit, definitely biases wavelength
     xcen[replacethese,*] = meanfit[replacethese,*]
   endif

   ;--------------------------------------------------------------------------
   ;  Now do final traceset fit
   ;--------------------------------------------------------------------------

   if (nlamp EQ 1) then ytmp = transpose(lamps.loglam * (dblarr(nfiber)+1)) $
    else ytmp = lamps.loglam # (dblarr(nfiber)+1)
   xy2traceset, transpose(double(xcen)), ytmp, $
    wset, func=func, ncoeff=ncoeff, $
    maxiter=nlamp, maxrej=1, /sticky, $
    xmin=0, xmax=npix-1, yfit=yfit

   print, 'Final arcfit complete'

   xmeasured = xcen

   ;---------------------------------------------------------------------------
   ; Quality Assurance
   ;---------------------------------------------------------------------------

   ; pixel positions derived from the traceset

   tset_pix = transpose( traceset2pix(wset, lamps.loglam) )

   xdif_tset = (xmeasured-tset_pix)  ; difference between measured line
                                     ; positions and fit positions

   splog, 'Time ',systime(1)-t_begin, ' seconds elapsed', $
    format='(A,F6.0,A)'

   lampfile = lampfilename ; Return the name of the lampfile used.
   lambda = lamps.lambda

   ; Replace the "measured" arc line positions (XNEW) with the fit positions
   ; Do this so that the sky-line fitting will use those fit positions for
   ; the arc lines
   xcen = tset_pix

   return
end
;------------------------------------------------------------------------------
