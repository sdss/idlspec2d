;+
; NAME:
;   fitarcimage
;
; PURPOSE:
;   Determine wavelength calibration from arclines
;
; CALLING SEQUENCE:
;   fitarcimage, arc, arcivar, xcen, ycen, wset, [wfirst=, $
;    color=, lampfile=, fibermask=, func=, aset=, ncoeff=, $
;    thresh=, row=, nmed=, maxdev=, /gauss, wrange=, $
;    lambda=, rejline=, /twophase, xdif_tset=, $
;    bestcorr=, _EXTRA=KeywordsForArcfit_guess
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
;                default to NFIBER/2
;   nmed       - Number of rows around ROW to median filter for initial
;                wavelengths solution; default to 5
;   maxdev     - max deviation in log lambda to allow (default 2d-5=13.8 km/s)
;   gauss      - Use gaussian profile fitting for final centroid fit
;   wrange     - Wavelength range [Ang vacuum] for searching for arc lines and
;                for issuing warnings about big wavelength gaps in
;                solution
;   twophase   - Set this keyword for BOSS red-side 2-phase readout format
;   _EXTRA     - Keywords for ARCFIT_GUESS, specifically ACOEFF,DCOEFF
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
;   rejline    - String set to non-zero for any rejected arc lines
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
;   25-Jan-2011  Added "twophase" keyword, A. Bolton, U. of Utah
;-
;------------------------------------------------------------------------------
function fitarc_maskbadfib, wset, outmask
   maxline = max( total(outmask,1) ) ; max number of good arc lines in 1 fiber
   ibadfib = where( total(outmask,1) LT 0.70 * maxline, nbadfib)
   if (nbadfib GT 0) then begin
      for i=0L, nbadfib-1L do $
       splog, 'Replacing all centroids in fiber ', ibadfib[i]
      imask = lonarr(size(wset.coeff,/dimens))
      imask[*,ibadfib] = 1
      wset.coeff = djs_maskinterp(wset.coeff, imask, iaxis=1, /const)
   endif

   return, wset
end
;------------------------------------------------------------------------------
pro fitarcimage, arc, arcivar, xcen, ycen, wset, wfirst=wfirst, $
 color=color, lampfile=lampfile, fibermask=fibermask, $
 func=func, aset=aset, ncoeff=ncoeff, thresh=thresh, $
 row=row, nmed=nmed, maxdev=maxdev, gauss=gauss, wrange=wrange, $
 lambda=lambda, rejline=rejline, twophase=twophase, $
 xdif_tset=xdif_tset, bestcorr=bestcorr, _EXTRA=KeywordsForArcfit_guess

   ;---------------------------------------------------------------------------
   ; Preliminary stuff
   ;---------------------------------------------------------------------------

   if (NOT keyword_set(aset)) then begin
      if (color NE 'blue' AND color NE 'red') then $
       message, "SIDE must be set to 'blue' or 'red' if ANS not specified"
   endif
   if (NOT keyword_set(func)) then func = 'legendre'
   if (NOT keyword_set(ans)) then ans = 0
   if NOT keyword_set(maxdev) then maxdev = 2.0d-5
   minlamp = 6 ; Minimum number of lamp lines

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

   if (NOT arg_present(row)) then row = nfiber/2
   if (NOT keyword_set(nmed)) then nmed = 5

   if (NOT keyword_set(thresh)) then begin
      if (color EQ 'blue') then thresh = 200
      if (color EQ 'red') then thresh = 500
   endif

   if (NOT keyword_set(ncoeff)) then ncoeff = 5
   splog, 'Setting ncoeff= ', ncoeff

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

      if nmed GT 1 then begin
         ii = ii[ sort(abs(ii-row)) ]
         ii = ii[0:(nmed<ngfiber)] ; At most NGFIBER good fibers
         spec = djs_median(arc[*,ii], 2)
      end else begin
         spec = arc[*,row]
      end

      aset = arcfit_guess( spec, lamps.loglam, lamps.intensity, color=color, $
       bestcorr=bestcorr, _EXTRA=KeywordsForArcfit_guess)

      splog, 'Best correlation = ', bestcorr
      splog, 'ASET = ', aset.coeff, format='(a,99f9.5)'
   endif

   ; Return from routine if XCEN, YCEN and WSET are not to be returned
   if (N_params() LE 2) then return

   ; Trim lamp list to only those within the wavelength range
   ; and denoted as good in the LAMPS structure.

   xstart = traceset2pix(aset, lamps.loglam)
   qtrim = xstart GT 0 AND xstart LT npix-2 AND lamps.good
   if (keyword_set(wrange)) then $
    qtrim = qtrim AND lamps.loglam GT alog10(wrange[0]) $
     AND lamps.loglam LT alog10(wrange[1])
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

   nlamp = n_elements(lamps)
   splog, 'Tracing', nlamp, ' arc lines'
   xcen = trace_crude(arc, arcivar, yset=ycen, nave=1, nmed=1, xstart=xstart, $
    ystart=row, maxshifte=0.5d, maxshift0=2.0d)

   ; Now fit traceset to trace_crude, this will "interpolate" or bad fibers
   ; as well as extend trace off CCD if need be.

   xy2traceset, ycen, xcen, crudeset, yfit=xcrudefit, ncoeff=ncoeff

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

   ;----------
   ; Reject bad lines that are saturated or very noisy

   ; Reject any arc line with more than 10% of the "good" fibers have bad arcs.
   ; Bad fibers are any with an infinite error (ARCIVAR=0) within 1 pixel
   ; of the central wavelength.  Note that saturated lines should then
   ; show up as bad.

   nmatch = N_elements(xstart) ; Number of lamp lines traced
   igfiber = where(fibermask EQ 0, ngfiber) ; Number of good fibers
   rejline = strarr(nlamp)

   for i=0, nmatch-1 do begin
      xpix = round(xcen[*,i]) ; Nearest X position (wavelength) in all traces
      mivar = fltarr(ngfiber) + 1
      for ix=-1, 1 do begin  ; Loop +/- 1 pix in the wavelength direction
         mivar = mivar * arcivar[ (((xpix+ix)>0)<(npix-1))[igfiber], igfiber ]
      endfor
      junk = where(mivar EQ 0, nbad)
      fracbad = float(nbad) / ngfiber
      if (fracbad GT 0.10) then begin
         rejline[i] = 'Reject-fracbad'
         splog, 'Discarding line', i, ' at ', lamps[i].lambda, $
          ' Ang: FRACBAD=', fracbad
      endif

      djs_iterstat, xdiff[*,i], sigma=sigma
      if (sigma GT 0.2) then begin
         rejline[i] = 'Reject-RMS'
         splog, 'Discarding line', i, ', at ', lamps[i].lambda, $
          ' Ang: RMS=', sigma, ' pix'
      endif
   endfor

   igood = where(rejline EQ '', ngood)
   splog, 'Number of arc lines: ', ngood
   if (ngood LT minlamp) then begin
      splog, 'WARNING: Reject arc image too few lines'
      wset = 0
      return
   endif

   ; Trim linelist
;   xcen = xcen[*,igood]
;   ycen = ycen[*,igood]
;   lamps = lamps[igood]
;   xerr = xerr[*,igood]

   ;----------
   ; Mask all bad centers

   ; Mask all centers with xerr = 999.0 (from trace_fweight)

   inmask = (xerr LT 990)

   if (nlamp EQ 1) then ytmp = transpose(lamps.loglam * (dblarr(nfiber)+1)) $
    else ytmp = lamps.loglam # (dblarr(nfiber)+1)

   ;------------------------------------------------------------------------
   ; Iterate the fit, rejecting the worst-fit line each time

; Include 2-phase readout arguments if necessary:
   if keyword_set(twophase) then begin
      xjumplo = npix/2 - 1 - 0.5
      xjumphi = npix/2 + 0.5
      xjumpval = 1./3.
   endif else begin
      xjumplo = 0
      xjumphi = 0
      xjumpval = 0
   endelse

   iiter = 0
   while (iiter LT nlamp) do begin
      splog, 'Iteration ', iiter
      xy2traceset, transpose(double(xcen)), ytmp, wset, $
       func=func, ncoeff=ncoeff, maxiter=nlamp, maxrej=nlamp, maxdev=maxdev, $
       xmin=0, xmax=npix-1, $
       inmask=transpose(inmask)*rebin(rejline EQ '',nlamp,nfiber), $
       outmask=outmask, yfit=yfit, $
       xjumplo=xjumplo, xjumphi=xjumphi, xjumpval=xjumpval
      ydiff = ytmp - yfit

      if (total(outmask) EQ 0) then begin
         splog, 'WARNING: All centroids rejected'
         wset = 0
         return
      endif

      ; For each lamp, compute median offset, RMS, and slope
      offsets = fltarr(nlamp)
      badper = lonarr(nlamp)
      for i=0, nlamp-1 do begin
         indx = where(outmask[i,*], ct)
         if (ct GT 2) then begin
            djs_iterstat, ydiff[i,indx], median=med1, sigma=sig1
            res1 = polyfitw(findgen(nfiber), ydiff[i,*], outmask[i,*], $
             2, yfit1)
            offsets[i] = sig1 > abs(med1) > abs(max(yfit1))
            ; For each lamp, compute max # of bad fibers per bundle
            badper[i] = max(total(reform(1L-outmask[i,*], 20, nfiber/20),1))
         endif else begin
            yfit1 = 0
            rejline[i] = 'Reject-offset'
            badper[i] = 20
         endelse
      endfor

      ; Preferentially mask any lines with any bundles of more than 5 bad,
      ; then after that mask any lines with bad deviations
      thismax = max(offsets*(rejline EQ '') $
       + ((badper-5)>0)*(rejline EQ ''), iworst)
      if (thismax GT maxdev) then begin
         splog, 'Rejecting line at ', 10d0^lamps[iworst].loglam, ' Ang ' $
          + rejline[iworst]
         rejline[iworst] = 'Reject'
         iiter++
      endif else begin
         splog, 'Iterations done'
         iiter = nlamp
      endelse
   endwhile

   ; Replace all centroids in any bad fibers
; Should this be in the loop above ???
   wset = fitarc_maskbadfib(wset, outmask)

   if (total(rejline EQ '') LT minlamp) then begin
      splog, 'WARNING: Reject arc image too few lines'
      wset = 0
      return
   endif

   ;----------
   ; Replace the higher-order coeff terms with polynomial (cubic) fits of those
   ; coefficients over fiber numbers
   ; The first 3 terms are fit fiber-by-fiber, the others other smoothed.

   ; It gets really ugly here in order to cope with the internal traceset
   ; coeff hacks, including the possibility of the 2-phase discontinuity.

   ; First we make a version of "xcen" with the jump included, if necessary:
   if tag_exist(wset, 'XJUMPVAL') then begin
      jfrac = (((xcen - wset.xjumplo) / (wset.xjumphi - wset.xjumplo)) > 0.) < 1.
      xncen = xcen + jfrac * wset.xjumpval
   endif else xncen = xcen

   if (ncoeff GT 3) then csmooth = lindgen(ncoeff-3) + 3
   if (keyword_set(csmooth)) then begin
      ; get rid of any jump-related tags, if present:
      wset_tmp = struct_selecttags(wset, select_tags=['FUNC', 'XMIN', 'XMAX', 'COEFF'])
      for i=0, n_elements(csmooth)-1 do begin
         poly_iter, lindgen(nfiber), wset_tmp.coeff[csmooth[i],*], 3, 3.0, $
          ycoeff1
         wset_tmp.coeff[csmooth[i],*] = ycoeff1
      endfor
      wset_tmp2 = wset_tmp
      wset_tmp2.coeff[csmooth,*] = 0
      ; Compute the offsets from only the (smoothed) higher-order terms
      xoffset = traceset2pix(wset_tmp, alog10(lamps.lambda)) $
       - traceset2pix(wset_tmp2, alog10(lamps.lambda))
      ; Re-fit, subtracting off those higher-order terms first
      xy2traceset, transpose(double(xncen)) - xoffset, ytmp, wset_lo, $
       func=func, ncoeff=3, maxiter=0, xmin=0, xmax=npix-1, $
       inmask=outmask, yfit=yfit
      wset.coeff[3:ncoeff-1,*] = wset_tmp.coeff[3:ncoeff-1,*]
      wset.coeff[0:2,*] = wset_lo.coeff[0:2,*]
   endif

   for i=0, nlamp-1 do $
    splog, 'Line ', i, lamps[i].lambda, offsets[i], ' '+rejline[i]

   ; Replace all centroids in any bad fibers
   wset = fitarc_maskbadfib(wset, outmask)

   ;----------
   ; Look for large gaps before the first arc line, between any consecutive
   ; arc lines, are after the last arc line.

   ; Get the starting and ending wavelengths on the image, but LOGWMIN
   ; and LOGWMAX will be switched if wavelengths are descending; that's
   ; why I sort that list.

   traceset2xy, wset, xx, yy
   indx = where(arcivar GT 0)
   logwmin = min(yy[indx], max=logwmax)
   if (keyword_set(wrange)) then begin
      logwmin = logwmin > alog10(wrange[0])
      logwmax = logwmax < alog10(wrange[1])
   endif
   logwlist = [logwmin, alog10(lamps[where(rejline EQ '')].lambda), logwmax]
   logwlist = logwlist[sort(logwlist)]
   n = n_elements(logwlist)
   logwdiff = logwlist[1:n-1] - logwlist[0:n-2] ; Allow for endpoints.

   for i=0, N_elements(logwdiff)-1 do begin
      if (logwdiff[i] GT 0.10) then begin
         ; Abort for gaps as large as [4000,5035] or [8000,10071] Angstroms.
         splog, 'WARNING: Reject arc. Big wavelength gap from ', 10^logwlist[i], $
          ' to ', 10^logwlist[i+1], ' Ang'
         wset = 0
         return
      endif else if (logwdiff[i] GT 0.05) then begin
         ; Warning for gaps as large as [3650,4095] or [8000,8976] Angstroms.
         splog, 'WARNING: Big wavelength gap from ', 10^logwlist[i], $
          ' to ', 10^logwlist[i+1], ' Ang'
      endif
   endfor

   ; Specifically test for the Cd 3610 line (ticket #641)
   if (color EQ 'blue') then begin
      junk = where( abs(lamps[where(rejline EQ '')].lambda - 3610.5) LT 2, ct)
      if (ct EQ 0) then $
       splog, 'WARNING: Cd I 3610 line missing (lamps not warm?)'
   endif

   ;---------------------------------------------------------------------------
   ; Quality Assurance
   ;---------------------------------------------------------------------------

   ; pixel positions derived from the traceset

   tset_pix = transpose( traceset2pix(wset, lamps.loglam) )

   xdif_tset = xcen - tset_pix  ; difference between measured line
                                ; positions and fit positions

   splog, 'Time ',systime(1)-t_begin, ' seconds elapsed', $
    format='(A,F6.0,A)'

   lampfile = lampfilename ; Return the name of the lampfile used.
   lambda = lamps.lambda

   ; Replace the "measured" arc line positions with the fit positions
   ; Do this so that the sky-line fitting will use those fit positions for
   ; the arc lines
   xcen = tset_pix

   return
end
;------------------------------------------------------------------------------
