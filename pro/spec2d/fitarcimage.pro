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
;    xdif_lfit=xdif_lfit, xdif_tset=xdif_tset, bestcorr=bestcorr ]
;
; INPUTS:
;   arc        - Extracted arc spectra with dimensions [NY,NFIBER]
;   arcivar    - Inverse variance of ARC
;
; OPTIONAL KEYWORDS:
;   color      - 'red' or 'blue'; not required if ANS is set
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in $IDLSPEC2D_DIR/etc.
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
;   xnew       - pixel position of lines [nfiber, nlambda]
;   ycen       - fiber number [nfiber, nlambda]
;   wset       - traceset (pix -> lambda)
;
; OPTIONAL OUTPUTS:
;   lampfile   - Modified from input to include full path name of file
;   lambda     - returns alog10(wavelength) of good lamp lines
;   fibermask  - (Modified)
;   xdif_lfit  - fit residual for individual arclines
;   xdif_tset  - fit residual of traceset
;   bestcorr   - Correlation coefficient with simulated arc spectrum
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
;   finalarcfit
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

pro fitarcimage, arc, arcivar, xnew, ycen, wset, $
 color=color, lampfile=lampfile, fibermask=fibermask, xcen=xcen, $
 func=func, aset=aset, ncoeff=ncoeff, lambda=lambda, thresh=thresh, $
 row=row, nmed=nmed, xdif_lfit=xdif_lfit, xdif_tset=xdif_tset, bestcorr=bestcorr

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

      wset = arcfit_guess( spec, lamps.loglam, lamps.intensity, color=color, $
       bestcorr=bestcorr )

      aset = wset

      splog, 'Best correlation = ', bestcorr

   endif else begin

      wset = aset

   endelse

   ; Return from routine if XCEN, YCEN and WSET are not to be returned
   if (N_params() LE 2) then return

   ; Trim lamp list to only those within the wavelength range
   ; and denoted as good in the LAMPS structure.

   xstart = traceset2pix(wset, lamps.loglam)
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
   ; but only 0.3 pixels while tracing

   splog, 'Tracing', N_elements(lamps), ' arc lines'
   xcen = trace_crude(arc, yset=ycen, nave=1, nmed=1, xstart=xstart, $
    ystart=row, maxshifte=0.3d, maxshift0=2.0d)

   ; Iterate the flux-weighted centers
   splog, 'Iterating flux-weighted centers'
   xnew = trace_fweight(arc, xcen, ycen)
   xnew = trace_fweight(arc, xnew, ycen)
   xnew = trace_fweight(arc, xnew, ycen, radius=2.0, invvar=arcivar, xerr=xerr)
   xcen = xnew

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
   xerr = xerr[*,igood]

   ;---------------------------------------------------------------------
   ; Junk all centers with xerr = 999.0
   ;
   xweight = (xerr LT 990)

   ;---------------------------------------------------------------------------
   ; Do the first traceset fit
   ;---------------------------------------------------------------------------

; ??? Let maxdev be a parameter; should be about 3.0d-5 = 20 km/s
maxdev = 3.0d-5
maxsig = 3.0

   nlamp = N_elements(lamps)
   xy2traceset, transpose(double(xnew)), lamps.loglam # (dblarr(nfiber)+1), $
     wset, invvar=transpose(xweight), func=func, ncoeff=ncoeff, $
     maxdev=maxdev, maxiter=nlamp, /singlerej, $
     xmask=xmask, xmin=0, xmax=npix-1, yfit=yfit

   print, 'Pass 1 complete'

   ;---------------------------------------------------------------------------
   ; Do the second traceset fit
   ;---------------------------------------------------------------------------

   ; Keep only "good" lines.
   ; The following logic means that an arc line is rejected if any
   ; bundle has more than 3 bad centers.

   ; do not count bad fibers in fibermask


   badfibers = where(fibermask EQ 0)
   if (badfibers[0] NE -1) then xmask[*,badfibers] = 1

   if (nfiber NE 320) then $
    message, 'Not 320 fibers -- Cannot figure out bundle test'
   testg = reform(xmask, nlamp, 20, 16)   
   gind = where(total(total(testg EQ 0,2) GT 3, 2) EQ 0, nlamp)
   if (nlamp EQ 0) then $
    message, 'No good arcs common to all fiber bundles'



   xnew = xnew[*,gind]
   xweight = xweight[*,gind]
   ycen = ycen[*,gind]
   lamps = lamps[gind]

   fixabove = 2

   arcfibermask = fibermask
   badfweight = where(xweight LE 0)
   if (badfweight[0] NE -1) then arcfibermask[badfweight mod nfiber] = 0

   finalarcfit, xnew, lamps.loglam, wset, ncoeff, fixabove, $
              fibermask=arcfibermask, xweight=xweight, func=func, $
              maxdev=maxdev/2.0, maxiter=nlamp, /singlerej, $
              nsetcoeff=8, maxsig=maxsig

   print, 'Final arcfit complete'

   xmeasured = xnew

   ;---------------------------------------------------------------------------
   ; Quality Assurance
   ;---------------------------------------------------------------------------

   ; pixel positions derived from the traceset

   tset_pix = transpose( traceset2pix(wset, lamps.loglam) )

   xdif_tset = (xmeasured-tset_pix)  ; difference between measured line 
                                     ;  positions and fit positions
   xdif_lfit = (xmeasured-xnew)      ; dif between measured line positions
                                     ;  and best fit for each line

   splog, '', /noname
   splog, 'Arcline fit summary'
   splog, 'All sigma values are in millipixels'
   splog, format='(71("-"))'
   for k=0, nlamp-1 do $
      splog,'Arcline',k,':  lambda =',lamps[k].lambda, $
            '    median_tset =', median(1e3*xdif_tset[*,k]), $
            '    sig_tset =', djsig(1e3*xdif_tset[*,k]), $
            format='(A,I3,A,F8.2,A,F7.2,A,F7.2)'

   highones = where(lamps.lambda GT 8000., highct)
   if (nlamp LT 6) then $
     splog, 'WARNING: only '+string(nlamp)+ ' good arclines found'

   splog, 'Found ', nlamp, ' good arc lines'
   if (highct GT 0) then splog, '----', highct, ' are above 8000 A'
   splog, 'Time ',systime(1)-t_begin, ' seconds elapsed', $
    format='(A,F6.0,A)'

   lampfile = lampfilename
   lambda = lamps.lambda

   ; Replace the "measured" arc line positions (XNEW) with the fit positions
   ; Do this so that the sky-line fitting will use those fit positions for
   ; the arc lines
   xnew = tset_pix

   return
end
;------------------------------------------------------------------------------
