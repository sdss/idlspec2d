;+
; NAME:
;   fitarcimage
;
; PURPOSE:
;   Determine wavelength calibration from arclines
;
; CALLING SEQUENCE:
;   fitarcimage, arc, arcivar, xnew, ycen, wset, $
;    color=color, lampfile=lampfile, $
;    func=func, aset=aset, ncoeff=ncoeff, lambda=lambda, $
;    thresh=thresh, row=row, $
;    xdif_lfit=xdif_lfit, xdif_tset=xdif_tset
;
; INPUTS:
;   arc        - Extracted arc spectra with dimensions [NY,NFIBER]
;   arcivar    - Inverse variance of ARC
;
; OPTIONAL KEYWORDS:
;   color      - 'red' or 'blue'; not required if ANS is set
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in the IDL path.
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
;   lampfile   - Modified from input to include full path name of file
;   lambda     - returns alog10(wavelength) of good lamp lines
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
;   Need to tag dead fibers - esp. for sky-subtraction.
;
; INTERNAL PROCEDURES:
;   fitmeanx()
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
; LAMBDA = log10-wavelength

function fitmeanx, wset, lambda, xpos, nord=nord

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

pro fitarcimage, arc, arcivar, xnew, ycen, wset, $
 color=color, lampfile=lampfile, $
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

   t_begin = systime(1)

   ndim = size(arc, /n_dim)
   if (ndim NE 2) then $
    message, 'Expecting 2-D arc image'
   dims = size(arc, /dim)
   npix = dims[0]
   nfiber = dims[1]
   if (total(dims NE size(arcivar, /dim))) then $
    message, 'ARC and ARCIVAR must have same dimensions'

   if (NOT keyword_set(row)) then row = (nfiber-30)/2

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
      lampdefault = 'lamphgcdne.dat'
      lampfilename = (djs_locate_file(lampdefault))[0]
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

      ; Extract one spectrum from the 5 spectra around fiber number ROW
      ; by taking the median value at each wavelength

      spec = djs_median(arc[*,row-2:row+2], 2)

      wset = arcfit_guess( spec, lamps.loglam, lamps.intensity, color=color )

   endif else begin

      wset = aset

   endelse

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

   ; Make use of the errors??? - Seems to just mess things up???
   ; Well... the reason for that is satured lines, which return infinite errors
;   xnew = trace_fweight(arc, xcen, ycen, invvar=arcivar, xerr=xerr)

   ;---------------------------------------------------------------------------
   ; Reject bad (i.e., saturated) lines
   ;---------------------------------------------------------------------------

   ; Reject any arc line with more than 10% of the fibers are bad.
   ; Bad fibers are any with an infinite error (ARCIVAR=0) within 1 pixel
   ; of the central wavelength.  Note that saturated lines should then
   ; show up as bad.

   nmatch = N_elements(xstart) ; Number of lamp lines traced
   qgood = bytarr(nmatch)

   for i=0, nmatch-1 do begin
      xpix = round(xnew[*,i]) ; Nearest X position (wavelength) in all traces
      mivar = fltarr(nfiber) + 1
      for ix=-1, 1 do begin
         mivar = mivar * arcivar[ ((xpix+ix)>0)<(npix-1), * ]
      endfor
      junk = where(mivar EQ 0, nbad)
      fracbad = float(nbad) / nfiber
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

   ; Fit arc lines subtracting out scatter term
   xmeasured = xnew
   xnew = fitmeanx(wset, lamps.loglam, xmeasured)

   ; In this final fit, do no rejection

   xy2traceset, transpose(xnew), lamps.loglam # (dblarr(nfiber)+1), wset, $
    func=func, ncoeff=ncoeff, maxdev=0, maxiter=nlamp, $
    xmin=0, xmax=npix-1

   print, 'Pass 3 complete'

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
            '    sig_lfit =', djsig(1e3*xdif_lfit[k,*]), $
            '    sig_tset =', djsig(1e3*xdif_tset[k,*]), $
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

   return
end
;------------------------------------------------------------------------------
