;+
; NAME:
;   spflatten2
;
; PURPOSE:
;   Create pixel-to-pixel flat-field from a stack of SDSS spectral flats.
;
; CALLING SEQUENCE:
;   spflatten2, flatname, arcname, allflats, [ pixflat, sigrej=, maxiter=, $
;    oldflat=, outfile=, indir=, outdir=, tmpdir=, $
;    pixspace=, nord=, lower=, upper= ]
;
; INPUTS:
;   flatname   - Name of flat image for tracing arc
;   arcname    - Name of arc image
;   allflats   - Name(s) of raw SDSS flat-field image(s).
;                Note that many flats from many nights can be combined.
;
; OPTIONAL INPUTS:
;   sigrej     - Sigma rejection level; default to 1, 1, 1.1, 1.3, 1.6 or 1.9
;                for 1,2,3,4,5 or 6 flats.  For more then 6 flats, default
;                to 2.0.
;   maxiter    - Number of rejection iterations; default to 2.
;   oldflat    - Name of old flat-field from which to select pixels to mask;
;                if not set, then read the masks in this source distribution.
;   outfile    - Write the image PIXFLAT to this file.
;   indir      - Input directory for FLATNAME; default to './'
;   outdir     - Output directory for OUTFILE; default to './'
;   tmpdir     - Directory for temporary files; default to same as OUTDIR
;   pixspace   - Approximate spacing in pixels for break points in the
;                spline fits to individual fibers; default to 50 pixels.
;                Note that this spacing need be small enough to fit out
;                the flux variations of multiple fibers crossing a single
;                column, but it need not fit out the actual flat-field
;                variations.
;
; PARAMETERS FOR SLATEC_SPLINEFIT:
;   nord       - Default to 4
;   lower      - Default to 2
;   upper      - Default to 2
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   pixflat    - Image containing all the information about pixel-to-pixel
;                variations.  Illumination variations are removed.
;
; COMMENTS:
;   This program writes 2*nflat temporary files to disk to save internal memory.
;   But it still creates an two arrays (FLATARR (float) and OUTMASK (byte))
;   that is as large as all of the input flats.
;   Thus, if you are using ten 16 MB flat images,
;   that array will be 10*(16+4)=200 MB (in addition to a few other images).
;
; EXAMPLES:
;
; BUGS:
;   Not sure what to do if exactly 2 frames are passed.
;   Pass UPPER and LOWER.
;   Look for bad fibers or outlier fibers to exclude from superflat.
;   Right now, we can get bad spline fits across the large hole in the
;     b1 flat-field if the break points fall too close to the edges of
;     that hole.  Try to keep more data points between the last break point
;     and the beginning of the hole, or do a linear interpolation across
;     the hole.
;
; PROCEDURES CALLED:
;   djs_avsigclip()
;   extract_image
;   genflatmask()
;   readfits()
;   rmfile
;   slatec_bvalu()
;   slatec_splinefit()
;   sdssproc
;   superflat
;   trace320crude()
;   traceset2xy
;   writefits
;   xy2traceset
;
; REVISION HISTORY:
;   13-Oct-1999  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
pro spflatten2, flatname, arcname, allflats, pixflat, $
 sigrej=sigrej, maxiter=maxiter, $
 oldflat=oldflat, outfile=outfile, indir=indir, outdir=outdir, tmpdir=tmpdir, $
 pixspace=pixspace, nord=nord, lower=lower, upper=upper

   if (NOT keyword_set(indir)) then indir = './'
   if (NOT keyword_set(outdir)) then outdir = './'
   if (NOT keyword_set(tmpdir)) then tmpdir = outdir

   if (N_elements(pixspace) EQ 0) then pixspace = 50
   if (N_elements(nord) EQ 0) then nord = 4
   if (N_elements(lower) EQ 0) then lower = 2
   if (N_elements(upper) EQ 0) then upper = 2

   nflat = N_elements(allflats)
   ngrow = 2

   if (NOT keyword_set(sigrej)) then begin
      if (nflat LE 2) then sigrej = 1.0 $ ; Irrelevant for only 1 or 2 flats
       else if (nflat EQ 3) then sigrej = 1.1 $
       else if (nflat EQ 4) then sigrej = 1.3 $
       else if (nflat EQ 5) then sigrej = 1.6 $
       else if (nflat EQ 6) then sigrej = 1.9 $
       else sigrej = 2.0
   endif
   if (NOT keyword_set(maxiter)) then maxiter = 2

   tmpname1 = filepath('tmp.flatimg.'+strtrim(string(indgen(nflat)),2)+'.fits',$
    root_dir=tmpdir)
   tmpname2 = filepath('tmp.ymodel.'+strtrim(string(indgen(nflat)),2)+'.fits',$
    root_dir=tmpdir)

   ;---------------------------------------------------------------------------
   ; First find the wavelength solution
   ;---------------------------------------------------------------------------

   ;------
   ; Read flat-field image that corresponds to the arc

   splog, 'Reading flat ', flatname
   sdssproc, flatname[0], flatimg, flativar, indir=indir, hdr=flathdr, $
    spectrographid=spectrographid, color=color

   dims = size(flatimg, /dimens)
   nx = dims[0]
   ny = dims[1]

   ;------
   ; Create spatial tracing from flat-field image

   splog, 'Tracing 320 fibers in ',  flatname
   xsol = trace320crude(flatimg, yset=ycen, maxdev=0.15)

   splog, 'Fitting traces in ',  flatname
   xy2traceset, ycen, xsol, tset, ncoeff=5, maxdev=0.1
   traceset2xy, tset, ycen, xsol

   ;------
   ; Read the arc

   splog, 'Reading arc ', arcname
   sdssproc, arcname, arcimg, arcivar, indir=indir, hdr=archdr, $
    color=color

   ;------
   ; Extract the arc image

   splog, 'Extracting arc image with simple gaussian'
   sigma = 1.0
   proftype = 1 ; Gaussian
   highrej = 15
   lowrej = 15
   npoly = 1 ; maybe more structure
   wfixed = [1,1] ; Just fit the first gaussian term

   extract_image, arcimg, arcivar, xsol, sigma, flux, fluxivar, $
    proftype=proftype, wfixed=wfixed, $
    highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1
arcimg = 0
arcivar = 0

   ;------
   ; Compute wavelength calibration for arc lamp

   arccoeff = 5

   splog, 'Searching for wavelength solution'
   fitarcimage, flux, fluxivar, xpeak, ypeak, wset, ncoeff=arccoeff, $
    color=color, lampfile=lampfile, bestcorr=corr

   ;---------------------------------------------------------------------------
   ; Construct wavelength image
   ;---------------------------------------------------------------------------

   ;------
   ; Compute wavelength at every pixel on the CCD
   ; WAVEIMG will be in units of log10(lambda)

   waveimg = fltarr(nx,ny)

   traceset2xy, wset, xx, loglam

   xy2traceset, transpose(xsol), transpose(loglam), tmpset, $
    func='legendre', ncoeff=5, xmin=0, xmax=nx-1, maxsig=2.0
   xtmp = 0
   traceset2xy, tmpset, xtmp, waveimg

   ;---------------------------------------------------------------------------
   ; Read or create a mask image
   ;---------------------------------------------------------------------------

   maskimg = genflatmask(oldflat, spectrographid=spectrographid, color=color, $
    indir=indir, tmpdir=tmpdir)

   ; If a column contains fewer than 5 good points, mask that entire column
   for i=0, nx-1 do begin
      junk = where(maskimg[i,*] EQ 0, ngind)
      if (ngind LT 5) then begin
         splog, 'Masking column ', i
         maskimg[i,*] = 1
      endif
   endfor

   junk = where(maskimg NE 0, ct)
   splog, 'Number of pixels in bad pixel mask = ', ct

   ;---------------------------------------------------------------------------
   ; Construct the flat
   ;---------------------------------------------------------------------------

   ; Always select the same break points in log-wavelength for all fibers

   nbkpts = fix(ny / pixspace) + 2
   bkpt = (findgen(nbkpts)-1) * (max(waveimg) - min(waveimg)) / (nbkpts-1) $
    + min(waveimg)

   for iflat=0, nflat-1 do begin

      ;----------------------
      ; Read flat-field image

      sdssproc, allflats[iflat], flatimg, flativar, indir=indir, hdr=flathdr

      ;----------------------
      ; Create spatial tracing from flat-field image

      xsol = trace320crude(flatimg, yset=ycen, maxdev=0.15)

      xy2traceset, ycen, xsol, tset, ncoeff=5, maxdev=0.1
      traceset2xy, tset, ycen, xsol

      ;----------------------
      ; Extract the flat-field vectors
      ; Apply the bad pixel mask, MASKIMG, so that we avoid those bad regions
      ; when computing the supersky vector.

      sigma = 1.0
      proftype = 1 ; Gaussian
      highrej = 15
      lowrej = 15
      npoly = 1  ; just fit flat background to each row
      wfixed = [1,1] ; Just fit the first gaussian term

      extract_image, flatimg, flativar * (1-maskimg), $
       xsol, sigma, flux, fluxivar, $
       proftype=proftype, wfixed=wfixed, $
       highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1

      ;----------------------
      ; Construct the "superflat" vector for this particular frame

      superflat, flux, fluxivar, wset, afullbkpt, acoeff, $
       fibermask=fibermask, minval=0.0, lower=lower, upper=upper
flux = 0
fluxivar = 0
      fitimg = slatec_bvalu(waveimg, afullbkpt, acoeff)

      ;----------------------
      ; Divide by the superflat image

      flatimg = flatimg / fitimg
      flativar = flativar * fitimg^2
fitimg = 0

; Test extraction...
;extract_image, flatimg, flativar, xsol, sigma, tmpflux, tmpivar, $
; proftype=proftype, wfixed=wfixed, $
; highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1
;splot, median(tmpflux[*,0], 11)
;for i=0,15 do soplot, median(tmpflux[*,i*15], 11)

      ;----------------------
      ; Create the array of preliminary flats

      if (iflat EQ 0) then pixflatarr = fltarr(nx,ny,nflat)

      ;----------------------
      ; Determine YMODEL image

      splog, 'Solving for YMODEL'

;splot,transpose(10^waveimg[1000,*]),transpose(flatimg[1000,*])
;for i=0,5 do soplot, transpose(10^waveimg[i*10,*]),transpose(flatimg[i*10,*])

      ymodel = 0.0 * flatimg
      for i=0, nx-1 do begin
         ; Burles counter for column number...
         print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])

         indx = where(maskimg[i,*] EQ 0 AND flativar[i,*] GT 0, ngind)

         ; If fewer than 5 good points, then fit all points in this column.
         ; Also, change the values of FLATIVAR to prevent SLATEC_SPLINEFIT()
         ; from crashing.
         if (ngind LT 5) then begin
            indx = lindgen(ny)
            flativar[i,*] = 1.0
         endif

         ;------
         ; The following spline fit chooses breaks points evenly separated
         ; in log-wavelength.

;         yaxis = waveimg[i,*]
;         fullbkpt = slatec_splinefit(yaxis[indx], flatimg[i,indx], coeff, $
;          invvar=flativar[i,indx], bkpt=bkpt, nord=nord, $
;          lower=lower, upper=upper, maxiter=3)

         ;------
         ; The following spline fit chooses breaks points evenly separated
         ; in row number.

         yaxis = findgen(ny)

         fullbkpt = slatec_splinefit(yaxis[indx], flatimg[i,indx], coeff, $
          invvar=flativar[i,indx], bkspace=pixspace, nord=nord, $
          lower=4, upper=4, maxiter=3)
         ymodel[i,*] = slatec_bvalu(yaxis, fullbkpt, coeff)

      endfor

      pixflatarr[*,*,iflat] = (flatimg > 1) / (ymodel > 1)

      ;----------------------
      ; Write FLATIMG and YMODEL to disk
      if (nflat GT 1) then begin
         writefits, tmpname1[iflat], flatimg
         writefits, tmpname2[iflat], ymodel
      endif
flatimg = 0
ymodel = 0

   endfor

   if (nflat EQ 1) then begin

      pixflat = temporary(pixflatarr)

   endif else begin
      ; Find deviant pixels in each pixflat
      meanimg = djs_avsigclip(pixflatarr, 3, sigrej=sigrej, maxiter=maxiter, $
       outmask=outmask)
meanimg = 0
pixflatarr = 0

      outmask = temporary(1-outmask) ; Change to 0=bad, 1=good

      flatimgsum = 0
      ymodelsum = 0
      for iflat=0, nflat-1 do begin

         junk = where(outmask[*,*,iflat] EQ 0, ct)
         splog, 'Number of pixels masked in ', allflats[iflat], ' = ', ct

         flatimgsum = flatimgsum + outmask[*,*,iflat] * readfits(tmpname1[iflat])
         ymodelsum = ymodelsum + outmask[*,*,iflat] * readfits(tmpname2[iflat])

         rmfile, tmpname1[iflat]
         rmfile, tmpname2[iflat]

      endfor

      pixflat = (flatimgsum > 1) / (ymodelsum > 1)
   endelse

   if (keyword_set(outfile)) then $
    writefits, outdir+outfile, pixflat

   return
end
;------------------------------------------------------------------------------
