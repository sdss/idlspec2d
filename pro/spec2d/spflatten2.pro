;+
; NAME:
;   spflatten2
;
; PURPOSE:
;   Create pixel-to-pixel flat-field from a stack of SDSS spectral flats.
;
; CALLING SEQUENCE:
;   spflatten2, flatname, arcname, allflats, [ pixflat, sigrej=, maxiter=, $
;    outfile=, indir=, outdir=, tmpdir=, $
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
;   outfile    - Write the image PIXFLAT to this file.
;   indir      - Input directory for FLATNAME; default to './'
;   outdir     - Output directory for OUTFILE; default to './'
;   tmpdir     - Directory for temporary files; default to same as OUTDIR
;   pixspace   - Approximate spacing in pixels for break points in the
;                spline fits to individual fibers; default to 50 pixels
;
; PARAMTERS FOR SLATEC_SPLINEFIT:
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
;   UNTESTED!  PROBABLY EVEN WORKS WORSE THAN SPFLATTEN.
;
; PROCEDURES CALLED:
;   djs_avsigclip()
;   extract_image
;   readfits()
;   slatec_bvalu()
;   slatec_splinefit()
;   sdssproc
;   writefits
;
; REVISION HISTORY:
;   13-Oct-1999  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------

pro spflatten2, flatname, arcname, allflats, pixflat, $
 sigrej=sigrej, maxiter=maxiter, $
 outfile=outfile, indir=indir, outdir=outdir, tmpdir=tmpdir, $
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

   tmpname1 = tmpdir+'/tmp.flatimg.'+strtrim(string(indgen(nflat)),2)+'.fits'
   tmpname2 = tmpdir+'/tmp.ymodel.'+strtrim(string(indgen(nflat)),2)+'.fits'

   ;---------------------------------------------------------------------------
   ; First find the wavelength solution
   ;---------------------------------------------------------------------------

   ;------
   ; Read first flat-field image

   splog, 'Reading flat ', flatname
   sdssproc, flatname[0], image, invvar, indir=indir, hdr=flathdr

   ;------
   ; Create spatial tracing from flat-field image

   splog, 'Tracing 320 fibers in ',  flatname
   xsol = trace320crude(image, yset=ycen, maxdev=0.15)

   splog, 'Fitting traces in ',  flatname
   xy2traceset, ycen, xsol, tset, ncoeff=5, maxdev=0.1
   traceset2xy, tset, ycen, xsol

   ;---------------------------------------------------------------------
   ; Read the arc
   ;---------------------------------------------------------------------

   splog, 'Reading arc ', arcname
   sdssproc, arcname, image, invvar, indir=indir, hdr=archdr, $
    color=color

   ;------
   ; Extract the arc image

   splog, 'Extracting arc image with simple gaussian'
   sigma = 1.0
   proftype = 1 ; Gaussian
   highrej = 15
   lowrej = 15
   nPoly = 1 ; maybe more structure
   wfixed = [1,1] ; Just fit the first gaussian term

   extract_image, image, invvar, xsol, sigma, flux, fluxivar, $
    proftype=proftype, wfixed=wfixed, $
    highrej=highrej, lowrej=lowrej, nPoly=nPoly, relative=1

   ;------
   ; Compute wavelength calibration for arc lamp

   arccoeff = 5

   splog, 'Searching for wavelength solution'
   fitarcimage, flux, fluxivar, xpeak, ypeak, wset, ncoeff=arccoeff, $
    color=color, lampfile=lampfile, lambda=lambda, bestcorr=corr

   ;------
   ; Compute wavelength at every pixel on the CCD
   ; WAVEIMG will be in units of log10(lambda)

   waveimg = 0 * image

   traceset2xy, wset, xx, yy

   dims = size(image, /dimens)
   nx = dims[0]
   ny = dims[1]

   xy2traceset, transpose(xsol), transpose(yy), tmpset, $
    func='legendre', ncoeff=5, xmin=0, xmax=ncol-1, maxsig=2.0
   xtmp = 0
   traceset2xy, tmpset, xtmp, waveimg

   ; Always select the same break points in log-wavelength for all fibers
   nbkpts = fix(ny / pixspace) + 2
   bkpt = findgen(nbkpts) * (max(waveimg) - min(waveimg)) / (nbkpts-1) $
    + min(waveimg)


   ;---------------------------------------------------------------------------
   ; Construct the flat
   ;---------------------------------------------------------------------------

   for iflat=0, nflat-1 do begin

      ;----------------------
      ; Read flat-field image

      sdssproc, allflats[iflat], flatimg, flativar, indir=indir, hdr=flathdr

      ;----------------------
      ; Create the array of preliminary flats

      if (iflat EQ 0) then pixflatarr = fltarr(nx,ny,nflat)

      ;----------------------
      ; Extract the flat-field image

      sigma = 1.0
      proftype = 1 ; Gaussian
      highrej = 10
      lowrej = 15
      nPoly = 4

      ; Determine YMODEL image
print, 'Working on file ', fullname[0]

      ymodel = 0.0 * flatimg
      for i=0, nx-1 do begin
         print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])
         yaxis = waveimg[i,*]
         fullbkpt = slatec_splinefit(yaxis, flatimg[i,*], coeff, $
          invvar=flativar[i,*], bkpt=bkpt, nord=nord, $
          lower=lower, upper=upper, maxiter=3)
         ymodel[i,*] = slatec_bvalu(yaxis, fullbkpt, coeff)
;plot,10^waveimg[i,*],flatimg[i,*]
;djs_oplot,10^waveimg[i,*],ymodel[i,*],color='red'
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
