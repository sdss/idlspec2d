;+
; NAME:
;   spflatten
;
; PURPOSE:
;   Create pixel-to-pixel flat-field from a stack of SDSS spectral flats.
;
; CALLING SEQUENCE:
;   spflatten, flatname, [ pixflat, sigrej=, maxiter=, $
;    outfile=, indir=, outdir=, tmpdir=, $
;    bkspace=, nord=, lower=, upper= ]
;
; INPUTS:
;   flatname   - Name(s) of raw SDSS flat-field image(s).
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
;
; PARAMETERS FOR SLATEC_SPLINEFIT:
;   bkspace
;   nord
;   lower
;   upper
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

pro spflatten, flatname, pixflat, sigrej=sigrej, maxiter=maxiter, $
 outfile=outfile, indir=indir, outdir=outdir, tmpdir=tmpdir, $
 bkspace=bkspace, nord=nord, lower=lower, upper=upper

   if (NOT keyword_set(indir)) then indir = './'
   if (NOT keyword_set(outdir)) then outdir = './'
   if (NOT keyword_set(tmpdir)) then tmpdir = outdir

   if (N_elements(bkspace) EQ 0) then bkspace = 15
   if (N_elements(nord) EQ 0) then nord = 4
   if (N_elements(lower) EQ 0) then lower = 10
   if (N_elements(upper) EQ 0) then upper = 10

   nflat = N_elements(flatname)
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

   for iflat=0, nflat-1 do begin

      ;----------------------
      ; Read flat-field image

      sdssproc, flatname[iflat], flatimg, flativar, indir=indir, hdr=flathdr

      dims = size(flatimg, /dimens)
      nx = dims[0]
      ny = dims[1]

      ;----------------------
      ; Create the array of preliminary flats

      if (iflat EQ 0) then pixflatarr = fltarr(nx,ny,nflat)

      ;----------------------
      ; Determine YMODEL image
print, 'Working on file ', fullname[0]

      ; There are still systematics in the result from EXTRACT_IMAGE,
      ; so instead fit down each column of the image (which is slowly
      ; varying) by just doing a median filter or spline fit.
      ymodel = 0.0 * flatimg
      yaxis = findgen(ny)
      for i=0, nx-1 do begin
         print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])
;         ymodel[i,*] = median(transpose(flatimg[i,*]), 25)
         fullbkpt = slatec_splinefit(yaxis, flatimg[i,*], coeff, $
          invvar=flativar[i,*], bkspace=bkspace, nord=nord, $
          lower=lower, upper=upper, maxiter=3)
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
