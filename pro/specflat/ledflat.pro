;+
; NAME:
;   ledflat
;
; PURPOSE:
;   Compute pixel-flats from LED flat-field images
;
; CALLING SEQUENCE:
;   ledflat, expnum, [ docams=, indir=, threshrms=, threshlo=, nuse= ]
;
; INPUTS:
;   expnum     - Exposure number(s)
;
; OPTIONAL KEYWORDS:
;   docams     - Cameras to analyze; default to ['b1','b2','r1','r2'].
;   indir      - Input directory for files; default to searching for
;                files in $BOSS_SPECTRO_DATA/*.  If $BOSS_SPECTRO_DATA is not set,
;                then it is assumed to be /data/spectro.
;   threshrms  - Threshold for masking pixels whose fractional RMS (after
;                outlier-rejection) is larger than this value; default 0.05
;   threshlo   - Threshold for low pixels compared to local median;
;                default to 0.05
;   nuse       - Maximum number of images to use in computing statistics
;                at each pixel; default to 5 or the number of exposures,
;                whichever is smaller
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Output file names are pixflat-$DOCAMS-$MJD.fits, where the MJD
;   is from the first input file header.
;
;   This algorithm is designed to work best if there are several LED
;   flats observed with the LED moved around such that it is not always
;   illuminating the same way.  For any given pixel, the flat value
;   is determined from the image that has the lowest variance near that pix.
;
; EXAMPLES:
;   Solve for pixel flats from LED images using data from both
;   MJD 55141 and MJD 55142.  Select only those images with mean >5000 e-.
;   > ledflat, [103045+lindgen(5),103075+lindgen(3)], docams=['b1','b2']
;   > ledflat, [103070+lindgen(8)], docams=['r1','r2']
;
; BUGS:
;   Images should be bias-subtracted first???
;   Non-linearities should be applied first ???
;   Requires commenting-out the masking of the top+bottom of images in SDSSPROC!
;   Low-flux images should get automatically discarded.
;
;   This implementation has Poisson noise at the level of a single image.
;   A better implementation would perhaps first combine the exposures
;   from the same MJD or MJD+exptime first to beat down the Poisson noise.
;
; REVISION HISTORY:
;   18-Mar-2010  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro ledflat, expnum, docams=docams, indir=indir1, threshrms=threshrms1, $
 threshlo=threshlo1, nuse=nuse1

   defaultcams = ['b1','b2','r1','r2']
   if (NOT keyword_set(docams)) then docams = defaultcams

   ;----------
   ; Set defaults

   if (keyword_set(indir1)) then begin
      indir = indir1
   endif else begin
      indir = getenv('BOSS_SPECTRO_DATA')
      if (NOT keyword_set(indir)) then $
       indir = '/data/spectro'
      indir = indir + '/*'
   endelse

   if (keyword_set(threshrms1)) then threshrms = threshrms1 $
    else threshrms = 0.05
   if (keyword_set(threshlo1)) then threshlo = threshlo1 $
    else threshlo = 0.05
   if (keyword_set(nuse1)) then nuse = nuse1 $
    else nuse = (5 < n_elements(expnum))
   minval = 1 ; clip images at this min value

   ;----------
   ; If DOCAMS is an array, then call this routine recursively

   ncam = n_elements(docams)
   if (ncam GT 1) then begin
      for icam=0, ncam-1 do begin
         ledflat, expnum, docams=docams[icam], indir=indir, $
          threshrms=threshrms, threshlo=threshlo, nuse=nuse
      endfor
      return
   endif

   ;----------

   nfile = n_elements(expnum)
   filtsz = 25
   kern = fltarr(filtsz,filtsz) + 1

   for i=0, nfile-1 do begin
      ; Locate the input file (either compressed or un-compressed)
      filename1 = 'sdR-' + docams + '-' + string(expnum[i], format='(i8.8)') $
       + '.fit*'
      filename1 = (findfile(djs_filepath(filename1, root_dir=indir), $
       count=ct1))[0]
      if (ct1 EQ 0) then $
       message, 'File not found for EXPNUM=' + string(expnum[i])

      splog,'Reading file ',i+1,' of ', nfile,': ',filename1
      sdssproc, filename1, rawimg, rawivar, hdr=hdr, /applybias
      medimg = median(rawimg, filtsz) > minval
      filtimg = rawimg / medimg
      ; Compute the local variance around each pixel, but clip
      ; crazy values to cliprange
      cliprange = [0.90,1.10]
      varimg = convol( (((filtimg>cliprange[0])<cliprange[1])-1)^2, $
       kern, /center, /edge_truncate)

      ; Determine in which pixels this flat-field has lower local variance
      if (i EQ 0) then begin
         dims = size(rawimg,/dimens)
         imgarr = fltarr([dims,nfile])
         vararr = fltarr([dims,nfile])
         maskarr = bytarr([dims,nfile])
      endif
      imgarr[*,*,i] = filtimg
      vararr[*,*,i] = varimg
      maskarr[*,*,i] = rawivar EQ 0 OR rawimg LT minval ; set =0 for good pix

      if (i EQ 0) then thismjd = sxpar(hdr,'MJD')
   endfor
   rawimg = 0 ; clear memory
   rawivar = 0 ; clear memory
   medimg = 0 ; clear memory
   filtimg = 0 ; clear memory
   varimg = 0 ; clear memory

   ; Look for bad pixels
   ; Compute the stdev in each pix, after some outlier-rejection...
   if (nuse LE 2) then sigrej = 1.0 $ ; Irrelevant for only 1 or 2 flats
    else if (nuse EQ 3) then sigrej = 1.1 $
    else if (nuse EQ 4) then sigrej = 1.3 $
    else if (nuse EQ 5) then sigrej = 1.6 $
    else if (nuse EQ 6) then sigrej = 1.9 $
    else sigrej = 2.0

   ;----------
   ; Loop through each pixel, computing the mean and RMS from the nuse
   ; images with the lowest local variance (and preferring unmasked values)

   meanperpix = fltarr(dims)
   devperpix = fltarr(dims)
   totmask = lonarr(dims)
   for ix=0L, dims[0]-1L do begin
print,'Row ',ix,dims[0],string(13b),format='(a,i,i,a,$)'
   maxval = max(vararr) + 1
   for iy=0L, dims[1]-1L do begin
      isort = sort(vararr[ix,iy,*] + maskarr[ix,iy,*]*maxval)
      djs_iterstat, imgarr[ix,iy,isort[0:nuse-1]], sigrej=sigrej, $
       mean=mn1, sigma=sig1, mask=mask1
      meanperpix[ix,iy] = mn1
      devperpix[ix,iy] = sig1
      totmask[ix,iy] = total(mask1) ; number of good points
   endfor
   endfor


;   res = djs_avsigclip(imgarr, 3, sigrej=sigrej, maxiter=2, $
;    inmask=maskarr, outmask=outmask)
;   totmask = total(1-outmask,3)
;   ; Subtract the mean from each pix in imgarr...
;   meanperpix = total(imgarr*(1-outmask),3) / (totmask>1)
;   for i=0, nfile-1 do imgarr[*,*,i] -= meanperpix
;   devperpix = sqrt(total((imgarr*(1-outmask))^2,3) / (totmask>1))
;
;   ; Test the RMS per pix as a fractional RMS of the mean of that pixel
;   masklo = (stackimg LT threshlo) OR (totmask EQ 0) ; =1 bad
;   maskrms = (devperpix/(stackimg*(1-masklo)+masklo) GT threshrms) $
;    OR masklo ; =1 bad

   ; Test the RMS per pix as a fractional RMS of the mean of that pixel
   masklo = (meanperpix LT threshlo) OR (totmask EQ 0) ; =1 bad
   maskrms = (devperpix/(meanperpix*(1-masklo)+masklo) GT threshrms) $
    OR masklo ; =1 bad

   ; Mask those pixels at the edge of the field as bad in apermask,
   ; and unmask those same pixels in maskrms
   case docams of
      'b1': maxr = 2800
      'b2': maxr = 2600
      'r1': maxr = 2600
      'r2': maxr = 2600
   endcase
   rr = sqrt( (djs_laxisgen(dims,iaxis=0) - dims[0]/2. - 0.5)^2 $
    + (djs_laxisgen(dims,iaxis=1) - dims[1]/2. - 0.5)^2 )
   apermask = rr GT maxr ; =1 bad
   rr = 0 ; clear memory
   maskrms *= (1-apermask)

   ; Search for bad columns, and mask them to the end of the amplifier read.
   ; A bad column is assumed to begin whenever we have at least 40 out of 51
   ; pixels in a single column that are bad.
   ; Do this separately on the bottom and top half of the CCD, since the
   ; pixels are read out in different directions.
   ; But don't search for bad columns at a distance more than maxr from center
   badcolmask = bytarr(dims)
   kern1d = fltarr(1,51) + 1
   for ix=0, dims[0]-1 do begin
      ; Bottom half of CCD...
      bad1d = convol(maskrms[ix,0:dims[1]/2-1], kern1d, $
       /center, /edge_truncate) GE 40
      indx = where(bad1d, ct)
      if (ct GT 0) then begin
         badcolmask[ix,indx[0]:dims[1]/2-1] = 1B
         splog, 'Bad column on bottom ', ix
      endif
      ; Top half of CCD...
      bad1d = convol(maskrms[ix,dims[1]/2:dims[1]-1], kern1d, $
      /center, /edge_truncate) GE 40
      indx = reverse(where(bad1d, ct))
      if (ct GT 0) then begin
         badcolmask[ix,dims[1]/2:dims[1]/2+indx[0]] = 1B
         splog, 'Bad column on top ', ix
      endif
   endfor

   finalimg = meanperpix * (1-apermask) * (1-maskrms) * (1-badcolmask)

   ; Report fraction of bad pixels
   splog, 'Fraction of pix beyond optics aperture =', $
    mean(apermask)
   splog, 'Fraction of pix lost to bad columns =', $
    total(badcolmask*(1-apermask)) / total(1-apermask)
   splog, 'Fraction of pix lost to bad RMS =', $
    total(maskrms(1-badcolmask)*(1-apermask)) $
    / total((1-badcolmask)*(1-apermask))
   rms_report = [0.01,0.02,0.03,0.04,0.05,0.50]
   for i=0, n_elements(rms_report)-1 do $
    print, 'Fraction of pix deviant by ', rms_report[i], ' = ', $
     total((abs(1-finalimg) GT rms_report[i]) $
      * (1-apermask)*(1-maskrms)*(1-badcolmask)) $
      / total((1-apermask)*(1-apermask)*(1-maskrms)*(1-badcolmask))

   ; Use the MJD from the first file in the output name
   mjdstr = string(thismjd,format='(i5.5)')
   outfile = 'pixflat-'+docams+'-'+mjdstr+'.fits'
   splog, 'Writing '+outfile
   mwrfits, finalimg, outfile, hdr, /create

   return
end
;------------------------------------------------------------------------------
