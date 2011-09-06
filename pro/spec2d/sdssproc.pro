;+
; NAME:
;   sdssproc
;
; PURPOSE:
;   Read in raw SDSS files, and process with opConfig, opECalib, opBC par files.
;
; CALLING SEQUENCE:
;   sdssproc, infile, [image, invvar, indir=, $
;    outfile=, nsatrow=, fbadpix=, $
;    hdr=hdr, configfile=, ecalibfile=, bcfile=, $
;    /applybias, /applypixflat, /silent, /do_lock, minflat=, maxflat=, $
;    spectrographid=, color=, camname=, /applycrosstalk, ccdmask= ]
;
; INPUTS:
;   infile     - Raw SDSS file name
;
; OPTIONAL KEYWORDS:
;   indir      - Input directory for INFILE
;   outfile    - Calibrated 2D frame, HDU#0 is the image and #1 is invvar;
;                if set as /OUTFILE, then default name is constructed from
;                INFILE as sdProc-XX-XXXXXXXX.fits, excluding any path info
;   nsatrow    - Number of saturated rows, assuming that a row is saturated
;                if at least 20 of its pixels are above saturation level
;   fbadpix    - Fraction of bad pixels, not including bad columns
;   configfile - Default to "opConfig*par", selecting the file with the
;                appropriate MJD.
;   ecalibfile - Default to "opECalib*par", selecting the file with the
;                appropriate MJD.
;   bcfile     - Default to "opBC*par", selecting the file with the
;                appropriate MJD.
;   applybias  - Apply 2-D bias image.
;   applypixflat- Apply 2-D pixel-to-pixel flat (after subtracting bias).
;   silent     - If set, then don't output any text.
;   do_lock    - If set, then lock the "sdHdrFix-$MJD.par" file
;                using DJS_LOCKFILE().
;   minflat    - Minimum values allowed for pixflat; pixels with the
;                flat out of range are masked; default to 0.
;   maxflat    - Maximum values allowed for pixflat; pixels with the
;                flat out of range are masked; default to 1e10.
;   applycrosstalk - choose whether or not to apply crosstalk correction
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   image      - Processed 2d image
;   invvar     - Associated inverse variance
;   hdr        - Processed FITS header
;   spectrographid - Return spectrograph ID (1 or 2)
;   color      - Return spectrograph color ('red' or 'blue')
;   camname    - Return camera name: 'b1', 'r1', 'b2', or 'r2'
;   ccdmask    - Image with 'NODATA' bit set for parts of BOSS data
;
; COMMENTS:
;   Only the header is read from the image if IMAGE, INVVAR, OUTFILE and
;   VARFILE are all not set.
;
;   Required header keywords: EXPTIME.
;
;   The returned image is in electrons, not ADU.
;
;   The signal-to-noise is limited to never exceed 100, by adding 1.e-4
;   times the flux squared to the variance term.
;
;   Change the CAMERAS keyword to the camera as specified by the file name.
;
;   Rename 'target' to 'science', and 'calibration' to 'arc' in the
;   header keyword FLAVOR.
;
;   Determine the exposure number from the file name itself.
;
; BUGS:
;   The open-shutter correction SMEARIMG will include smeared data from
;   any cosmic rays, which is wrong.  At the minimum, I could interpolate
;   over A/D saturations (in ADMASK) before constructing SMEARIMG.
;
; PROCEDURES CALLED:
;   djs_filepath()
;   djs_iterstat
;   fileandpath()
;   findopfile()
;   fits_purge_nans
;   fits_wait()
;   headfits()
;   idlspec2d_version()
;   idlutils_version()
;   lookforgzip()
;   mwrfits
;   rdss_fits()
;   sphdrfix
;   splog
;   sxaddpar
;   sxpar()
;   yanny_free
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;   make_badcolumn_mask()
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/opConfig*par
;   $IDLSPEC2D_DIR/examples/opECalib*par
;   $IDLSPEC2D_DIR/examples/opBC*par
;   $SPECFLAT_DIR/biases/pixbias*.fits*
;   $SPECFLAT_DIR/flats/pixflat*.fits*
;
; REVISION HISTORY:
;   13-May-1999  Written by Scott Burles & David Schlegel, Apache Point.
;   08-Sep-1999  Modified to read Yanny param files instead of FITS
;                versions of the same (DJS).
;   01-Dec-1999  Added version stamping (DJS).
;   07-Dec-1999  Mask neighbors of pixels that saturated the A/D converter.
;                Identify blead trails and mask from that row up (DJS).
;   10-Dec-1999  Test if the shutter was open during readout, and try
;                to correct the light for that (DJS).
;   04-Feb-2000  Declare that the shutter was open if it is a >640 sec
;                exposure taken before MJD=51570 (DJS).
;   26-Jul-2000  Added fix for "dropped pixel" problem for data on or after
;                MJD=51688 (23 May 2000).  Should disable this code for later
;                MJD's once this problem is fixed in the electronics.
;   26-Jul-2000  Added fix for more severe "shifted row" electronics problem
;                for data taken on MJD=51578 to 51580 (3 to 5 Feb 2000).
;   26-Aug-2000  Horrible kludge for this night MJD=51779 (22/23 Aug 2000).
;                Add a noise term of 100 ADU to the left amplifier of r2.
;                That amplifier had random bits being set incorrectly,
;                in particular the 32 bit and 256 bit.
;   04-Nov-2000  Measure the bias values using DJS_ITERSTAT instead of MEDIAN,
;                since the median is always only an integer value.
;   31-Jan-2001  Determine the exposure number from the file name itself,
;                since the counting got off by one on MJD=51882.
;   08-Aug-2011: Changing bias-subtraction recipe for BOSS (A. Bolton, Utah)
;-
;------------------------------------------------------------------------------
pro sdssproc_badformat, image, camname=camname, mjd=mjd

   dims = size(image,/dimens)
   nx = dims[0]
   ny = dims[1]

   case camname of
   'b1': begin
       xs = [77,4274,77,4274]
       xw = [1,1,1,1]
       thresh = [15,15,0,0]
       end
   'b2': begin
       xs = [77,4274,77,4274]
       xw = [1,1,1,1]
       thresh = [5,10,10,0]
       end
   'r1': begin
       ; The transients on the right of r1 are really 2 columns wide
       xs = [111,4240,111,4240]
       xw = [1,2,1,2]
       thresh = [30,20,60,15]
       end
   'r2': begin
       if (mjd LT 55300) then begin
          xs = [111,4240,112,4239]
          xw = [1,1,1,1]
          thresh = [20,10,5,5]
       endif else begin
          xs = [112,4239,111,4240]
          xw = [1,1,1,1]
          thresh = [5,10,12,5]
       endelse
       end
   endcase

   nsub = 7
   reg0 = image[xs[0]:xs[0]+nsub-1,0:ny/2-1]
   reg1 = rotate(image[xs[1]-nsub+1:xs[1],0:ny/2-1],5)
   reg2 = rotate(image[xs[2]:xs[2]+nsub-1,ny/2:ny-1],7)
   reg3 = rotate(image[xs[3]-nsub+1:xs[3],ny/2:ny-1],2)

   if (xw[0] GT 1) then $
    reg0 = convol(reg0, [intarr(xw[0]-1),1,intarr(xw[0]-1)+1], /edge_wrap)
   if (xw[1] GT 1) then $
    reg1 = convol(reg1, [intarr(xw[1]-1),1,intarr(xw[1]-1)+1], /edge_wrap)
   if (xw[2] GT 1) then $
    reg2 = convol(reg2, [intarr(xw[2]-1),1,intarr(xw[2]-1)+1], /edge_wrap)
   if (xw[3] GT 1) then $
    reg3 = convol(reg3, [intarr(xw[3]-1),1,intarr(xw[3]-1)+1], /edge_wrap)

   shiftvec = lonarr(ny/2) ; Pixel shift of each row (same in all quadrants)
   nshift = lonarr(nsub) ; Keep track of # of rows shifted each pixel amount

   peak = intarr(4)
   for iy=0L, ny/2-1L do begin
      ; Find the peak value relative to where it should be in this row
      max0 = max(reg0[0:nsub-1-xw[0],iy], peak0)
      max1 = max(reg1[0:nsub-1-xw[1],iy], peak1)
      max2 = max(reg2[0:nsub-1-xw[2],iy], peak2)
      max3 = max(reg3[0:nsub-1-xw[3],iy], peak3)
      ; Test that the peak is at least "thresh" larger than the next pix,
      ; where that threshold is different in each amplifier
      peak0 *= (reg0[peak0,iy] GE reg0[peak0+xw[0],iy] + thresh[0])
      peak1 *= (reg1[peak1,iy] GE reg1[peak1+xw[1],iy] + thresh[1])
      peak2 *= (reg2[peak2,iy] GE reg2[peak2+xw[2],iy] + thresh[2])
      peak3 *= (reg3[peak3,iy] GE reg3[peak3+xw[3],iy] + thresh[3])
      peak = [peak0,peak1,peak2,peak3]
      bpeak = peak[(sort(peak))[1]]
      ; If at least 3 quadrants appear shifted the same amount, then
      ; use that shift
      if (total(peak EQ bpeak) GE 3) then begin
         shiftvec[iy] = bpeak
         nshift[bpeak]++
      endif else begin
         ; If no shifts detected, then assume it is the same as the previous row
         if (iy GT 0) then shiftvec[iy] = shiftvec[iy-1]
      endelse
   endfor

   ; Do the actual image shifts
   for iy=0L, ny/2-1L do begin
      if (shiftvec[iy] GT 0) then begin
         image[0:nx/2-1,iy] = shift(image[0:nx/2-1,iy], -shiftvec[iy])
         image[0:nx/2-1,ny-1-iy] = shift(image[0:nx/2-1,ny-1-iy], -shiftvec[iy])
         image[nx/2:nx-1,iy] = shift(image[nx/2:nx-1,iy], shiftvec[iy])
         image[nx/2:nx-1,ny-1-iy] = shift(image[nx/2:nx-1,ny-1-iy], shiftvec[iy])
      endif
   endfor

   nunknown = ny - long(2*total(nshift))
   if (nunknown GT 0) then $
    splog, 'Warning: Serial transient not detected in ', nunknown, ' rows'

   for i=1, nsub-1 do $
    if (nshift[i] GT 0) then $
     splog, 'WARNING: Electronics shifted ', 2*nshift[i], ' rows by ', i, ' pix'

   return
end
;------------------------------------------------------------------------------
;  Create the bad column mask (1 for a masked pixel) with image size nc,nr
;  If the operation doesn't work just return 0 for no masked pixels

function make_badcolumn_mask, bcfile, camrow, camcol, nc=nc, nr=nr, $
    silent=silent
    
  if NOT keyword_set(nc) then nc=2048L
  if NOT keyword_set(nr) then nr=2048L
  
  yanny_read, bcfile, pdata
  if (size(pdata,/tname)) EQ 'INT' then begin
    if (NOT keyword_set(silent)) then $
      splog, 'WARNING: Could not read BC file ' + fileandpath(bcfile)
    return, 0
  endif
  
  bc = *pdata[0]
  yanny_free, pdata
  
  ibc = where(bc.camrow EQ camrow AND bc.camcol EQ camcol, nbc)
  if (NOT keyword_set(nbc)) then begin
    if (NOT keyword_set(silent)) then $
      splog,'WARNING: Could not find this camera info in BC file ' $
      + fileandpath(bcfile)
    return, 0
  endif
  
  bc = bc[ibc]
  bcmask = bytarr(nc, nr)
  
  ; Mask out bad columns
  
  if (nbc GT 0) then begin
    bcsc = (bc.dfcol0 > 0) < (nc-1)
    bcec = (bc.dfcol0 + bc.dfncol - 1 < (nc-1)) > bcsc
    bcsr = (bc.dfrow0 > 0) < (nr-1)
    bcer = (bc.dfrow0 + bc.dfnrow - 1 < (nr-1)) > bcsr
    
    for i=0, nbc-1 do bcmask[bcsc[i]:bcec[i],bcsr[i]:bcer[i]] = 1
  endif
  
  return, bcmask
end

;------------------------------------------------------------------------------
pro sdssproc, infile1, image, invvar, indir=indir, $
 outfile=outfile1, nsatrow=nsatrow, fbadpix=fbadpix, $
 hdr=hdr, configfile=configfile, ecalibfile=ecalibfile, bcfile=bcfile, $
 applybias=applybias, applypixflat=applypixflat, silent=silent, $
 do_lock=do_lock, minflat=minflat, maxflat=maxflat, $
 spectrographid=spectrographid, color=color, camname=camname, $
 applycrosstalk=applycrosstalk, ccdmask=ccdmask

  common com_sdssproc, vers2d, versutils, versflat, verslog
    
  if (N_params() LT 1) then begin
    doc_library, 'sdssproc'
    return
  endif

  infile = infile1[0]
  readimg = arg_present(image) OR keyword_set(outfile1)
  readivar = arg_present(invvar) OR keyword_set(outfile1) $
    OR arg_present(nsatrow) OR arg_present(fbadpix)

  fullname = djs_filepath(infile, root_dir=indir)
  fullname = (lookforgzip(fullname, count=ct))[0]
  ;   fullname = (findfile(fullname, count=ct))[0]
  if (ct NE 1) then $
    message, 'Cannot find image ' + infile
    
  if (keyword_set(outfile1)) then begin
     if (size(outfile1,/tname) EQ 'STRING') then begin
        outfile = outfile1
     endif else begin
        outfile = 'sdProc-' + strmid(fileandpath(infile),4,11)+'.fits'
     endelse
  endif

  if (readimg OR readivar) then $
    rawdata = rdss_fits(fullname, hdr, /nofloat, silent=silent) $
  else $
    hdr = headfits(fullname)

  sxdelpar, hdr, 'BZERO'
  sxdelpar, hdr, 'BSCALE'
  sxdelpar, hdr, 'CHECKSUM'
  sxdelpar, hdr, 'DATASUM'

  ;-----------
  ; Determine the exposure number from the file name itself.
  ; Very bad form, but this information is sometimes wrong in the header.
  ; In particular, the counting got off by 1 on MJD=51882.
  
  i = strpos(infile, '-', /reverse_search)
  expnum = long( strmid(infile, i+1, 8) )
  if (NOT keyword_set(expnum)) then $
   message, 'Cannot determine exposure number from file name ' + infile

  ;-----------
  ; Fix the headers with any hand-edits that we have determined.
    
  ; Fix for BOSS test data
  if (expnum GE 100053 AND expnum LE 100114) then sxaddpar, hdr, 'MJD', 55050
  if (expnum GE 100115 AND expnum LE 100160) then sxaddpar, hdr, 'MJD', 55051
  if (expnum GE 100162 AND expnum LE 100186) then sxaddpar, hdr, 'MJD', 55052
  if (expnum GE 100272 AND expnum LE 100287) then sxaddpar, hdr, 'MJD', 55061
  if (expnum GE 100288 AND expnum LE 100298) then sxaddpar, hdr, 'MJD', 55062
  if (expnum GE 100299 AND expnum LE 100304) then sxaddpar, hdr, 'MJD', 55063
  if (expnum GE 100371 AND expnum LE 100397) then sxaddpar, hdr, 'MJD', 55065
  if (expnum GE 100398 AND expnum LE 100460) then sxaddpar, hdr, 'MJD', 55066
  if (expnum GE 100461 AND expnum LE 100468) then sxaddpar, hdr, 'MJD', 55067
  if (expnum GE 100469 AND expnum LE 100532) then sxaddpar, hdr, 'MJD', 55068
  if (expnum GE 100533 AND expnum LE 100567) then sxaddpar, hdr, 'MJD', 55069
  if (expnum GE 100568 AND expnum LE 100653) then sxaddpar, hdr, 'MJD', 55070
  if (expnum GE 100655 AND expnum LE 100705) then sxaddpar, hdr, 'MJD', 55071
  if (expnum GE 100706 AND expnum LE 100745) then sxaddpar, hdr, 'MJD', 55072
  if (expnum GE 100746 AND expnum LE 100781) then sxaddpar, hdr, 'MJD', 55073

  sphdrfix, infile, hdr, silent=silent, do_lock=do_lock
  
  ;-----------
  ; Replace exposure number with that found in the file name.
  
  hdrexp = sxpar(hdr, 'EXPOSURE')
  if (expnum NE hdrexp) then begin
    if (NOT keyword_set(silent)) then $
      splog, 'WARNING: Exposure number in header (' $
      + strtrim(string(hdrexp),2) + ') disagrees w/filename (' $
      + strtrim(string(expnum),2) + ') !!'
    sxaddpar, hdr, 'EXPOSURE', expnum
  endif
  
  ;-----------
  ; Determine which CCD from the file name itself, using either the
  ; numbering scheme (01,02,03,04) or naming scheme (b1,r2,b2,r1).
  ; Very bad form, but this information is not in the header since
  ; the CAMERAS keyword is sometimes wrong.
  
  i = strpos(infile, '-', /reverse_search)
  if (i[0] EQ -1 OR i-2 LT 0) then $
    message, 'Cannot determine CCD number from file name ' + infile
    
  camnames = ['b1', 'r2', 'b2', 'r1']
  camnums = ['01', '02', '03', '04']
  
  ; First try to match a camera name (e.g., 'b1'), then try to match
  ; a camera number (e.g., '01').  If both fail, then abort.
  indx = where(strmid(infile, i-2, 2) EQ camnames, ct)
  if (ct NE 1) then $
    indx = where(strmid(infile, i-2, 2) EQ camnums, ct)
  if (ct NE 1) then $
    message, 'Cannot determine CCD number from file name ' + infile
    
  ; Do not read the camera from the CAMERAS keyword, since this was often
  ; wrong in the early days!
  ;   cameras = strtrim( sxpar(hdr, 'CAMERAS'), 2 )
  camname = camnames[indx[0]]
  case camname of
    'b1': begin
      spectrographid = 1
      color = 'blue'
    end
    'r1': begin
      spectrographid = 1
      color = 'red'
    end
    'b2': begin
      spectrographid = 2
      color = 'blue'
    end
    'r2': begin
      spectrographid = 2
      color = 'red'
    end
  endcase
  camcol = indx[0] + 1
  camrow = 0

  ; Cache these versions in a common block, since the calls to spawn are slow
  if (NOT keyword_set(vers2d)) then begin
     vers2d = idlspec2d_version()
     versutils = idlutils_version()
     spawn, 'speclog_version', verslog, err, /noshell
     if (NOT keyword_set(verslog)) then verslog = 'Unknown'
     spawn, 'specflat_version', versflat, err, /noshell
     if (NOT keyword_set(versflat)) then versflat = 'Unknown'
  endif
  
  sxaddpar, hdr, 'CAMROW', camrow
  sxaddpar, hdr, 'CAMCOL', camcol
  sxaddpar, hdr, 'TELESCOP', 'SDSS 2.5-M', ' Sloan Digital Sky Survey'
  sxaddpar, hdr, 'AUTHOR', 'Scott Burles & David Schlegel'
  sxaddpar, hdr, 'VERSIDL', !version.release, ' Version of IDL'
  sxaddpar, hdr, 'VERSUTIL', versutils, ' Version of idlutils'
  sxaddpar, hdr, 'VERSREAD', vers2d, $
    ' Version of idlspec2d for pre-processing raw data', after='VERSUTIL'
  sxaddpar, hdr, 'VERSLOG', verslog[0], ' Version of SPECLOG product', $
    after='VERSREAD'
  sxaddpar, hdr, 'VERSFLAT', versflat[0], ' Version of SPECFLAT product', $
    after='VERSFLAT'
    
  ;-----------
  ; Rename 'target' -> 'science', and 'calibration' -> 'arc'
    
  mjd = sxpar(hdr, 'MJD')
  flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)
  if (flavor EQ 'target') then flavor = 'science'
  if (flavor EQ 'calibration') then flavor = 'arc'
  if (mjd GT 51576) then begin
    if (sxpar(hdr, 'COLBIN') GT 1 OR $
      sxpar(hdr, 'ROWBIN') GT 1) then flavor = 'unknown'
  endif
  sxaddpar, hdr, 'FLAVOR', flavor
  sxaddpar, hdr, 'CAMERAS', camname

  ;--------------
  ; Flag to trigger bolton bias subtraction for survey-quality BOSS data:
  if (mjd ge 55170) then bossgood = 1B else bossgood = 0B

  ; Mark when the BOSS red CCDs switched from 1-phase to 2-phase readout
  if (mjd GE 55415 AND strmatch(camname,'r*')) then $
   sxaddpar, hdr, 'TWOPHASE', 'T' $
  else $
   sxaddpar, hdr, 'TWOPHASE', 'F'
  
  ;-----------
  ; Deal with the first light BOSS data on MJD 55052

if (mjd GE 55052) then begin
   ; Declare any BOSS exposures 'excellent' if that keyword missing,
   ; unless it was a Hartmann exposure
   junk = sxpar(hdr, 'QUALITY', count=ct)
   if (ct EQ 0) then begin
      hartmann = strtrim(sxpar(hdr, 'HARTMANN'),2)
      if (hartmann EQ 'Left' OR hartmann EQ 'Right') then $
       sxaddpar, hdr, 'QUALITY', 'test' $
      else $
       sxaddpar, hdr, 'QUALITY', 'excellent'
   endif
   if (readimg OR readivar) then begin
      rawdata += 32768.
      sdssproc_badformat, rawdata, camname=camname, mjd=mjd

      ; Rotate these images to put fiber #1 on the left,
      ; and have wavelength running up.
      if (mjd LT 55113) then rawdata = rotate(rawdata,2)

      case strmid(camname,0,1) of
      'b': begin
         case spectrographid of
            1: gain = [1.048, 1.048, 1.018, 1.006] ; b1 gain
            2: gain = [1.040, 0.994, 1.002, 1.010] ; b2 gain
         end
         ; Do bolton bias subtraction for survey-quality BOSS dates:
         ; (Note that these lines are identical between b and r cams.)
         if bossgood then begin
            if (keyword_set(silent) EQ 0) then $
               splog, 'BOSS survey-quality MJD: applying pixbias whether you like it or not!'
            pp = filepath('', root_dir=getenv('SPECFLAT_DIR'), subdirectory='biases')
            pixbiasname = findopfile('boss_pixbias-*-'+camname+'.fits*', mjd, pp, $
                                     silent=silent)
            image = bolton_biassub(rawdata, pp+pixbiasname, rnoise=rnoise, $
                                   cam=camname, sigthresh=3.0)
            xwid = (size(image))[1]
            ywid = (size(image))[2]
            rdnoise = rnoise[*] * 1.015 * gain ; account for sigma-clipping
            if (keyword_set(silent) EQ 0) then $
               splog, 'Read noise = ', rdnoise, ' electrons'
            image[0:xwid/2-1,0:ywid/2-1] *= gain[0]
            image[xwid/2:xwid-1,0:ywid/2-1] *= gain[1]
            image[0:xwid/2-1,ywid/2:ywid-1] *= gain[2]
            image[xwid/2:xwid-1,ywid/2:ywid-1] *= gain[3]
         endif else begin
            image = fltarr(4096,4112)
            bias = [ median(rawdata[10:67,56:2111]), $
                     median(rawdata[4284:4340,56:2111]), $
                     median(rawdata[10:67,2112:4167]), $
                     median(rawdata[4284:4340,2112:4167]) ]
         ; Subtract the bias before computing noise for numerical reasons
            rdnoise = [ djsig(rawdata[10:67,56:2111] - bias[0], sigrej=3.0), $
                        djsig(rawdata[4284:4340,56:2111] - bias[1], sigrej=3.0), $
                        djsig(rawdata[10:67,2112:4167] - bias[2], sigrej=3.0), $
                        djsig(rawdata[4284:4340,2112:4167] - bias[3], sigrej=3.0) ]
            rdnoise *= 1.015 * gain ; account for sigma-clipping
            if (keyword_set(silent) EQ 0) then $
               splog, 'Read noise = ', rdnoise, ' electrons'
            image[0:2047,0:2055] = gain[0] * (rawdata[128:2175,56:2111] - bias[0])
            image[2048:4095,0:2055] = gain[1] * (rawdata[2176:4223,56:2111] - bias[1])
            image[0:2047,2056:4111] = gain[2] * (rawdata[128:2175,2112:4167] - bias[2])
            image[2048:4095,2056:4111] = gain[3] * (rawdata[2176:4223,2112:4167] - bias[3])
         endelse
         if (readivar) then begin
            invvar = 0.*image
            invvar[0:2047,0:2055] = $
             1. / (rdnoise[0]^2 + (image[0:2047,0:2055]>0))
            invvar[0:2047,2056:4111] = $
             1. / (rdnoise[1]^2 + (image[0:2047,2056:4111]>0))
            invvar[2048:4095,0:2055] = $
             1. / (rdnoise[2]^2 + (image[2048:4095,0:2055]>0))
            invvar[2048:4095,2056:4111] = $
             1. / (rdnoise[3]^2 + (image[2048:4095,2056:4111]>0))
            mask = 0.*image
            mask[*,696:3516] = 1
            invvar *= mask
         endif
         end
      'r': begin
         case spectrographid of
            1: if (mjd LT 55131) then gain = [2.63, 2.45, 2.30, 2.72] $
             else if (mjd LT 55800) then gain = [1.966, 1.566, 1.542, 1.546] $
             else gain = [1.9253, 1.5122, 1.4738, 1.5053]  ; R1 replaced summer 2011
            2: if (mjd LT 55131) then gain = [1.89, 1.51, 1.40, 1.44] $
             else if (mjd LT 55141) then gain = [2.66, 2.53, 2.02, 3.00] $
             else if (mjd LT 55300) then gain = [1.956, 1.618, 1.538, 1.538] $
             else gain = [1.598, 1.656, 1.582, 1.594] ; R2 replaced April 2011
         end
          ; Do bolton bias subtraction for survey-quality BOSS dates:
         ; (Note that these lines are identical between b and r cams.)
         if bossgood then begin
            if (keyword_set(silent) EQ 0) then $
               splog, 'BOSS survey-quality MJD: applying pixbias whether you like it or not!'
            pp = filepath('', root_dir=getenv('SPECFLAT_DIR'), subdirectory='biases')
            pixbiasname = findopfile('boss_pixbias-*-'+camname+'.fits*', mjd, pp, $
                                     silent=silent)
            image = bolton_biassub(rawdata, pp+pixbiasname, rnoise=rnoise, $
                                   cam=camname, sigthresh=3.0)
            xwid = (size(image))[1]
            ywid = (size(image))[2]
            rdnoise = rnoise[*] * 1.015 * gain ; account for sigma-clipping
            if (keyword_set(silent) EQ 0) then $
               splog, 'Read noise = ', rdnoise, ' electrons'
            image[0:xwid/2-1,0:ywid/2-1] *= gain[0]
            image[xwid/2:xwid-1,0:ywid/2-1] *= gain[1]
            image[0:xwid/2-1,ywid/2:ywid-1] *= gain[2]
            image[xwid/2:xwid-1,ywid/2:ywid-1] *= gain[3]
         endif else begin
            image = fltarr(4114,4128)
            bias = [ median(rawdata[10:100,48:2111]), $
                     median(rawdata[4250:4340,48:2111]), $
                     median(rawdata[10:100,2112:4175]), $
                     median(rawdata[4250:4340,2112:4175]) ]
            ; Subtract the bias before computing noise for numerical reasons
            rdnoise = [ djsig(rawdata[10:100,48:2111] - bias[0], sigrej=3.0), $
                        djsig(rawdata[4250:4340,48:2111] - bias[1], sigrej=3.0), $
                        djsig(rawdata[10:100,2112:4175] - bias[2], sigrej=3.0), $
                        djsig(rawdata[4250:4340,2112:4175] - bias[3], sigrej=3.0) ]
            rdnoise *= 1.015 * gain ; account for sigma-clipping
            if (keyword_set(silent) EQ 0) then $
               splog, 'Read noise = ', rdnoise, ' electrons'
            if (mjd LT 55067) then begin
               image[0:2055,0:2063] = gain[0] * (rawdata[120:2175,48:2111] - bias[0])
               image[2056,0:2063] = gain[0] * (rawdata[0,48:2111] - bias[0])
               image[2058:4113,0:2063] = gain[1] * (rawdata[2176:4231,48:2111] - bias[1])
               image[2057,0:2063] = gain[1] * (rawdata[4351,48:2111] - bias[1])
               image[0:2055,2064:4127] = gain[2] * (rawdata[120:2175,2112:4175] - bias[2])
               image[2056,2064:4127] = gain[2] * (rawdata[0,2112:4175] - bias[2])
               image[2058:4113,2064:4127] = gain[3] * (rawdata[2176:4231,2112:4175] - bias[3])
               image[2057,2064:4127] = gain[3] * (rawdata[4351,2112:4175] - bias[3])
            endif else begin
               ; Craig fixed the r1,r2 image formatting on MJD 55067
               image[0:2056,0:2063] = gain[0] * (rawdata[119:2175,48:2111] - bias[0])
               image[2057:4113,0:2063] = gain[1] * (rawdata[2176:4232,48:2111] - bias[1])
               image[0:2056,2064:4127] = gain[2] * (rawdata[119:2175,2112:4175] - bias[2])
               image[2057:4113,2064:4127] = gain[3] * (rawdata[2176:4232,2112:4175] - bias[3])
            endelse
         endelse
         if (readivar) then begin
            invvar = 0.*image
            invvar[0:2056,0:2063] = $
             1. / (rdnoise[0]^2 + (image[0:2056,0:2063]>0))
            invvar[0:2056,2064:4127] = $
             1. / (rdnoise[1]^2 + (image[0:2056,2064:4127]>0))
            invvar[2057:4113,0:2063] = $
             1. / (rdnoise[2]^2 + (image[2057:4113,0:2063]>0))
            invvar[2057:4113,2064:4127] = $
             1. / (rdnoise[3]^2 + (image[2057:4113,2064:4127]>0))
            mask = 0.*image
            mask[*,28:3668] = 1
            invvar *= mask
         endif
         end
      endcase

      if (arg_present(ccdmask)) then begin
         ccdmask = lonarr(size(image,/dimens))
         ccdmask += (mask EQ 0) * sdss_flagval('SPPIXMASK','NODATA')
      endif

      ; Add read-noise to header
      sxaddpar, hdr, 'RDNOISE0', rdnoise[0], 'CCD read noise amp 0 [electrons]'
      sxaddpar, hdr, 'RDNOISE1', rdnoise[1], 'CCD read noise amp 1 [electrons]'
      sxaddpar, hdr, 'RDNOISE2', rdnoise[2], 'CCD read noise amp 2 [electrons]'
      sxaddpar, hdr, 'RDNOISE3', rdnoise[3], 'CCD read noise amp 3 [electrons]'

      ; Identify CRs, and grow by 1 pix
      if (keyword_set(invvar)) then begin
;         psfvals = [0.625,0.391] ; for FWHM=2.5 pix
         psfvals = [0.496,0.246] ; for FWHM=2.0 pix
         reject_cr, image, invvar, psfvals, rejects, c2fudge=0.8, niter=6, $
          nrejects=nrejects
         if (nrejects GT 0) then begin
            crmask = 0. * image
            crmask[rejects] = 1
            crmask = dilate(crmask, replicate(1,3,3))
            invvar *= (crmask EQ 0)
         endif
      endif
      if (arg_present(fbadpix)) then fbadpix = mean(invvar EQ 0 AND mask EQ 1)
   end

endif else begin ; SDSS-I data

  ;-----------
  ; If MJD <= 51813, then set QUALITY=excellent unless this is over-written
  ; by SPHDRFIX.  This is because for the early data, this keyword was
  ; arbitrarily set to 'unknown', 'bad', 'acceptable', or 'excellent'.
  
  if (mjd LE 51813) then sxaddpar, hdr, 'QUALITY', 'excellent'
  
  ;-----------
  ; If the OBSCOMM keyword contains the words "dithered" or "focus",
  ; then assume this is test data (set QUALITY=test).
  
  obscomm = sxpar(hdr,'OBSCOMM')
  if (strmatch(obscomm,'*dithered*') OR strmatch(obscomm,'*focus*')) then $
    sxaddpar, hdr, 'QUALITY', 'test'
    
  ;-----------
  ; Check that this was not a non-test exposure taken during the daytime.

  if (keyword_set(silent) EQ 0 AND strmatch(flavor,'science*')) then $ 
   warn_daytime, hdr

  ;-----------
  ; Dispose of other invalid header keywords
  
  fits_purge_nans, hdr, verbose=(keyword_set(silent) EQ 0)
  
  ;-----------
  ; Return if only the header (and no data) was requested.
  
  if (NOT readimg AND NOT readivar) then return
  
  ;-----------
  ; Fix the shifted-row problem in the electronics that appeared
  ; on MJD=51578 to 51580 (3-5 Feb 2000) for spectrograph-2.
  ; Note that the bad rows need to be identified on the red frame,
  ; so if we are reducing b2, we need to read in the image for r2
  ; to identify the bad rows.
  
  if ((mjd GE 51578 AND mjd LE 51580) $
    AND (camname EQ 'b2' OR camname EQ 'r2') $
    AND (readimg OR readivar)) then begin
    if (camname EQ 'b2') then begin
      ;         i1 = strpos(infile,'b2',/reverse_search) ; IDL 5.3 command
      i1 = strpos(infile,'b2',/reverse_search)
      if (i1 EQ -1) then $
        message, 'Unable to parse corresponding red file for '+infile
      redfile = infile
      strput, redfile, 'r2', i1
      redfile = (lookforgzip(djs_filepath(redfile, root_dir=indir)))[0]
      
      if (fits_wait(redfile, deltat=1, tmax=1)) then $
        reddata = rdss_fits(redfile, /nofloat, silent=silent) $
      else $
        if (NOT keyword_set(silent)) then $
        splog, 'Warning: Could not read corresponding red file ' + redfile
    endif else begin
      reddata = rawdata
    endelse
    
    if (keyword_set(reddata)) then $
      ibad = where( reddata[20,*] LT median(reddata[20,*]) - 100 , nbad ) $
    else $
      nbad = 0
      
    if (nbad GT 0) then begin
      if (NOT keyword_set(silent)) then $
        splog, 'WARNING: Fixing ', nbad, ' shifted rows (from electronics)'
        
      ; For unrecoverable data, set KILLDATA=0
      killdata = byte(0*rawdata) + 1b
      
      medval = median(rawdata[22:39,*])
      for ii=0, nbad-1 do begin
        nshift = $
          (reverse(where(rawdata[20:39,ibad[ii]] GT medval + 100)))[0]
          
        ; If nshift is -1, we will also zero the whole row
          
        if (nshift GE 0 AND nshift LT 10 ) then begin
          rawdata[0:1063-nshift,ibad[ii]] = rawdata[nshift:1063,ibad[ii]]
          killdata[1063-nshift+1:1063,ibad[ii]] = 0
        endif else begin
          killdata[*,ibad[ii]] = 0
        endelse
        
      endfor
    endif
    
    reddata = 0  ; Free memory
  endif
  
  ;-----------
  ; Fix the "dropped pixel" problem in the electronics that appeared
  ; on MJD=51688 (23 May 2000) for spectrograph-2.
  ; Note that the bad rows need to be identified on the red frame,
  ; so if we are reducing b2, we need to read in the image for r2
  ; to identify the bad rows.
  ; Reference e-mail discussion with JEG on 01-Jun-2000.
  
  if ((mjd GE 51688 AND mjd LT 55025) $
    AND (camname EQ 'b2' OR camname EQ 'r2') $
    AND (readimg OR readivar)) then begin
    if (camname EQ 'b2') then begin
      i1 = strpos(infile,'b2',/reverse_search) ; IDL 5.3 command
      if (i1 EQ -1) then $
        message, 'Unable to parse corresponding red file for '+infile
      redfile = infile
      strput, redfile, 'r2', i1
      redfile = (lookforgzip(djs_filepath(redfile, root_dir=indir)))[0]
      
      if (fits_wait(redfile, deltat=1, tmax=1)) then $
        reddata = rdss_fits(redfile, /nofloat, silent=silent) $
      else $
        if (NOT keyword_set(silent)) then $
        splog, 'Warning: Could not read corresponding red file ' + redfile
    endif else begin
      reddata = rawdata
    endelse
    
    if (keyword_set(reddata)) then $
      ibad = where( reddata[20,*] LT median(reddata[20,*]) - 100 $
      AND reddata[2107,*] GT median(reddata[2107,*]) + 100, nbad ) $
    else $
      nbad = 0
      
    if (nbad GT 0) then begin
      if (NOT keyword_set(silent)) then $
        splog, 'WARNING: Fixing ', nbad, $
        ' dropped-pixel rows (from electronics)'
      rawdata[1:1063,ibad] = rawdata[0:1062,ibad]
      rawdata[20,ibad] = median( rawdata[20,ibad] )
      rawdata[2108:2127,ibad] = rawdata[2107:2126,ibad]
      rawdata[2107,ibad] = median( rawdata[2107,ibad] )
    endif
    
    reddata = 0  ; Free memory
  endif
  
  ;-----------
  ; Find names of the configurations files
  
  config_dir = filepath('', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
    
  if (NOT keyword_set(configfile)) then $
    configfile = findopfile('opConfig*par', mjd, config_dir, $
    /abort_notfound, silent=silent)
  if (NOT keyword_set(ecalibfile)) then $
    ecalibfile = findopfile('opECalib*par', mjd, config_dir, $
    /abort_notfound, silent=silent)
  if (NOT keyword_set(bcfile)) then $
    bcfile = findopfile('opBC*par', mjd, config_dir, $
    /abort_notfound, silent=silent)
    
  naxis1 = sxpar(hdr,'NAXIS1')
  naxis2 = sxpar(hdr,'NAXIS2')
  ;   if (naxis1 NE 2128 OR naxis2 NE 2069 AND NOT keyword_set(silent)) then $
  ;    splog, 'WARNING: Expecting 2128x2069, found '+string(naxis1)+'x'$
  ;     +string(naxis2)
  
  ;------
  ; Read in opConfig.par file
  ; Take the first entry for the configuration of each CCD in the event
  ; that there are several.
  
  yanny_read, filepath(configfile, root_dir=config_dir), pdata
  config = *pdata[0]
  yanny_free, pdata
  i = where(config.camrow EQ camrow AND config.camcol EQ camcol)
  config = config[i[0]]
  
  if (naxis1 NE config.ncols OR naxis2 NE config.nrows $
    AND NOT keyword_set(silent)) then $
    splog, 'WARNING! Config file dimensions do not match raw image'
    
  qexist = [config.amp0, config.amp1, config.amp2, config.amp3]
  
  ;------
  ; Define the "overscan" regions
  
  sover = [config.soverscan0, config.soverscan1, config.soverscan2, $
    config.soverscan3]
  nover = [config.noverscan0, config.noverscan1, config.noverscan2, $
    config.noverscan3]
    
  ;------
  ; Define the "mapped overscan" regions
    
  smapover = [config.smapoverscan0, config.smapoverscan1, $
    config.smapoverscan2, config.smapoverscan3]
  nmapover = [config.nmapoverscan0, config.nmapoverscan1, $
    config.nmapoverscan2, config.nmapoverscan3]
    
  ;------
  ; Define the "overscan rows" (at the bottom of the CCD)
    
  soverrow = [config.soverscanrows0, config.soverscanrows1, $
    config.soverscanrows2, config.soverscanrows3]
  noverrow = [config.noverscanrows0, config.noverscanrows1, $
    config.noverscanrows2, config.noverscanrows3]
    
  ;------
  ; Data position in the original image
    
  sdatarow = [config.sdatarow0, config.sdatarow1, $
    config.sdatarow2, config.sdatarow3]
  sdatacol = [config.sdatasec0, config.sdatasec1, config.sdatasec2, $
    config.sdatasec3]
  nrow = [config.ndatarow0, config.ndatarow1, $
    config.ndatarow2, config.ndatarow3]
  ncol = [config.ndatasec0, config.ndatasec1, config.ndatasec2, $
    config.ndatasec3]
    
  ;------
  ; Data position in the final (trimmed) image
    
  srow = [config.sccdrowsec0, config.sccdrowsec1, $
    config.sccdrowsec2, config.sccdrowsec3]
  scol = [config.sccdcolsec0, config.sccdcolsec1, $
    config.sccdcolsec2, config.sccdcolsec3]
    
  if (naxis2 EQ 2049) then begin
    if (NOT keyword_set(silent)) then $
      splog, 'WARNING: NROWS is 2049, adjusting config entries'
    sdatarow = sdatarow - 20
    noverrow = 1
  endif
  
  ;------
  ; Read in ECalib File
  
  yanny_read, filepath(ecalibfile, root_dir=config_dir), pdata
  ecalib = *pdata[0]
  yanny_free, pdata
  ecalib = ecalib[ where(ecalib.camrow EQ camrow AND ecalib.camcol EQ camcol) ]
  
  gain = [ecalib.gain0, ecalib.gain1, ecalib.gain2, ecalib.gain3]
  rnoise_expect = [ecalib.readnoiseDN0, ecalib.readnoiseDN1, $
    ecalib.readnoiseDN2, ecalib.readnoiseDN3]
  rnoise_measure = fltarr(4)
  fullWellDN = [ecalib.fullWellDN0, ecalib.fullWellDN1, $
    ecalib.fullWellDN2, ecalib.fullWellDN3]
    
  ;------
  ; Construct the final image
    
  igood = where(qexist)
  nr = max((srow+nrow)[igood])
  nc = max((scol+ncol)[igood])
  if ((size(image))[0] NE 2) then image = fltarr(nc, nr) $
  else if ((size(image))[1] NE nc OR (size(image))[2] NE nr OR $
  (size(image))[3] NE 4) then image = fltarr(nc, nr)
  
   ;------
   ; Test to see if the shutter was open during readout if the exposure
   ; was longer than 640 seconds.
  
   exptime = sxpar(hdr, 'EXPTIME')
   qshutter = 0

   ; Toggle the variable QSHUTTER if the observation was taken before
   ; MJD=51570 and this was not a bias or dark exposure.

   if (exptime GT 640 AND (readimg OR readivar) $
     AND mjd GT 0 AND mjd LT 51570 $
     AND flavor NE 'bias' AND flavor NE 'dark') then qshutter = 1
  
   ;------
   ; Construct IMAGE
  
  for iamp=0, 3 do begin
    if (qexist[iamp] EQ 1) then begin
      if (readimg OR readivar) then begin
        if (nover[iamp] NE 0) then begin
          ; Use the "overscan" region
          biasreg = rawdata[sover[iamp]:sover[iamp]+nover[iamp]-1, $
            sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1]
        endif else if (nmapover[iamp] NE 0) then begin
          ; Use the "mapped overscan" region
          biasreg = rawdata[smapover[iamp]:smapover[iamp]+nmapover[iamp]-1, $
            sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1]
        endif
        
        ; Dispose of the first 5 rows of the bias when computing stats,
        ; since those are the rows with transient effects (especially
        ; in the right amplifier of r1).
        ; (If we subtracted the pixel bias image first, we could
        ; probably keep these first 5 rows.)
        ; Also dispose of the first and last columns, since those
        ; sometimes show transients.
        bdimen = size(biasreg, /dimens)
        biasreg = biasreg[1:bdimen[0]-2,5:bdimen[1]-1]
        
        ; Compute statistics of the bias region that only reject
        ; the 0.5% of smallest and largest values.  This will then
        ; report a large standard deviation for very skew distributions,
        ; whereas using DJS_ITERSTAT might just clip away those values.
        ;            djs_iterstat, biasreg, sigrej=3.0, mean=biasval, sigma=readoutDN
        isort = sort(biasreg)
        nn = n_elements(biasreg)
        ii = isort[long(0.005*nn) : long(0.995*nn)]
        biasval = mean(biasreg[ii])
        ; The factor of 1.04 below is to account for clipping the
        ; lowest and highest 0.5% of all pixel values.
        rnoise_measure[iamp] = 1.04 * stddev(biasreg[ii], /double)
        
        if (NOT keyword_set(silent)) then begin
          splog, 'Measured read-noise in DN for amp#', iamp, ' = ', $
            rnoise_measure[iamp]
          splog, 'Measured bias value in DN for amp#', iamp, ' = ', biasval
          splog, 'Applying gain for amp#', iamp, ' = ', gain[iamp]
        endif
        
        ; Trigger warning message if measured read noise is > 1 ADU
        ; larger than that expected from the op files.
        if (rnoise_measure[iamp] GT rnoise_expect[iamp]+1.0) then $
          splog, 'WARNING: Amp #', iamp, $
          ' expected read noise = ', rnoise_expect[iamp], $
          ', measured = ', rnoise_measure[iamp], ' DN', $
          format='(a,i1,a,f5.2,a,f8.2,a)'
          
        ; Trigger a warning message if the bias region has differences
        ; between the 16th-percentile and 84th-percentile that are
        ; either too small or too large compared to the expected noise,
        ; or between the 2.3 and 97.7-percentile (e.g., 2-sigma),
        ; or between the 0.15 and 99.85-percentile (e.g., 3-sigma).
        for siglevel=1,3 do begin
          ; The value of BIASDIFF should be 2*readnoise*siglevel.
          fnormal = erf(sqrt(2.)/2.*siglevel) ; =0.682,0.954,0.997
          biasdiff = biasreg[isort[long(nn*(0.5+0.5*fnormal))]] $
            - biasreg[isort[long(nn*0.5*(1.-fnormal))]]
          qbiaslo = biasdiff LT 2*siglevel*(0.9 * rnoise_expect[iamp] - 1)
          qbiashi = biasdiff GT 2*siglevel*(1.1 * rnoise_expect[iamp] + 1)
          if (qbiaslo OR qbiashi) then $
            splog, 'WARNING: Amp #', iamp, ' bias region difference at ', $
            100.*fnormal, 'th-percentile = ', $
            biasdiff, ' DN, expected ', 2*siglevel*rnoise_expect[iamp], $
            format='(a,i1,a,f5.2,a,f7.1,a,f5.1)'
        endfor
        
        ; Compute the standard deviation in the bias region again,
        ; but using a weaker 5-sigma clipping.  This is done solely
        ; for the purpose of identifying any electronics problems,
        ; such as that on the left half of r2 for plate 356 on MJD 51779.
        ; Trigger a warning message if above 10 ADU.
        ; --> This is now deprecated by the above tests.
        
        ;            djs_iterstat, biasreg, sigrej=5.0, sigma=testrms
        ;            if (NOT keyword_set(silent)) then begin
        ;               if (testrms LT 10.) then $
        ;                splog, 'Std. dev. in bias region for amp#', iamp, $
        ;                 ' = ', testrms, ' DN (5-sig clip)' $
        ;               else $
        ;                splog, 'WARNING: Std. dev. in bias region for amp#', iamp, $
        ;                 ' = ', testrms, ' DN (5-sig clip)'
        ;            endif
        
        ; In the data region, not more than 3e-7 of the pixels
        ; should be more than 5-sigma below the bias level.
        ; Report a warning message if this fraction is > 0.01 percent.
        ; (Dispose of the first 5 rows of the bias when computing stats,
        ; as we do in all the above tests too.)
        
        lovalue = (biasval - 5.*rnoise_expect[iamp])
        lopixval = $
          rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
          sdatarow[iamp]+5:sdatarow[iamp]+nrow[iamp]-1] $
          LT lovalue
        numlopix = total(lopixval)
        fraclopix = numlopix / n_elements(lopixval)
        lopixval = 0 ; clear memory
        if (NOT keyword_set(silent)) then $
          splog, 'Amp #', iamp, ' Number of pix below bias-5*sigma=', $
          lovalue, ' DN is ', numlopix, $
          format='(a,i1,a,f6.0,a,i7)'
        if (100*fraclopix GT 0.01) then $
          splog, 'WARNING: Amp #', iamp, ' way too many pixels (', $
          100*fraclopix, '%) below bias-5*sigma=', lovalue, ' DN', $
          format='(a,i1,a,f6.2,a,f6.0,a)'
          
        ; Copy the data for this amplifier into the final image
        ; Subtract the bias (in DN), and then multiply by the gain
        ; Now image is in electrons
          
        image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
          srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
          (rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
          sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] - biasval) $
          * gain[iamp]
          
        ; Add to the header
        sxaddpar, hdr, 'GAIN'+string(iamp,format='(i1)'), $
          gain[iamp], ' Gain in electrons per ADU'
        sxaddpar, hdr, 'RDNOISE'+string(iamp,format='(i1)'), $
          gain[iamp]*rnoise_measure[iamp], ' Readout noise in electrons'
      endif
    endif
  endfor

  ;------
  ; If the shutter was open during readout, try to correct for that.
  ; Construct an image of the smeared light on the CCD during readout.
  ; Note that we work from IMAGE, which already has the bias removed
  ; and is multiplied by the gain.
  
  if (qshutter) then begin
    if (NOT keyword_set(silent)) then $
      splog, 'WARNING: Correcting for open shutter during readout '
      
    t1 = exptime ; Read time for entire frame
    t2 = 0.026976 ; Read time for one row of data (from Connie Rockosi)
    
    smearimg = 0 * image
    smearimg[*,0] = image[*,0]
    ny = (size(image,/dimens))[1]
    for i=1, ny-1 do begin
      ; Burles counter of row number...
      ;print, format='($, ".",i4.4,a5)', i, string([8b,8b,8b,8b,8b])
    
      smearimg[*,i] = smearimg[*,i-1] + image[*,i]
    endfor
    smearimg = (t2/t1) * smearimg
    image = image - smearimg
    
    if (NOT keyword_set(silent)) then begin
      splog, 'Median value of open-shutter contamination = ', median(smearimg)
      splog, 'Max value of open-shutter contamination = ', max(smearimg)
    endif
  endif
  
  ;------
  ; Construct INVVAR
  
  if (readivar) then begin
    if ((size(invvar))[0] NE 2) then $
      invvar = fltarr(nc, nr) $
    else if ((size(invvar))[1] NE nc OR (size(invvar))[2] NE nr OR $
    (size(invvar))[3] NE 4) then $
    invvar = fltarr(nc, nr)
    
  ;------
  ; SATMASK = Mask for saturated the detector, 0=bad
  ; ADMASK = Mask for saturating the A/D converter (at 65535), 1=bad
  ; BCMASK = Mask for bad columns, 1=bad
    
  satmask = bytarr(nc, nr)
  admask = bytarr(nc, nr)
  bcmask = bytarr(nc, nr)
  
  for iamp=0, 3 do begin
    if (qexist[iamp] EQ 1) then begin
    
      satmask[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
        srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
        rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
        sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] LT $
        fullWellDN[iamp]
        
      admask[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
        srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
        rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
        sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] EQ 65535
        
      ; Flux term below
      expr1 = abs(image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
        srow[iamp]:srow[iamp]+nrow[iamp]-1])
      ; Add to the variance image from the open shutter
      ; by adding 1% of that signal^2 to the variance.
      ; This says that the uncertainty in this subtracted
      ; quantity is about 10%.
      if (qshutter) then $
        expr2 = abs(smearimg[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
        srow[iamp]:srow[iamp]+nrow[iamp]-1]) $
      else expr2 = 0
      ; Read noise term below
      expr3 = (rnoise_measure[iamp]*gain[iamp])^2
      ; Term below to limit best S/N to under 100
      expr4 = 1.e-4 * expr1^2
      
      ; Horrible kludge for this night MJD=51779 (22/23 Aug 2000).
      ; Add a noise term of 100 ADU to the left amplifier of r2
      if (mjd EQ 51779 AND camname EQ 'r2' AND iamp EQ 2) then $
        expr4 = expr4 + 100.^2
        
      invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
        srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
        1.0/(expr1 + expr2 + 0.01 * expr2^2 + expr3 + expr4)
        
      if (keyword_set(killdata)) then $
        invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
        srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
        invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
        srow[iamp]:srow[iamp]+nrow[iamp]-1] * $
        killdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
        sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1]
    endif
  endfor
  
  ;------
  ; Look for blead trails and mask them.
  ; At present, SATMASK is set to 0 for saturated pixels.  Look for any
  ; column with >=11 saturated pixels in a row.  In that case, mask
  ; all pixels from that row number up to the top of the CCD.
  
  kern = transpose(fltarr(11) + 1.0)
  mask1 = fix( convol(satmask+0.0, kern, /center, /edge_truncate) ) EQ 0
  ; 1=bad
  qblead = total(mask1, 2) GT 0 ; =1 for each column with a blead trail
  iblead = where(qblead, nblead)
  sxaddpar, hdr, 'NBLEAD', nblead, ' Number of columns with blead trails'
  if (nblead GT 0) then begin
    if (NOT keyword_set(silent)) then $
      splog, 'Number of bleading columns = ', nblead
    for i=0, nblead-1 do begin
      icol = iblead[i] ; Column number for this blead trail
      irow = (where(mask1[icol,*]))[0] ; First bad row in this column
      satmask[icol,irow:nr-1] = 0
    endfor
  endif
  
  ;------
  ; Count the number of rows (wavelengths) that we think are saturated.
  ; A row is considered to be saturated if at least 20 of its pixels are.
  ; Note that this counting is done before masking out bad columns.
  ; What we really should do is ignore bad columns when computing this.
  
  if (arg_present(nsatrow)) then begin
    totsat = total((1-satmask), 1)
    junk = where(totsat GE 20, nsatrow)
    if (NOT keyword_set(silent)) then $
      splog, 'Number of saturated rows = ', nsatrow
  endif
  
  ;------
  ; Mask out pixels that saturated the A/D converter, plus mask
  ; all neighbors within 1 pixel
  
  ngrow = 1
  width = 2*ngrow + 1
  admask = smooth(admask * width^2, width) GT 0 ; 1=bad
  
  ;------
  ; Mask out bad columns
  
  bcmask = make_badcolumn_mask( $
    filepath(bcfile,root_dir=config_dir), camrow, camcol, silent=silent, $
    nr=nr,nc=nc)
    
  ;------
  ; For masked pixels, set INVVAR=0
    
  invvar = invvar * satmask * (1-admask) * (1-bcmask)
  
  ;------
  ; Count the fraction of bad pixels, not including bad columns
  
  if (arg_present(fbadpix)) then begin
    junk = where((satmask EQ 0 OR admask EQ 1) AND (bcmask EQ 0), njunk)
    fbadpix = float(njunk) / (float(nc) * float(nr))
  endif
endif ; endif from what statement ???

satmask = 0 ; clear memory
admask = 0 ; clear memory
; bcmask = 0 ; clear memory

; Generate the CCD pixel mask (all zeros for SDSS-I data)
if (arg_present(ccdmask) AND readimg) then $
 ccdmask = lonarr(size(image,/dimens))

endelse ; End toggle between reading BOSS or SDSS-I images

;---------------------------------------------------------------------------
; Correct image with bias image
; (If BOSSGOOD then this has already been done -- ASB 2011aug.)
;---------------------------------------------------------------------------

if (keyword_set(applybias) AND readimg AND (bossgood eq 0)) then begin
  pp = filepath('', root_dir=getenv('SPECFLAT_DIR'), subdirectory='biases')
  ; First search for files "pixbiasave-*.fits*", and if not found then
  ; look for "pixbias-*.fits".
  pixbiasname = findopfile('pixbiasave-*-'+camname+'.fits*', mjd, pp, $
    silent=silent)
  if (NOT keyword_set(pixbiasname)) then $
    pixbiasname = findopfile('pixbias-*-'+camname+'.fits*', mjd, pp, $
    silent=silent)
    
  if (NOT keyword_set(pixbiasname)) then begin
    if (NOT keyword_set(silent)) then $
      splog, 'WARNING: Bias image not found for camera ' + camname
  endif else begin
    if (NOT keyword_set(silent)) then $
      splog, 'Correcting with bias image ' + pixbiasname
    pixbiasimg = mrdfits(djs_filepath(pixbiasname, root_dir=pp), $
      silent=silent)
      
    if (total(size(pixbiasimg,/dimens) NE size(image,/dimens)) GT 0) then $
     message, 'Dimensions of image and bias image differ!'

    image -= pixbiasimg
    pixbiasimg = 0 ; clear memory
    
    ; Add pixbiasname to header since it has just been applied
    sxaddpar, hdr, 'PIXBIAS', pixbiasname
  endelse
endif

;---------------------------------------------------------------------------
; Correct image with pixel-to-pixel flat-field
;---------------------------------------------------------------------------

if (keyword_set(applypixflat) AND (readimg OR readivar)) then begin
  pp = filepath('', root_dir=getenv('SPECFLAT_DIR'), subdirectory='flats')
  ; First search for files "pixflatave-*.fits*", and if not found then
  ; look for "pixflat-*.fits*".
  pixflatname = findopfile('pixflatave-*-'+camname+'.fits*', mjd, pp, $
    silent=silent)
  if (NOT keyword_set(pixflatname)) then $
    pixflatname = findopfile('pixflat-*-'+camname+'.fits*', mjd, pp, $
    silent=silent)
   
  if (NOT keyword_set(pixflatname)) then begin
    if (NOT keyword_set(silent)) then $
      splog, 'WARNING: Pixel flat not found for camera ' + camname
  endif else begin
    if (NOT keyword_set(silent)) then $
      splog, 'Correcting with pixel flat ' + pixflatname
      
    pixflatimg = mrdfits(djs_filepath(pixflatname, root_dir=pp), $
     /fscale, silent=silent)
    if (total(size(pixflatimg,/dimens) NE size(image,/dimens)) GT 0) then $
     message, 'Dimensions of image and pixel flat differ!'

    ; now get bad pixel mask
    badpixelname = findopfile('badpixels-*-'+camname+'.fits*', mjd, pp, $
     silent=silent)
   
    if (NOT keyword_set(badpixelname)) then begin
      if (NOT keyword_set(silent)) then $
       splog, 'WARNING: Badpixels not found for camera ' + camname
    endif else begin
      if (NOT keyword_set(silent)) then $
      splog, 'Correcting with badpixels ' + badpixelname
  
      badpixelimg = mrdfits(djs_filepath(badpixelname, root_dir=pp), $
       /fscale, silent=silent)
      if (total(size(badpixelimg,/dimens) NE size(image,/dimens)) GT 0) then $
       message, 'Dimensions of image and badpixels differ!'
    
      ; include badpixels into pixflat
      badpixuse=where(badpixelimg ne 0,ct)
      if (ct ne 0) then pixflatimg[badpixuse]=0.0
    endelse

    if (readimg) then image /= (pixflatimg + (pixflatimg LE 0))
    if (NOT keyword_set(minflat)) then minflat = 0.0
    if (NOT keyword_set(maxflat)) then maxflat = 1.0e10
    if (readivar) then $
      invvar = invvar * pixflatimg^2 * (pixflatimg GT minflat) $
      * (pixflatimg LT maxflat)
    pixflatimg = 0 ; clear memory
    badpixelimg = 0 ; clear memory

    ; Add pixflatname to header since it has just been applied
    sxaddpar, hdr, 'PIXFLAT', pixflatname
    if (keyword_set(badpixelname)) then $
     sxaddpar, hdr, 'BADPIXEL', badpixelname

  endelse
endif

;---------------------------------------------------------------------------
; Check for NaN's
;---------------------------------------------------------------------------

; This should never happen, but just in case...

if (readimg OR readivar) then begin
  inan = where(finite(image) EQ 0, nnan)
  if (nnan GT 0) then begin
    if (NOT keyword_set(silent)) then $
      splog, 'WARNING: Replacing ', nnan, ' NaN values'
    image[inan] = 0
    invvar[inan] = 0
  endif
endif

;---------------------------------------------------------------------------
; Write output files
;---------------------------------------------------------------------------

if keyword_set(bcfile) then sxaddpar, hdr, 'OPBC', bcfile
if keyword_set(configfile) then sxaddpar, hdr, 'OPCONFIG', configfile
if keyword_set(ecalibfile) then sxaddpar, hdr, 'OPECALIB', ecalibfile
sxdelpar, hdr, 'UNSIGNED'

if (keyword_set(outfile)) then begin
  mwrfits, image, outfile, hdr, /create
  mwrfits, invvar, outfile
endif

return
end
;------------------------------------------------------------------------------
