;+
; NAME:
;   spgain
;
; PURPOSE:
;   Measure the gain + read noise for SDSS spectroscopic images
;
; CALLING SEQUENCE:
;   spgain, flatfile1, flatfile2, [ biasfile1, biasfile2, indir=, $
;    xskip=, yskip=, xave, yave=, /simulate, gain=, rnoise=]
;
; INPUTS:
;   flatfile1  - File name for flat #1
;   flatfile2  - File name for flat #2
;   biasfile1  - File name for bias #1
;   biasfile2  - File name for bias #2
;
; OPTIONAL KEYWORDS:
;   indir      - Input directory for all files
;   xskip      - Number of columns to ignore at beginning and end of each
;                amplifier; default to 50
;   yskip      - Number of rows to ignore at beginning and end of each
;                amplifier; default to 5
;   xave       - Number of columns to analyze for each calculation of gain and
;                read noise; default to 100
;   yave       - Number of rows to analyze for each calculation of gain and
;                read noise; default to 10
;   simulate   - If set, then replace the images with simulated images with
;                a gain of 1.3 e-/ADU and read noise of 3.5 ADU.
;
; OUTPUTS:
;   gain       - Gain in electrons/ADU for each amplifier (array)
;   rnoise     - Read noise in electrons for each amplifier (array)
;
; COMMENTS:
;
; EXAMPLES:
;   Compute the read nosie and gain using the bias frames and flats
;   taken for this very purpose on MJD 53114 (18/19 April 2004).
;   For the b1 CCD:
;     IDL> spgain, 'sdR-r2-00026145.fit.gz', 'sdR-r2-00026146.fit.gz', $
;          'sdR-r2-00026147.fit.gz', 'sdR-r2-00026148.fit.gz'
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_iterstat
;   sdssproc
;
; REVISION HISTORY:
;   21-Nov-1999  Written by D. Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro spgain, flatfile1, flatfile2, biasfile1, biasfile2, indir=indir, $
 xskip=xskip, yskip=yskip, xave=xave, yave=yave, gain=gain, rnoise=rnoise, $
 simulate=simulate

   if (N_params() NE 2 AND N_params() NE 4) then $
    message, 'Must specify either 2 or 4 file names'

   if (NOT keyword_set(xskip)) then xskip = 50
   if (NOT keyword_set(yskip)) then yskip = 5
   if (NOT keyword_set(xave)) then xave = 100
   if (NOT keyword_set(yave)) then yave = 10

   namp = 2 ; 2 amplifiers

   gain = fltarr(namp)
   rnoise = fltarr(namp)
   gain_rms = fltarr(namp)
   rnoise_rms = fltarr(namp)

   ; Read flats
   sdssproc, flatfile1, flatimg1, indir=indir, hdr=hdr
   sdssproc, flatfile2, flatimg2, indir=indir

   ; Read biases
   if (N_params() EQ 4) then begin
      sdssproc, biasfile1, biasimg1, indir=indir
      sdssproc, biasfile2, biasimg2, indir=indir
   endif else begin
      biasmean1 = 0
      biasmean2 = 0
      biasdifsig = 0
   endelse

   ; Take out the gain so that all images are still in ADU...
   config_dir = filepath('', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
   ecalibfile = findopfile('opECalib*par', sxpar(hdr,'MJD'), config_dir, $
    /abort_notfound, /silent)
   ecalib = yanny_readone(filepath(ecalibfile, root_dir=config_dir))
   i = where(ecalib.camrow EQ sxpar(hdr,'CAMROW') $
    AND ecalib.camcol EQ sxpar(hdr,'CAMCOL'))
   gain_input = [fltarr(1024,2048)+ecalib[i].gain2, $
    fltarr(1024,2048)+ecalib[i].gain3]
   flatimg1 = flatimg1 / gain_input
   flatimg2 = flatimg2 / gain_input
   if (N_params() EQ 4) then begin
      biasimg1 = biasimg1 / gain_input
      biasimg2 = biasimg2 / gain_input
   endif

   if (keyword_set(simulate)) then begin
      simgain = 1.3
      simnoise = 3.5
      fakeimg = flatimg1 > 0
;      fakeimg = smooth(randomu(123456,2048,2048) * 100., 5)
;      fakeimg = fltarr(2048,2048) + 1000.
      flatimg1= fakeimg + randomn(123,2048,2048) * simnoise $
       + randomn(987,2048,2048) * sqrt(fakeimg) / sqrt(simgain)
      flatimg2= fakeimg + randomn(234,2048,2048) * simnoise $
       + randomn(888,2048,2048) * sqrt(fakeimg) / sqrt(simgain)
      biasimg1 = randomn(345,2048,2048) * simnoise
      biasimg2 = randomn(456,2048,2048) * simnoise
   endif

   dims = size(flatimg1, /dimens)
   nx = dims[0]
   ny = dims[1]

   xstart = [0, 1024] + xskip
   xend = [1023, 2047] - xskip

   nxblock = fix((nx-2*xskip)/xave)
   nyblock = fix(ny/yave)
   gainarr = fltarr(nxblock,nyblock)
   rnoisearr = fltarr(nxblock,nyblock)
   corrfac = sqrt(xave*yave / (xave*yave-1.)) ; Correction factor for measured sigmas

   for ix=0, nxblock-1 do begin
   for iy=0, nyblock-1 do begin

      x1 = ix * xave + xskip
      x2 = x1 + xave - 1
      y1 = iy * yave
      y2 = y1 + yave - 1

      flatsub1 = flatimg1[x1:x2,y1:y2]
      flatsub2 = flatimg2[x1:x2,y1:y2]

      ; Compute statistics for flats
      djs_iterstat, flatsub1, sigrej=sigrej, maxiter=maxiter, $
       mean=flatmean1
      djs_iterstat, flatsub2, sigrej=sigrej, maxiter=maxiter, $
       mean=flatmean2
      djs_iterstat, flatsub2 - flatsub1, sigrej=sigrej, maxiter=maxiter, $
       sigma=flatdifsig
      flatdifsig = corrfac * flatdifsig

      ; Compute statistics for biases
      if (N_params() EQ 4) then begin
         biassub1 = biasimg1[x1:x2,y1:y2]
         biassub2 = biasimg2[x1:x2,y1:y2]

         djs_iterstat, biassub1, sigrej=sigrej, maxiter=maxiter, $
          mean=biasmean1
         djs_iterstat, biassub2, sigrej=sigrej, maxiter=maxiter, $
          mean=biasmean2
         djs_iterstat, biassub2 - biassub1, sigrej=sigrej, maxiter=maxiter, $
          sigma=biasdifsig
         biasdifsig = corrfac * biasdifsig
      endif

      gainarr[ix,iy] = (flatmean1 + flatmean2 - biasmean1 - biasmean2) / $
       (flatdifsig^2 - biasdifsig^2)
      rnoisearr[ix,iy] = biasdifsig / sqrt(2.)

   endfor
   ; Burles counter of row number...
   print, format='($, ".",i4.4,a5)', ix, string([8b,8b,8b,8b,8b])
   endfor

   ; Compute the median gain + read noise for each amplifier
   for iamp=0, namp-1 do begin
      djs_iterstat, gainarr[nxblock/2*[iamp,iamp+1]+[0,-1],*], $
       median=med1, sigma=sig1
      gain[iamp] = med1
      gain_rms[iamp] = sig1
      if (N_params() EQ 4) then begin
         djs_iterstat, rnoisearr[nxblock/2*[iamp,iamp+1]+[0,-1],*], $
          median=med1, sigma=sig1
         rnoise[iamp] = med1
         rnoise_rms[iamp] = sig1
      endif
      splog, 'Amplifier #', iamp, ' Gain=', gain[iamp], ' +/-', $
       gain_rms[iamp], ' Rnoise=', rnoise[iamp], ' +/-', rnoise_rms[iamp], ' DN'
   endfor

   return
end
;------------------------------------------------------------------------------
