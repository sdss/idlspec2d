; spgain,'sdR-01-00001367.fit','sdR-01-00001369.fit', $
;  gain=gain,rnoise=rnoise, yave=9

;+
; NAME:
;   spgain
;
; PURPOSE:
;   Measure the gain + read noise for SDSS spectroscopic images
;
; CALLING SEQUENCE:
;   spgain, flatfile1, flatfile2, [ biasfile1, biasfile2, indir=indir, $
;    xskip=xskip, yave=yave, gain=gain, rnoise=rnoise ]
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
;   yave       - Number of rows to analyze for each calculation of gain and
;                read noise; default to 15
;
; OUTPUTS:
;   gain       - Gain in electrons/ADU for each amplifier (array)
;   rnoise     - Read noise in electrons for each amplifier (array)
;
; COMMENTS:
;
; BUGS:;
;
; PROCEDURES CALLED:
;   djs_iterstat
;   sdssproc
;
;; REVISION HISTORY:
;   21-Nov-1999  Written by D. Schlegel, Princeton.
;-
;------------------------------------------------------------------------------

pro spgain, flatfile1, flatfile2, biasfile1, biasfile2, indir=indir, $
 xskip=xskip, yave=yave, gain=gain, rnoise=rnoise

   if (N_params() NE 2 AND N_params() NE 4) then $
    message, 'Must specify either 2 or 4 file names'

   if (NOT keyword_set(xskip)) then xskip = 50
   if (NOT keyword_set(yave)) then yave = 15

   namp = 2 ; 2 amplifiers

   gain = fltarr(namp)
   rnoise = fltarr(namp)

   ; Read flats
   sdssproc, flatfile1, flatimg1, indir=indir
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

   dims = size(flatimg1, /dimens)
   nx = dims[0]
   ny = dims[1]

   xstart = [0, 1024] + xskip
   xend = [1023, 2047] - xskip

   nyblock = fix(ny/yave)
   gainarr = fltarr(nx,nyblock)
   rnoisearr = fltarr(nx,nyblock)
   corrfac = sqrt(yave / (yave-1)) ; Correction factor for measured sigmas

   for ix=0, nx-1 do begin
   for iy=0, nyblock-1 do begin

      y1 = iy * yave
      y2 = y1 + yave - 1

      flatsub1 = flatimg1[ix,y1:y2]
      flatsub2 = flatimg2[ix,y1:y2]

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
         biassub1 = biasimg1[ix,y1:y2]
         biassub2 = biasimg2[ix,y1:y2]

         djs_iterstat, biassub1, sigrej=sigrej, maxiter=maxiter, $
          mean=biasmean1
         djs_iterstat, biassub2, sigrej=sigrej, maxiter=maxiter, $
          mean=biasmean2
         djs_iterstat, biassub2 - biassub1, sigrej=sigrej, maxiter=maxiter, $
          sigma=biasdifsig
         biasdifsig = biasfac * flatdifsig
      endif

      gainarr[ix,iy] = (flatmean1 + flatmean2 - biasmean1 - biasmean2) / $
       (flatdifsig^2 - biasdifsig^2)
      rnoisearr[ix,iy] = gainarr[ix,iy] * biasdifsig / sqrt(2.)

   endfor
   ; Burles counter of column number...
   print, format='($, ".",i4.4,a5)', ix, string([8b,8b,8b,8b,8b])
   endfor

   ; Compute the median gain + read noise for each amplifier
   for iamp=0, namp-1 do begin
      gain[iamp] = median( gainarr[xstart[iamp]:xend[iamp],*] )
      if (N_params() EQ 4) then $
       rnoise[iamp] = median( rnoisenarr[xstart[iamp]:xend[iamp],*] )
   endfor

stop
   return
end
;------------------------------------------------------------------------------
