;+
; NAME:
;   spgain
;
; PURPOSE:
;   Measure the gain + read noise for SDSS spectroscopic images
;
; CALLING SEQUENCE:
;   spgain, flatfile1, flatfile2, [ biasfile1, biasfile2, indir=, $
;    xskip=, yskip=, /simulate, gain=, rnoise=]
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
 xskip=xskip, yskip=yskip, gain=gain, rnoise=rnoise, $
 simulate=simulate

   if (N_params() NE 2 AND N_params() NE 4) then $
    message, 'Must specify either 2 or 4 file names'

   if (NOT keyword_set(xskip)) then xskip = 50
   if (NOT keyword_set(yskip)) then yskip = 5

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
   ximg = djs_laxisgen([nx, ny], iaxis=0)
   ampimg = 1 + (ximg GT nx/2)
   ampimg[0:xskip-1,*] = 0
   ampimg[nx-xskip:nx-1,*] = 0
   ampimg[*,0:yskip-1] = 0
   ampimg[*,ny-yskip:ny-1] = 0

   nloop = 100L
   ngoodpix = lonarr(nloop,2)
   gainarr = fltarr(nloop,2)
   rnoisearr = fltarr(nloop,2)

   ; Loop over different count levels
   for iloop=0, nloop-1 do begin
      for iamp=0, 1 do begin
         flux1 = (iloop+1) * 200. + 1000.
         flux2 = (iloop+2) * 200. + 1000.
         indx = where(flatimg1 GE flux1 AND flatimg1 LT flux2 $
          AND ampimg EQ iamp+1, ct)
         ngoodpix[iloop,iamp] = ct

         if (ct GT 1000) then begin ; At least 1000 pixels...

            ; Correction factor for measured sigmas...
            corrfac = sqrt(ct / (ct-1.))

            flatsub1 = flatimg1[indx]
            flatsub2 = flatimg2[indx]

            ; Compute statistics for flats
            djs_iterstat, flatsub1, sigrej=sigrej, mean=flatmean1
            djs_iterstat, flatsub2, sigrej=sigrej, mean=flatmean2
            djs_iterstat, flatsub2 - flatsub1, sigrej=sigrej, $
             sigma=flatdifsig
            flatdifsig = corrfac * flatdifsig

            ; Compute statistics for biases
            if (N_params() EQ 4) then begin
               biassub1 = biasimg1[indx]
               biassub2 = biasimg2[indx]

               djs_iterstat, biassub1, sigrej=sigrej, mean=biasmean1
               djs_iterstat, biassub2, sigrej=sigrej, mean=biasmean2
               djs_iterstat, biassub2 - biassub1, sigrej=sigrej, $
                sigma=biasdifsig
               biasdifsig = corrfac * biasdifsig
            endif

            gainarr[iloop,iamp] = $
             (flatmean1 + flatmean2 - biasmean1 - biasmean2) / $
             (flatdifsig^2 - biasdifsig^2)
            rnoisearr[iloop,iamp] = biasdifsig / sqrt(2.)

         endif
      endfor
   endfor

   ; Compute the median gain + read noise for each amplifier
   for iamp=0, namp-1 do begin
      ii = where(gainarr[*,iamp] GT 0)
      djs_iterstat, gainarr[ii,iamp], median=med1, sigma=sig1
      gain[iamp] = med1
      gain_rms[iamp] = sig1
      if (N_params() EQ 4) then begin
         djs_iterstat, rnoisearr[ii,iamp], median=med1, sigma=sig1
         rnoise[iamp] = med1
         rnoise_rms[iamp] = sig1
      endif
      splog, 'Amplifier #', iamp, ' Gain=', gain[iamp], ' +/-', $
       gain_rms[iamp], ' Rnoise=', rnoise[iamp], ' +/-', rnoise_rms[iamp], ' DN'
   endfor

   return
end
;------------------------------------------------------------------------------
