;+
; NAME:
;   quickwave
;
; PURPOSE:
;   Perform full wavelength calibration with boxcar extraction
;    requires trace from a previously output tsetfile
;
; CALLING SEQUENCE:
;   rstruct = quickwave( arcname, tsetfile, wsetfile, fflatfile, $
;       radius=radius, doplot=doplot )
;
; INPUTS:
;   arcname    - Raw SDSS Arclamp image to be processed
;   tsetfile   - Name of fits file which contains matched trace
;   wsetfile   - Name of fits file which will contain wavelength solution
;   fflatfile  - Name of fits file which will contain flat field vectors
;
; OPTIONAL INPUTS:
;   radius     - Radius for boxcar extraction (default 3)
;   doplot     - Used for debugging purposes
;
; OUTPUT:
;   rstruct    - Results to be added html file upon completion
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   extract_boxcar()
;   fiberflat()
;   fileandpath()
;   fitarcimage
;   qbadarc
;   mrdfits()
;   mwrfits
;   reject_arc()
;   sdssproc
;   traceset2xy
;
; REVISION HISTORY:
;   3-Apr-2000  Written by S. Burles & D. Schlegel, APO
;-
;------------------------------------------------------------------------------
function quickwave, arcname, tsetfile, wsetfile, fflatfile, radius=radius, $
    doplot=doplot

   if (n_elements(arcname) NE 1) then return, 0
   if (n_elements(wsetfile) NE 1) then return, 0
   if (n_elements(tsetfile) NE 1) then return, 0
   if (NOT keyword_set(radius)) then radius = 3.0

   ;----------
   ; Read in image

   sdssproc, arcname, arcimg, hdr=archdr, color=color, camname=camname, $
    nsatrow=nsatrow, fbadpix=fbadpix

   ;-----
   ; Decide if this arc is bad

   qbadarc = reject_arc(arcimg, archdr, nsatrow=nsatrow, fbadpix=fbadpix)
   if (qbadarc) then begin
      splog, 'ABORT: Unable to reduce arc'
      return, 0
   endif

   ;----------
   ; Read in the reduced flat

   tset = mrdfits(tsetfile,2)
   fibermask = mrdfits(tsetfile,4)
   traceset2xy, tset, ycen, xcen 

   ;----------
   ; Boxcar extract the arc

   flux = extract_boxcar(arcimg, xcen, radius=radius)

   ;----------
   ; Estimate inverse variance

   fluxivar = 1.0 / (abs(flux) + 10.0)

   ;----------
   ; Now find the wavelength solution

   fitarcimage, flux, fluxivar, xpeak, ypeak, wset, aset=aset, $
     fibermask=fibermask, bestcorr=bestcorr, lambda=lambda, $
     color=color, maxdev=2.0e-5, xdif_tset=xdif_tset

   if keyword_set(doplot) then qaplot_arcline, xdif_tset, wset, lambda, $
     color=color, title=' Arcline Fit for '+arcname

   if (NOT keyword_set(wset)) then return, 0
   traceset2xy, wset, xx, loglam

   ;----------
   ; Fit the wavelength dispersion, and trigger warnings if the spectrographs
   ; look out-of-focus in the wavelength dimension.

   nfitcoeff = color EQ 'red' ? 4 : 3
   dispset = fitdispersion(flux, fluxivar, xpeak, $
    sigma=1.0, ncoeff=nfitcoeff, xmin=0.0, xmax=2047.0, medwidth=medwidth)

   ;----------
   ; Compute fiber-to-fiber flat-field variations.
   ; First see if we've already done this.

   fflatexist = keyword_set(findfile(fflatfile))
   if (NOT fflatexist) then begin
      flat_flux = mrdfits(tsetfile,0)
      flat_ivar = mrdfits(tsetfile,1)
      fflat = fiberflat(flat_flux, flat_ivar, wset, fibermask=fibermask, $
       /dospline)
      flat_flux = 0 ; clear memory
      flat_ivar = 0 ; clear memory
      mwrfits, fflat, fflatfile, /create
      mwrfits, fibermask, fflatfile
      fflat = 0 ; clear memory
   endif

   ;----------
   ; Write out wavelength solution

   mwrfits, wset, wsetfile, /create

   wavemin = 10^(min(loglam))
   wavemax = 10^(max(loglam))
   nlamps = (size(xpeak,/dimens))[1]
   rstruct = create_struct('WSETFILE', fileandpath(wsetfile), $
                           'WAVEMIN', wavemin, $
                           'WAVEMAX', wavemax, $
                           'BESTCORR', bestcorr, $
                           'NLAMPS', nlamps, $
                           'DISPSIGMA_QUADRANT', medwidth, $
                           'DISPSIGMA', max(medwidth) )

   return, rstruct
end
;------------------------------------------------------------------------------
