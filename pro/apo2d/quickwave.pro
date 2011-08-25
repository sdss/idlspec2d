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
;     [ radius=, /doplot, /do_lock ] )
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
;   do_lock    - Optional keyword for SDSSPROC
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
;   apo_checklimits()
;   extract_boxcar()
;   fiberflat()
;   fileandpath()
;   findfile()
;   fitarcimage
;   fitdispersion()
;   qbadarc
;   mrdfits()
;   mwrfits
;   reject_arc()
;   rmfile
;   sdssproc
;   traceset2xy
;
; REVISION HISTORY:
;   3-Apr-2000  Written by S. Burles & D. Schlegel, APO
;-
;------------------------------------------------------------------------------
function quickwave, arcname, tsetfile, wsetfile, fflatfile, radius=radius, $
 doplot=doplot, do_lock=do_lock

   if (n_elements(arcname) NE 1) then return, 0
   if (n_elements(wsetfile) NE 1) then return, 0
   if (n_elements(tsetfile) NE 1) then return, 0
   if (NOT keyword_set(radius)) then radius = 3.0

   ;----------
   ; Dispose of any pre-existing wsetfile for this exposure.

   if (keyword_set(findfile(wsetfile))) then rmfile, wsetfile

   ;----------
   ; Read in image

   sdssproc, arcname, arcimg, hdr=archdr, color=color, camname=camname, $
    nsatrow=nsatrow, fbadpix=fbadpix, do_lock=do_lock

   configuration=obj_new('configuration', sxpar(archdr, 'MJD'))

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

   dims = size(xcen, /dimens)
   nrow = dims[0]
   ntrace = dims[1]

   ;----------
   ; Boxcar extract the arc

   flux = extract_boxcar(arcimg, xcen, radius=radius)

   ;----------
   ; Estimate inverse variance

   fluxivar = 1.0 / (abs(flux) + 10.0)

   ;----------
   ; Now find the wavelength solution

   arccoeff = configuration->spcalib_arccoeff()

   fitarcimage, flux, fluxivar, xpeak, ypeak, wset, ncoeff=arccoeff, $
     aset=aset, color=color, fibermask=fibermask, $ ; ?? maxdev=4.d-5, $
     bestcorr=bestcorr, lambda=lambda, xdif_tset=xdif_tset, $
     acoeff=configuration->spcalib_arcfitguess_acoeff(color), $
     dcoeff=configuration->spcalib_arcfitguess_dcoeff(color), $
     wrange=configuration->spcalib_fitarcimage_wrange(color)

   if keyword_set(doplot) then qaplot_arcline, xdif_tset, wset, lambda, $
     color=color, title=' Arcline Fit for '+arcname

   if (NOT keyword_set(wset)) then return, 0
   traceset2xy, wset, xx, loglam

   ;----------
   ; Fit the wavelength dispersion, and trigger warnings if the spectrographs
   ; look out-of-focus in the wavelength dimension.

   nfitcoeff = configuration->spcalib_ncoeff(color)
   nx = (size(arcimg,/dimens))[0]
   dispset = fitdispersion(flux, fluxivar, xpeak, $
    sigma=configuration->spcalib_sigmaguess(), ncoeff=nfitcoeff, $
    xmin=0.0, xmax=nx-1., $
    medwidth=medwidth, numbundles=ntrace/20, /quick) ; Hard-wires 20 fibers/bundle???

   if (apo_checklimits('arc', 'WSIGMA', camname, max(medwidth)) $
    EQ 'red') then $
    splog, 'WARNING: Median wavelength widths = ' $
    + string(medwidth,format='(4f5.2)') + ' pix (Left Bottom Top Right)';quadrupole
   ;----------
   ; Compute fiber-to-fiber flat-field variations.
   ; First see if we've already done this.

   fflatexist = keyword_set(findfile(fflatfile))
;   if (NOT fflatexist) then begin
      flat_flux = mrdfits(tsetfile,0)
      flat_ivar = mrdfits(tsetfile,1)
      fflat = fiberflat(flat_flux, flat_ivar, wset, fibermask=fibermask, $
       /dospline, $
       badflatfracthresh=configuration->spcalib_fiberflat_badflatfracthresh(),$
       minval=configuration->spcalib_fiberflat_minval(flux))
      flat_flux = 0 ; clear memory
      flat_ivar = 0 ; clear memory
      mwrfits, fflat, fflatfile, /create
      mwrfits, fibermask, fflatfile
      fflat = 0 ; clear memory
;   endif

   obj_destroy,configuration

   ;----------
   ; Write out wavelength solution

   if (sxpar(archdr,'quality') EQ 'excellent') then begin
      mwrfits, wset, wsetfile, /create
      mwrfits, flux, wsetfile
      mwrfits, fluxivar, wsetfile
      mwrfits, xpeak, wsetfile
      mwrfits, ypeak, wsetfile
   endif else begin
      splog, 'Quality is not excellent - do not write wsetfile'
   endelse

   ; Compute the wavelength between the central 2 fibers
   ; at their central 2 pixels.
   wavemid = mean( 10.^loglam[nrow/2-1:nrow/2,ntrace/2-1:ntrace/2] )

   nlamps = (size(xpeak,/dimens))[1]
   rstruct = create_struct('WSETFILE', fileandpath(wsetfile), $
                           'WAVEMID', double(wavemid), $
                           'BESTCORR', float(bestcorr), $
                           'NLAMPS', long(nlamps), $
                           'WSIGMA_QUADRANT', float(medwidth), $
                           'WSIGMA', float(max(medwidth)) )

   return, rstruct

end 
;------------------------------------------------------------------------------
