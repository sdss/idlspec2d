;+
; NAME:
;   spcoadd_v5
;
; PURPOSE:
;   Combine several reduced frames of the same objects
;
; CALLING SEQUENCE:
;   spcoadd_v5, spframes, outputname, $
;    [ mjd=, binsz=, zeropoint=, nord=, wavemin=, $
;    bkptbin=, window=, maxsep=, adderr=, plotsnfile=, combinedir= ]
;
; INPUTS:
;   spframes       - Name(s) of files to combine (written by SPREDUCE)
;   outputname     - Output file name
;
; OPTIONAL KEYWORDS:
;   mjd            - The MJD to put in the output header
;   binsz          - Bin size (in log-10 wavelength) in output spectra;
;                    default to 1d-4, which corresponds to 69.02977415 km/s.
;   zeropoint      - Log10(lambda) zero-point of the output spectra;
;                    the output wavelength bins are chosen such that
;                    one bin falls exactly on this value;
;                    default to 3.5D, which corresponds to 3162.27766 Ang.
;   nord           - Order for spline fit; default to 3 (cubic spline).
;   wavemin        - Log-10 wavelength of first pixel in output spectra;
;                    default to the nearest bin to the smallest wavelength
;                    of the input spectra.
;   bkptbin        - ???
;   window         - Window size for apodizing the errors of the spectrum
;                    from each individual frame;
;                    default to 100 pixels apodization on each end of the
;                    spectra.
;   maxsep         - ???
;   adderr         - Additional error to add to the formal errors, as a
;                    fraction of the flux.
;   combinedir     - Optional output directory
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine can combine data from multiple (different) plug maps.
;   Objects are matched based upon their positions agreeing to 2 arc sec.
;
;   All input files must have the same number of pixels per spectrum,
;   i.e. 2048 wavelength samplings, although those wavelengths can
;   be different.
;
;   The input files (FILENAMES) have their pixelmasks modified by this routine.
;
;   Flux-correction files are also read in, where they are assumed to
;   have the name spFluxcorr-EEEEEEEE-S.fits, where EEEEEEEE is the exposure
;   number and S is the spectrograph ID (1 or 2).
;
; EXAMPLES:
;
; BUGS:
;   Should only apodize starting with the first/last good pixel of a spectrum.
;
; PROCEDURES CALLED:
;   combine1fiber
;   correct_dlam
;   divideflat
;   djs_diff_angle()
;   idlspec2d_version()
;   mkhdr
;   modfits
;   mrdfits()
;   pixelmask_bits()
;   platesn
;   splog
;   sxaddpar
;   sxdelpar
;   sxpar()
;   traceset2xy
;   writefits
;
; INTERNAL SUPPORT PROCEDURES:
;   makelabel()
;
; REVISION HISTORY:
;   02-Jan-2000  Written by D. Schlegel; modified from COMBINE2DOUT
;-
;------------------------------------------------------------------------------
function makelabel, hdr

   camera = strtrim(sxpar(hdr, 'CAMERAS'),2)
   mjd = strtrim(string(sxpar(hdr, 'MJD')),2)
   expos =  strtrim(string(sxpar(hdr, 'EXPOSURE')),2)
   flat  =  strmid(sxpar(hdr, 'FLATFILE'),6,9)
   arc   =  strmid(sxpar(hdr, 'ARCFILE'),6,9)

   label = string(camera, "-", mjd, "-", expos, flat, arc, $
            format='(a2,a1,i5.5,a1,i8.8,a9,a9)')

   return, label
end

;------------------------------------------------------------------------------
pro add_iraf_keywords, hdr, wavemin, binsz

   sxaddpar, hdr, 'WAT0_001', 'system=linear'
   sxaddpar, hdr, 'WAT1_001', $
    'wtype=linear label=Wavelength units=Angstroms'
   sxaddpar, hdr, 'CRVAL1', wavemin, $
    ' Central wavelength (log10) of first pixel'
   sxaddpar, hdr, 'CD1_1', binsz, ' Log10 dispersion per pixel'
   sxaddpar, hdr, 'CRPIX1', 1, ' Starting pixel (1-indexed)'
   sxaddpar, hdr, 'CTYPE1', 'LINEAR'
   sxaddpar, hdr, 'DC-FLAG', 1, ' Log-linear flag'

   return
end

;------------------------------------------------------------------------------
pro spcoadd_v5, spframes, outputname, $
 mjd=mjd, binsz=binsz, zeropoint=zeropoint, nord=nord, wavemin=wavemin, $
 bkptbin=bkptbin, window=window, maxsep=maxsep, adderr=adderr, $
 docams=camnames, plotsnfile=plotsnfile, combinedir=combinedir

   if (NOT keyword_set(binsz)) then binsz = 1.0d-4 $
    else binsz = double(binsz)
   if (NOT keyword_set(zeropoint)) then zeropoint = 3.5D
   if (n_elements(window) EQ 0) then window = 100
   if (NOT keyword_set(combinedir)) then combinedir=''

   ;----------
   ; Sort filenames such that this list contains first the blue then the red
; Why should this matter to sort them ???

   nfiles = n_elements(spframes)
   if (nfiles EQ 0) then return

   filenames = spframes[sort(spframes)]

   ;---------------------------------------------------------------------------

   if NOT keyword_set(camnames) then camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)
   exptimevec = fltarr(ncam)

   ;---------------------------------------------------------------------------
   ; Loop through each 2D output and read in the data
   ;---------------------------------------------------------------------------
   for ifile=0, nfiles-1 do begin

      ;----------
      ; Read in all data from this input file.
      ; Reading the plug-map structure will fail if its structure is
      ; different between different files.

      splog, 'Reading file #', ifile, ': ', filenames[ifile]
      spframe_read, filenames[ifile], objflux=tempflux, objivar=tempivar, $
       mask=temppixmask, wset=tempwset, dispset=tempdispset, plugmap=tempplug, $
       skyflux=tempsky, hdr=hdr, adderr=adderr

      if (ifile EQ 0) then $
       hdrarr = ptr_new(hdr) $
      else $
       hdrarr = [hdrarr, ptr_new(hdr)]

      ;----------
      ; Read header info

      thismjd = sxpar(hdr, 'MJD')
      if (NOT keyword_set(mjdlist)) then mjdlist = thismjd $
       else mjdlist = [mjdlist, thismjd]
      cameras = strtrim(sxpar(hdr, 'CAMERAS'),2)
      expstr = string(sxpar(hdr, 'EXPOSURE'), format='(i8.8)')

      ;----------
      ; Solve for wavelength and lambda-dispersion at each pixel in the image

      traceset2xy, tempwset, junk, tempwave
      traceset2xy, tempdispset, junk, tempdispersion

      ;----------
      ; Here is the correct conversion from pixels to log-lambda dispersion.
      ; We are converting from the dispersion in units of spFrame pixel sizes
      ; to the dispersion in units of the new rebinned pixel size, which is
      ; BINSZ in log-lambda units.

      correct_dlam, tempdispersion, 0, tempwset, dlam=binsz, /inverse

      ;----------

      dims = size(tempflux, /dimens)
      npix = dims[0]
      nfib = dims[1]

      ;----------
      ; Make a map of the size of each pixel in delta-(log10-Angstroms),
      ; and re-normalize the flux to ADU/(dloglam)

      correct_dlam, tempflux, tempivar, tempwset, dlam=binsz
      correct_dlam, tempsky, 0, tempwset, dlam=binsz

      ;----------
      ; Determine if this is a blue or red spectrum

      icam = (where(cameras EQ camnames))[0]
      if (icam EQ -1) then $
       message, 'Unknown camera ' + cameras
      exptimevec[icam] = exptimevec[icam] + sxpar(hdr, 'EXPTIME')

      ;----------
      ; Apply spectro-photometric calibration

      expnum = sxpar(hdr, 'EXPOSURE')
      calibfile = djs_filepath(string(camnames[icam], expnum, $
       format='("spFluxcalib-", a2, "-", i8.8, ".fits")'), $
       root_dir=combinedir)
      calibfile = (findfile(calibfile+'*'))[0]

      calibfac = mrdfits(calibfile, 0, calibhdr, /silent)
      minval = 0.05 * mean(calibfac)
      divideflat, tempflux, invvar=tempivar, calibfac, minval=minval
      divideflat, tempsky, calibfac, minval=minval
      temppixmask = temppixmask $
       OR (calibfac LE minval) * pixelmask_bits('BADFLUXFACTOR')

      ;----------
      ; Apply flux-correction factor between spectro-photometric exposure
      ; and this exposure.  There's also an optional additive term.
      ; So the flux is first multiplied by HDU#0, then we add HDU#1.

      corrfile = djs_filepath(string(camnames[icam], expnum, $
       format='("spFluxcorr-", a2, "-", i8.8, ".fits")'), $
       root_dir=combinedir)
      thisfile = (findfile(corrfile+'*'))[0]
      if (NOT keyword_set(thisfile)) then $
       message,' Could not find flux-corr file ' + corrfile

      aterm = mrdfits(thisfile, 0, corrhdr, /silent)
      bterm = mrdfits(thisfile, 1)
      invertcorr = 1. / aterm
      minval = 0.05 * invertcorr
      divideflat, tempflux, invvar=tempivar, invertcorr, minval=minval
      tempflux = tempflux + bterm
      divideflat, tempsky, invertcorr, minval=minval
      temppixmask = temppixmask $
       OR (invertcorr LE minval) * pixelmask_bits('BADFLUXFACTOR')

      ;----------
      ; Apodize the errors
      ; Do this only for the dichroic overlap region, which are the first
      ; rows in both the blue and red CCD's.

      if (keyword_set(window)) then begin
         swin = window < npix
         indx = lindgen(swin)
         tempivar[indx,*] = $
          tempivar[indx,*] * (indx # replicate(1,nfib)) / swin
      endif

      ;----------
      ; Concatenate data from all images

      if (ifile EQ 0) then begin
         flux = tempflux
         fluxivar = tempivar
         wave = tempwave
         dispersion = tempdispersion
         pixelmask = temppixmask
         skyflux = tempsky

         camerasvec = cameras
         label = makelabel(hdr)
         filenum = lonarr(nfib) + ifile
         plugmap = tempplug
      endif else begin
         ; Append as images...
         flux = [[flux], [tempflux]]
         fluxivar = [[fluxivar], [tempivar]]
         wave = [[wave], [tempwave]]
         dispersion = [[dispersion], [tempdispersion]]
         pixelmask = [[pixelmask], [temppixmask]]
         skyflux = [[skyflux], [tempsky]]

         ; Append as vectors...
         camerasvec = [camerasvec, cameras]
         label = [label, makelabel(hdr)]
         filenum = [filenum, lonarr(nfib) + ifile]
         plugmap = [plugmap, tempplug]
      endelse

   endfor

   tempflux = 0
   tempivar = 0
   tempwave = 0
   tempdispersion = 0
   temppixmask = 0
   tempsky = 0

   ;----------
   ; Check how many exposures we have in each of the (4) cameras

   for icam=0, ncam-1 do begin
      junk = where(camerasvec EQ camnames[icam], nmatch)
      splog, 'Files for camera ' + camnames[icam] + ':', nmatch
      if (icam EQ 0) then nminfile = nmatch $
       else nminfile = nminfile < nmatch
   endfor
; ??? Should make this routine robust to fewer files!!!
   if (nminfile LT 2) then begin
      splog, 'ABORT: At least 2 files needed for each camera'
      return
   endif

   ;---------------------------------------------------------------------------
   ; Construct output data structures, including the wavelength scale
   ;---------------------------------------------------------------------------

   totalpix = (size(flux, /dimens))[0]

   nonzero = where(fluxivar GT 0.0)
   minfullwave = min(wave[nonzero])
   maxfullwave = max(wave[nonzero])

   ; Get max and min wavelength from good pixels

   if (NOT keyword_set(wavemin)) then begin
      spotmin = long((minfullwave - zeropoint)/binsz) + 1L
      spotmax = long((maxfullwave - zeropoint)/binsz)
      wavemin = spotmin * binsz + zeropoint
      wavemax = spotmax * binsz + zeropoint
   endif else begin
      spotmin = 0L
      if (NOT keyword_set(wavemax)) then begin
        spotmax = long((maxfullwave - wavemin)/binsz)
        wavemax = spotmax * binsz + wavemin
      endif else spotmax = long((wavemax - wavemin)/binsz)
   endelse

   nfinalpix = spotmax - spotmin + 1L
   finalwave = dindgen(nfinalpix) * binsz + wavemin

   nfiber = max(plugmap.fiberid)

   finalflux = fltarr(nfinalpix, nfiber)
   finalivar = fltarr(nfinalpix, nfiber)
   finalandmask = lonarr(nfinalpix, nfiber)
   finalormask = lonarr(nfinalpix, nfiber)
   finaldispersion = fltarr(nfinalpix, nfiber)
   finalsky = fltarr(nfinalpix, nfiber)
   finalplugmap = replicate(plugmap[0], nfiber)
   struct_assign, {fiberid: 0L}, finalplugmap ; Zero out all elements in this
                                              ; FINALPLUGMAP structure.

   ;----------
   ; Issue a warning about any object fibers with OBJTYPE = 'NA', which
   ; should be impossible, but the special plate 673 and possibly others
   ; had some such fibers.

   ibad = where(strtrim(plugmap.objtype,2) EQ 'NA', nbad)
   if (nbad GT 0) then $
    splog, 'WARNING: ', nbad, ' fibers have OBJTYPE=NA in the plug-map'

   ;---------------------------------------------------------------------------
   ; Combine each fiber, one at a time
   ;---------------------------------------------------------------------------

   for ifiber=0, nfiber-1 do begin

      ; Find the first occurance of fiber number IFIBER+1
      indx = (where(plugmap.fiberid EQ ifiber+1))[0]

      if (indx NE -1) then begin
         splog, 'Fiber', ifiber+1, ' ', plugmap[indx].objtype, $
          plugmap[indx].mag, format = '(a, i4.3, a, a, f6.2, 5f6.2)'

         finalplugmap[ifiber] = plugmap[indx]

         ; Identify all objects within 2 arcsec of this position, and
         ; combine all these objects into a single spectrum.
         ; If all pluggings are identical, then this will always be
         ; the same fiber ID.
         ; Also, insist that the object type is not 'NA', which would
         ; occur for unplugged fibers.

         adist = djs_diff_angle(plugmap.ra, plugmap.dec, $
          plugmap[indx].ra, plugmap[indx].dec, units='degrees')
         indx = where(adist LT 2./3600. AND strtrim(plugmap.objtype,2) NE 'NA')
      endif

      if (indx[0] NE -1) then begin
         temppixmask = pixelmask[*,indx]
         combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
          finalmask=temppixmask, indisp=dispersion[*,indx], $
          skyflux=skyflux[*,indx], $
          newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
          andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
          newsky=bestsky, $
          nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
          maxiter=50, upper=3.0, lower=3.0, maxrej=1

         finalflux[*,ifiber] = bestflux
         finalivar[*,ifiber] = bestivar
         finalandmask[*,ifiber] = bestandmask
         finalormask[*,ifiber] = bestormask
         finaldispersion[*,ifiber] = bestdispersion
         finalsky[*,ifiber] = bestsky

         ; The following adds the COMBINEREJ bit to the input pixel masks
         pixelmask[*,indx] = temppixmask
      endif else begin
         splog, 'Fiber', ifiber+1, ' NO DATA'
         finalandmask[*,ifiber] = pixelmask_bits('NODATA')
         finalormask[*,ifiber] = pixelmask_bits('NODATA')
      endelse
   endfor

   ;---------------------------------------------------------------------------
   ; Write the corrected spCFrame files.
   ; All the fluxes + their errors are calibrated.
   ; The wavelengths + dispersions are converted from trace sets to images.
   ; The pixel mask has the COMBINEREJ bit set.
; The plug map should have the modified magnitudes !!!???
   ;---------------------------------------------------------------------------

   for ifile=0, nfiles-1 do begin
      thisfile = fileandpath(filenames[ifile], path=thispath)
      thisfile = djs_filepath(repstr(thisfile,'spFrame','spCFrame'), $
       root_dir=thispath)
      splog, 'Writing file #', ifile, ': ', thisfile
      indx = where(filenum EQ ifile)
; Change the units in the headers to be erg/sec !!!???
      mwrfits, flux[*,indx], thisfile, *hdrarr[ifile], /create
      mwrfits, fluxivar[*,indx], thisfile
      mwrfits, pixelmask[*,indx], thisfile
      mwrfits, wave[*,indx], thisfile
      mwrfits, dispersion[*,indx], thisfile
      mwrfits, plugmap[indx], thisfile
      mwrfits, skyflux[*,indx], thisfile
   endfor

   ;----------
   ; Clear memory

   wave = 0
   flux = 0
   fluxivar = 0
   temppixmask = 0
   dispersion = 0
   skyflux = 0

   ;---------------------------------------------------------------------------
   ; Generate S/N plots
   ;---------------------------------------------------------------------------

   ; Modify the 1st file's header to use for the combined plate header.

   hdr = *hdrarr[0]

   platesn, finalflux, finalivar, finalandmask, finalplugmap, finalwave, $
     hdr=hdr, plotfile=djs_filepath(plotsnfile, root_dir=combinedir)

   ;---------------------------------------------------------------------------
   ; Create the output header
   ;---------------------------------------------------------------------------

   ;----------
   ; Remove header cards that were specific to this first exposure
   ; (where we got the header).

   ncoeff = sxpar(hdr, 'NWORDER')
   for i=2, ncoeff-1 do sxdelpar, hdr, 'COEFF'+strtrim(string(i),2)

   sxdelpar, hdr, ['SPA', 'IPA', 'IPARATE']
   sxdelpar, hdr, 'EXPOSURE'
   sxdelpar, hdr, 'SEQID'
   sxdelpar, hdr, 'DARKTIME'
   sxdelpar, hdr, 'CAMERAS'
   sxdelpar, hdr, 'PLUGMAPO'
   for i=1, 4 do sxdelpar, hdr, 'GAIN'+strtrim(string(i),2)
   for i=1, 4 do sxdelpar, hdr, 'RDNOISE'+strtrim(string(i),2)
   sxdelpar, hdr, ['CAMCOL', 'CAMROW']
   sxdelpar, hdr, ['AMPLL', 'AMPLR', 'AMPUL', 'AMPUR']
   sxdelpar, hdr, ['FFS', 'FF', 'NE', 'HGCD']
   sxdelpar, hdr, ['SPEC1', 'SPEC2']
   sxdelpar, hdr, 'NBLEAD'
   sxdelpar, hdr, 'PIXFLAT'
   sxdelpar, hdr, 'PIXBIAS'
   sxdelpar, hdr, 'FLATFILE'
   sxdelpar, hdr, 'ARCFILE'
   sxdelpar, hdr, 'OBJFILE'
   sxdelpar, hdr, 'FRAMESN2'

   ;----------
   ; Average together some of the fields from the individual headers.

   cardname = [ 'AZ', 'ALT', 'TAI', 'WTIME', 'AIRTEMP', 'DEWPOINT', $
    'DEWDEP', 'DUSTA', 'DUSTB', 'DUSTC', 'DUSTD', 'GUSTS', 'HUMIDITY', $
    'HUMIDOUT', 'PRESSURE', 'WINDD', 'WINDS', 'TEMP01', 'TEMP02', $
    'TEMP03', 'TEMP04', 'HELIO_RV', 'SEEING20', 'SEEING50', 'SEEING80', $
    'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
    'WSIGMA', 'XSIGMA' ]
   sxcombinepar, hdrarr, cardname, hdr, func='average'

   sxcombinepar, hdrarr, 'TAI-BEG', hdr, func='min'
   sxcombinepar, hdrarr, 'TAI-END', hdr, func='max'

   sxcombinepar, hdrarr, 'XCHI2', hdr, func='max', outcard='XCHI2MAX', $
    after='XCHI2'
   sxcombinepar, hdrarr, 'XCHI2', hdr, func='min', outcard='XCHI2MIN', $
    after='XCHI2'

   sxcombinepar, hdrarr, 'SKYCHI2', hdr, func='max', outcard='SCHI2MAX', $
    after='SKYCHI2'
   sxcombinepar, hdrarr, 'SKYCHI2', hdr, func='min', outcard='SCHI2MIN', $
    after='SKYCHI2'

   sxcombinepar, hdrarr, 'WSIGMA', hdr, func='max', outcard='WSIGMAX', $
    after='WSIGMA'
   sxcombinepar, hdrarr, 'WSIGMA', hdr, func='min', outcard='WSIGMIN', $
    after='WSIGMA'

   sxcombinepar, hdrarr, 'XSIGMA', hdr, func='max', outcard='XSIGMAX', $
    after='XSIGMA'
   sxcombinepar, hdrarr, 'XSIGMA', hdr, func='min', outcard='XSIGMIN', $
    after='XSIGMA'

   ; Add the NGUIDE keywords for all headers of one flavor of CAMERAS
   ; (e.g., for all the 'b1' exposures if the first frame is 'b1'.)
   cardname = 'NGUIDE'
   sxcombinepar, hdrarr[0], cardname, hdr, func='total'
   cameras0 = sxpar(*(hdrarr[0]), 'CAMERAS')
   for ihdr=1, n_elements(hdrarr)-1 do begin
      if (sxpar(*(hdrarr[ihdr]), 'CAMERAS') EQ cameras0) then $
       sxcombinepar, hdrarr[ihdr], cardname, hdr, func='total'
   endfor

   ;----------
   ; Use the MJD passed as a keyword, which will typically be for the most
   ; observation, and be consistent with the output file names

   if (keyword_set(mjd)) then $
    sxaddpar, hdr, 'MJD', mjd

   ; Get the list of MJD's used for these reductions, then convert to a string
   mjdlist = mjdlist[uniq(mjdlist, sort(mjdlist))]
   mjdlist = strtrim(strcompress(string(mjdlist,format='(99a)')),2)
   sxaddpar, hdr, 'MJDLIST', mjdlist, after='MJD'

   ;----------
   ; Add new header cards

   sxaddpar, hdr, 'VERSCOMB', idlspec2d_version(), $
    ' Version of idlspec2d for combining multiple spectra', after='VERS2D'
   sxaddpar, hdr, 'NEXP', nfiles, $
    ' Number of exposures in this file', before='EXPTIME'
   for ifile=0,nfiles-1 do $
    sxaddpar, hdr, string('EXPID',ifile, format='(a5,i2.2)'), label[ifile], $
     ' ID string for exposure '+strtrim(string(ifile),2), before='EXPTIME'

   sxaddpar, hdr, 'EXPTIME', min(exptimevec), $
    ' Minimum of exposure times for all cameras'
   for icam=0, ncam-1 do $
    sxaddpar, hdr, 'EXPT_'+camnames[icam], exptimevec[icam], $
     ' '+camnames[icam]+' camera exposure time (seconds)', before='EXPTIME'
   sxaddpar, hdr, 'SPCOADD', systime(), $
    ' SPCOADD finished', after='EXPTIME'

   sxaddpar, hdr, 'NWORDER', 2, ' Linear-log10 coefficients'
   sxaddpar, hdr, 'NWORDER', 2, ' Linear-log10 coefficients'
   sxaddpar, hdr, 'WFITTYPE', 'LOG-LINEAR', ' Linear-log10 dispersion'
   sxaddpar, hdr, 'COEFF0', wavemin, $
    ' Central wavelength (log10) of first pixel'
   sxaddpar, hdr, 'COEFF1', binsz, ' Log10 dispersion per pixel'

   sxaddpar, hdr, 'NAXIS1', n_elements(bestflux)
   sxaddpar, hdr, 'NAXIS2', nfiber

   spawn, 'uname -n', uname
   sxaddpar, hdr, 'UNAME', uname[0]

   ;----------
   ; Check for smear exposure used and place info in header

   smearused = total((finalandmask AND pixelmask_bits('SMEARIMAGE')) NE 0) $
    GT 0 ? 'T' : 'F'
   sxaddpar, hdr, 'SMEARUSE', smearused, ' Smear image used?'

   ;----------
   ; Compute the fraction of bad pixels in total, and on each spectrograph.
   ; Bad pixels are any with SKYMASK(INVVAR)=0, excluding those where
   ; the NODATA bit is set in the pixel mask.

   ifib1 = where(finalplugmap.spectrographid EQ 1, nfib1)
   ifib2 = where(finalplugmap.spectrographid EQ 2, nfib2)
   qbadpix = skymask(finalivar, finalandmask, finalormask) EQ 0 $
    AND (finalandmask AND pixelmask_bits('NODATA')) EQ 0
   if (nfib1 GT 0) then $
    fbadpix1 = total(qbadpix[*,ifib1]) / (nfib1 * nfinalpix)
   if (nfib2 GT 0) then $
    fbadpix2 = total(qbadpix[*,ifib2]) / (nfib2 * nfinalpix)
   if (nfib1 GT 0 AND nfib2 GT 0) then $
    fbadpix = total(qbadpix[*,[ifib1,ifib2]]) / ((nfib1+nfib2) * nfinalpix) $
   else if (nfib1 GT 0) then $
    fbadpix = fbadpix1 $
   else if (nfib2 GT 0) then $
    fbadpix = fbadpix1 $
   else $
    fbadpix = 0

   sxaddpar, hdr, 'FBADPIX', fbadpix, ' Fraction of bad pixels'
   sxaddpar, hdr, 'FBADPIX1', fbadpix1, ' Fraction of bad pixels on spectro-1'
   sxaddpar, hdr, 'FBADPIX2', fbadpix2, ' Fraction of bad pixels on spectro-2'

   ;----------
   ; Add keywords for IRAF-compatability

   add_iraf_keywords, hdr, wavemin, binsz

   mkhdr, hdrfloat, finalivar, /image, /extend
   add_iraf_keywords, hdrfloat, wavemin, binsz

   mkhdr, hdrlong, finalandmask, /image, /extend
   add_iraf_keywords, hdrlong, wavemin, binsz

   ;---------------------------------------------------------------------------
   ; Write combined output file
   ;---------------------------------------------------------------------------

   fulloutname = djs_filepath(outputname, root_dir=combinedir)

   ; 1st HDU is flux
   sxaddpar, hdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
   mwrfits, finalflux, fulloutname, hdr, /create

   ; 2nd HDU is inverse variance
   sxaddpar, hdrfloat, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
   mwrfits, finalivar, fulloutname, hdrfloat

   ; 3rd HDU is AND-pixelmask
   mwrfits, finalandmask, fulloutname, hdrlong

   ; 4th HDU is OR-pixelmask
   mwrfits, finalormask, fulloutname, hdrlong

   ; 5th HDU is dispersion map
   sxaddpar, hdrfloat, 'BUNIT', 'pixels'
   mwrfits, finaldispersion, fulloutname, hdrfloat

   ; 6th HDU is plugmap
   mwrfits, finalplugmap, fulloutname

   ; 7th HDU is the sky
   mwrfits, finalsky, fulloutname

   return
end
;------------------------------------------------------------------------------
