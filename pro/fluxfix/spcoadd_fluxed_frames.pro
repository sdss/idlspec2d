;+
; NAME:
;   spcoadd_frames
;
; PURPOSE:
;   Combine several reduced frames of the same objects
;
; CALLING SEQUENCE:
;   spcoadd_frames, spframes, outputname, $
;    fcalibprefix=, [ mjd=, binsz=, zeropoint=, nord=, wavemin=, $
;    bkptbin=, window=, maxsep=, adderr=, plotsnfile=, combinedir= ]
;
; INPUTS:
;   spframes       - Name(s) of files to combine (written by SPREDUCE)
;   outputname     - Output file name
;
; REQUIRED KEYWORDS:
;   fcalibprefix   - Prefix for flux-calibration files.
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
pro spcoadd_fluxed_frames, spframes, outputname, fcalibprefix=fcalibprefix, $
 mjd=mjd, binsz=binsz, zeropoint=zeropoint, nord=nord, wavemin=wavemin, $
 bkptbin=bkptbin, window=window, maxsep=maxsep, adderr=adderr, $
 docams=camnames, plotsnfile=plotsnfile, combinedir=combinedir, $
 tsobjname = tsobjname, smearname = smearname

   ;---------------------------------------------------------------------------

   if (NOT keyword_set(binsz)) then binsz = 1.0d-4 $
    else binsz = double(binsz)
   if (NOT keyword_set(zeropoint)) then zeropoint = 3.5D
   if (n_elements(window) EQ 0) then window = 100
   if (NOT keyword_set(combinedir)) then combinedir=''

   ;----------
   ; Sort filenames such that this list contains first the blue then the red

   nfiles = n_elements(spframes)
   if (nfiles EQ 0) then return

   filenames = spframes[sort(spframes)]

   if NOT keyword_set(camnames) then camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   ;------------------
   ; Get plate + mjd strings
   
   words = strsplit(fcalibprefix, '-', /extract)
   prefix = words[0]
   words = strsplit(outputname, '-', /extract)
   plate_str = words[1]
   mjd_str = words[2]

   ;---------------------------------------------------------------------------
   ; Loop through each 2D output and read in the data
   ;---------------------------------------------------------------------------
 
   spread_frames, filenames, window=window, binsz = binsz, $
     adderr=adderr, camnames=camnames, tsobjname = tsobjname, $
     flux = flux, ivar = fluxivar, wave = wave, dispersion = dispersion, $
     pixelmask = pixelmask, plugmap = plugmap, fibertag = fibertag, $
     camerasvec = camerasvec, label = label, filenum = filenum,  $
     expid = expid, sn2 = sn2, exptimevec = exptimevec, mjdlist = mjdlist, $
     hdrarr = hdrarr

   fibertag.plateid = plate_str
   fibertag.mjd = mjd_str

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

   ;-----------------------------------------------------------------------
   ; Do spectral typing of standards & derive sphoto correction for 
   ; each frame, then apply it
   ;-----------------------------------------------------------------------

   sphoto_err = fltarr(n_elements(camerasvec))

   ; Loop through cameras
   for icam=0, ncam-1 do begin
     camid = camnames[icam]
     camcol = strmid(camnames[icam], 0,1)
     specnum = strmid(camnames[icam], 1,1)
     frames = expid[where(camerasvec eq camid)]

     ; Create output flux calibration file names
     fcalfiles = djs_filepath(prefix + '-' + camid + '-' + frames + '.fits', $
                 root_dir=combinedir)
     stdstarfile = djs_filepath('spStd-' + plate_str + '-' + mjd_str +  $
                   '-' + specnum + '.fits', root_dir=combinedir)

     ;--------------------------------
     ; Use plugmap to find standard stars

     isphoto = where((strtrim(fibertag.objtype) EQ 'SPECTROPHOTO_STD' OR $
               strtrim(fibertag.objtype) EQ 'REDDEN_STD') AND $
               (fibertag.spectrographid eq specnum) AND $
               (fibertag.camcolor eq camcol))

     ;---------------------------------
     ; Compute spectral types and write flux calibration vectors

     sphoto_calib, wave[*,isphoto], flux[*,isphoto], fluxivar[*,isphoto], $
                 pixelmask[*,isphoto], fibertag[isphoto], $
                 fcalfiles, stdstarfile, stype = (camcol eq 'b'), $
                 pca_calibset = pca_calibset

     ;---------------------------------
     ; Apply sphoto calibration to all fibers in each frame

     for iframe = 0, n_elements(frames) - 1 do begin

       indx = where(fibertag.expid eq frames[iframe] AND $
                    fibertag.camcolor eq camcol AND $ 
                    fibertag.spectrographid eq specnum)

       junk = mrdfits(fcalfiles[iframe], 0, calibhdr, /silent)
       calibset = mrdfits(fcalfiles[iframe], 1)

       ; Store measures of spectrophotometry error to use later in 
       ; deciding which exposure is the best
       sphoto_err[where(camerasvec eq camid and expid eq frames[iframe])] $
                  = sxpar(calibhdr, 'SPHOTERR')

       cwavemin = sxpar(calibhdr, 'WAVEMIN')
       cwavemax = sxpar(calibhdr, 'WAVEMAX')
       calibfac = bspline_valu(wave[*,indx], calibset)
         
       ; Set to bad any pixels whose wavelength is outside the known
       ; flux-calibration region.
       ibad = where(wave[*,indx] LT alog10(cwavemin) OR $
                    wave[*,indx] GT alog10(cwavemax))
       if (ibad[0] NE -1) then calibfac[ibad] = 0

       tempflux = flux[*,indx]
       tempivar = fluxivar[*,indx]
       divideflat, tempflux, invvar=tempivar, calibfac, $
                   minval=0.005*mean(calibfac)

       flux[*,indx] = tempflux
       fluxivar[*,indx] = tempivar
       pixelmask[*,indx] = pixelmask[*,indx] $
         OR (calibfac LE 0.005*mean(calibfac)) * pixelmask_bits('BADFLUXFACTOR')
     endfor
   endfor

   ;-----------------------------------------------------------------------
   ; Correct each spectrum for exposure to exposure differences (formerly 
   ; using the smear exposure -- now using the best quality exposure) 
   ; Note -- this is done in blue-red pairs!
   ;-----------------------------------------------------------------------

   iframe_b1 = where(fibertag.camcolor eq 'b' and fibertag.spectrographid eq 1)
   iframe_r1 = where(fibertag.camcolor eq 'r' and fibertag.spectrographid eq 1)
   iframe_b2 = where(fibertag.camcolor eq 'b' and fibertag.spectrographid eq 2)
   iframe_r2 = where(fibertag.camcolor eq 'r' and fibertag.spectrographid eq 2)

   ;---------------
   ; Determine which exposure is best from the S/N and spectrophotometry errors
  
   ; add up the S/N^2 of the 4 frames per exposure (b1, b2, r1, r2)
   sn2_exp = sn2[where(camerasvec eq 'b1')] + sn2[where(camerasvec eq 'r1')] $
           + sn2[where(camerasvec eq 'b2')] + sn2[where(camerasvec eq 'r2')] 

   ; add up the spectrophotometry errors for the 4 frames per exposure
   sphoto_sn2 = (1 / sphoto_err)^2
   sphoto_exp = sphoto_err[where(camerasvec eq 'b1')] + $
                sphoto_err[where(camerasvec eq 'r1')] + $
                sphoto_err[where(camerasvec eq 'b2')] + $
                sphoto_err[where(camerasvec eq 'r2')] 

   ; Pick the best exposure (want the same one for spec 1 & 2)
   maxval = max(sn2_exp/median(sn2_exp) + $
                sphoto_exp/median(sphoto_exp), iframe_best)
   uniqexp = expid[where(camerasvec eq 'b1')]
   best_exp = uniqexp[iframe_best]
   splog, 'Best Exposure is: ' + best_exp

   ;------------------
   ; Compute the exposure-to-exposure corrections (1/fiber)

   corrfiles1 = 'spFluxcorr-' + uniqexp + '-1.fits'
   corrfiles2 = 'spFluxcorr-' + uniqexp + '-2.fits'

   frame_flux_tweak, wave[*, iframe_b1], wave[*, iframe_r1], $
                   flux[*, iframe_b1], flux[*, iframe_r1], $
                   fluxivar[*, iframe_b1], fluxivar[*, iframe_r1], $
                   best_exp, fibertag[iframe_b1], corrfiles1
  
   frame_flux_tweak, wave[*, iframe_b2], wave[*, iframe_r2], $
                   flux[*, iframe_b2], flux[*, iframe_r2], $
                   fluxivar[*, iframe_b2], fluxivar[*, iframe_r2], $
                   best_exp, fibertag[iframe_b2], corrfiles2
  
   ; save the magnitude of this correction? -- goodness of spec photo?
   ; COADD_QUALITY ==> sigma of poly1 coeff

   ;--------------------------------------
   ; Read back the corrections and apply

   ; Loop through each exposure & camera
   for iexp = 0, n_elements(expid) - 1 do begin
     corrfile = 'spFluxcorr-' + expid[iexp] + '-' + $
                strmid(camerasvec[iexp], 1, 1) + '.fits'
     corrset = mrdfits(corrfile, 1)

     indx = where(fibertag.expid eq expid[iexp] AND $
                  fibertag.camcolor eq strmid(camerasvec[iexp], 0, 1) AND $ 
                  fibertag.spectrographid eq strmid(camerasvec[iexp], 1, 1))

     traceset2xy, corrset, wave[*,indx], corrimg

     ; Don't let the flux correction be more than a factor of 5!!
     invertcorr = 1.0 / corrimg
     tempflux = flux[*,indx]
     tempivar = fluxivar[*,indx]
     divideflat, tempflux, invvar=tempivar, invertcorr, minval=0.25

     flux[*,indx] = tempflux
     fluxivar[*,indx] = tempivar
     pixelmask[*,indx] = pixelmask[*,indx] OR $
                         (corrimg GE 4) * pixelmask_bits('BADFLUXFACTOR')
     pixelmask[*,indx] = pixelmask[*,indx] OR $
                         (corrimg LE 0.25) * pixelmask_bits('BADFLUXFACTOR')

   endfor

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
   finalplugmap = replicate(plugmap[0], nfiber)
   struct_assign, {fiberid: 0L}, finalplugmap ; Zero out all elements in this
                                              ; FINALPLUGMAP structure.

   finalfibertag = replicate(fibertag[0], nfiber)
   struct_assign, {fiberid: 0L}, finalfibertag

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
         finalfibertag[ifiber] = fibertag[indx]

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
           newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
           andmask = bestandmask, ormask=bestormask, newdisp=bestdispersion, $
           nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
           maxiter=50, upper=3.0, lower=3.0, maxrej=1

         finalflux[*,ifiber] = bestflux
         finalivar[*,ifiber] = bestivar
         finalandmask[*,ifiber] = bestandmask
         finalormask[*,ifiber] = bestormask
         finaldispersion[*,ifiber] = bestdispersion

         ; The following adds the COMBINEREJ bit to the input pixel masks
         pixelmask[*,indx] = temppixmask

         ;-------------
         ; QA plot

         ;set_plot, 'x'

         ;plot, 10.0^finalwave, bestflux, xr=[3800, 9200], /xs
         ;for ii = 0, n_elements(indx) - 1 do $
         ;  djs_oplot, 10.0^wave[*,indx[ii]], djs_median(flux[*,indx[ii]], $
         ;             width=75, boundary='reflect'), color='red'

      endif else begin
         splog, 'Fiber', ifiber+1, ' NO DATA'
         finalandmask[*,ifiber] = pixelmask_bits('NODATA')
         finalormask[*,ifiber] = pixelmask_bits('NODATA')
      endelse
   endfor

   ;-------------
   ; Clear memory

   wave = 0
   flux = 0
   fluxivar = 0
   temppixmask = 0
   dispersion = 0

   ;--------------------------------------------------------------------------
   ; Measure the difference between the combined science and smear spectra
   ;--------------------------------------------------------------------------

   if keyword_set(smearname) then begin

     smearset = smear_compare(smearname, finalwave, finalflux, finalivar, $
       best_exp, plate_str, mjd_str, camnames = camnames, adderr=adderr, $
       combinedir = combinedir, tsobjname = tsobjname)

     smear_coeff = smearset.coeff
   endif else begin
     splog, 'No smear exposures found!'
     smear_coeff = dblarr(3, nfiber)
   endelse

   ;---------------------------------------------------------------------------
   ; Generate S/N plots
   ;---------------------------------------------------------------------------

   ; Modify the 1st file's header to use for the combined plate header.

   hdr = *hdrarr[0]
  
   ; Use fibertag instead of the plugmap b/c this contains updated 
   ; fibermags from the tsObj (if available)

   platesn, finalflux, finalivar, finalandmask, finalfibertag, finalwave, $
     hdr=hdr, plotfile=djs_filepath(plotsnfile, root_dir=combinedir), $
     snvec=snvec, synthmag=synthmag

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

   ifib1 = where(finalfibertag.spectrographid EQ 1, nfib1)
   ifib2 = where(finalfibertag.spectrographid EQ 2, nfib2)
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

   ; 7th HDU are S/N vectors for g,r,i
   mwrfits, snvec, fulloutname

   ; 8th HDU are synthetic magnitude vectors
   mwrfits, synthmag, fulloutname

   ; 9th HDU are the legendre coefficients of the smear correction vectors
   ; (3 x 640 array)
   mwrfits, smear_coeff, fulloutname

   ;---------------------------------------------------------------------------
   ; Write the modified pixel masks to the input files
   ;---------------------------------------------------------------------------

;   for ifile=0, nfiles-1 do begin
;      splog, 'Modifying file #', ifile, ': ', filenames[ifile]
;      indx = where(filenum EQ ifile)
;      djs_modfits, filenames[ifile], pixelmask[*,indx], exten_no=2
;   endfor

   return
end
;------------------------------------------------------------------------------
