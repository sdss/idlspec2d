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
;    bkptbin=, window=, maxsep=, adderr=, plotsnfile=, $
;    combinedir=, bestexpnum= ]
;
; INPUTS:
;   spframes       - Name(s) of spFrame files (written by SPREDUCE)
;   outputname     - Output file name
;   plotsnfile     - Name of output plot file
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
;   bkptbin        - Parameter for COMBINE1FIBER
;   window         - Window size for apodizing the errors of the spectrum
;                    from each individual frame;
;                    default to 100 pixels apodization on each end of the
;                    spectra.
;   maxsep         - Parameter for COMBINE1FIBER
;   adderr         - Additional error to add to the formal errors, as a
;                    fraction of the flux.
;   combinedir     - Optional output directory
;   bestexpnum     - Exposure number for best exposure, to which all other
;                    exposures are tied; this is only used in this procedure
;                    for logging in the FITS header
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   All input files must have the same number of pixels per spectrum,
;   i.e. 2048 wavelength samplings, although those wavelengths can
;   be different.
;
;   Flux-correction files are also read in, where they are assumed to
;   have the name spFluxcorr-EEEEEEEE-S.fits, where EEEEEEEE is the exposure
;   number and S is the spectrograph ID (1 or 2).
;
; EXAMPLES:
;
; BUGS:
;   This routine used to combine data from multiple (different) plug maps.
;   Objects are matched based upon their positions agreeing to 2 arc sec.
;   This is *not* true any longer, especially when applying the
;   flux-distortion solutions.
;
; PROCEDURES CALLED:
;   combine1fiber
;   correct_dlam
;   divideflat
;   djs_diff_angle()
;   fcalib_default()
;   fiber_rollcall
;   flux_distortion()
;   idlspec2d_version()
;   mkhdr
;   mrdfits()
;   mwrfits
;   pixelmask_bits()
;   platesn
;   splog
;   spframe_read
;   sxaddpar
;   sxdelpar
;   sxpar()
;   traceset2xy
;
; INTERNAL SUPPORT PROCEDURES:
;   makelabel()
;   add_iraf_keywords
;
; REVISION HISTORY:
;   02-Jan-2000  Written by D. Schlegel; modified from COMBINE2DOUT
;-
;------------------------------------------------------------------------------
function makelabel, hdr

   camera = strtrim(sxpar(hdr, 'CAMERAS'),2)
   expos =  strtrim(string(sxpar(hdr, 'EXPOSURE')),2)
   flat  =  strmid(sxpar(hdr, 'FLATFILE'),7,8)
   arc   =  strmid(sxpar(hdr, 'ARCFILE'),7,8)

   label = string(camera, expos, flat, arc, $
    format='(a2,"-",i8.8,"-",a8,"-",a8)')

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
 docams=camnames, plotsnfile=plotsnfile, combinedir=combinedir, $
 bestexpnum=bestexpnum

   if (NOT keyword_set(binsz)) then binsz = 1.0d-4 $
    else binsz = double(binsz)
   if (NOT keyword_set(zeropoint)) then zeropoint = 3.5D
   if (n_elements(window) EQ 0) then window = 100
   if (NOT keyword_set(combinedir)) then combinedir=''

   string1 = repstr(outputname,'spPlate','spFluxdistort')
   string2 = repstr(string1,'.fits','.ps')
   distortfitsfile = djs_filepath(string1, root_dir=combinedir)
   distortpsfile = djs_filepath(string2, root_dir=combinedir)

   ;----------
   ; Sort filenames such that this list contains first the blue then the red

   nfiles = n_elements(spframes)
   if (nfiles EQ 0) then return

   filenames = spframes[sort(spframes)]

   ;---------------------------------------------------------------------------

   if NOT keyword_set(camnames) then camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)
   nexpvec = lonarr(ncam)
   exptimevec = fltarr(ncam)

   ;---------------------------------------------------------------------------
   ; Loop through each 2D output and read in the data
   ;---------------------------------------------------------------------------

   ; Start by determining the size of all the files
   npixarr = lonarr(nfiles)
   for ifile=0, nfiles-1 do begin
      spframe_read, filenames[ifile], hdr=objhdr
      npixarr[ifile] = sxpar(objhdr,'NAXIS1')
   endfor
   npixmax = max(npixarr)
   nobj = sxpar(objhdr,'NAXIS2') ; Number of fibers per spectrograph

   for ifile=0, nfiles-1 do begin
      ;----------
      ; Read in all data from this input file.
      ; Reading the plug-map structure will fail if its structure is
      ; different between different files.

      splog, 'Reading file #', ifile, ': ', filenames[ifile], $
       prename=filenames[ifile]
      spframe_read, filenames[ifile], objflux=tempflux, objivar=tempivar, $
       mask=temppixmask, wset=tempwset, dispset=tempdispset, plugmap=tempplug, $
       skyflux=tempsky, ximg=tempximg, superflat=tempsuperflat, $
       hdr=hdr, adderr=adderr

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

      ;need to reset junk since the array lengths change
      junk=0
      traceset2xy, tempwset, junk, tempwave
      traceset2xy, tempdispset, junk, tempdispersion
      ;----------
      ; Here is the correct conversion from pixels to log-lambda dispersion.
      ; We are converting from the dispersion in units of spFrame pixel sizes
      ; to the dispersion in units of the new rebinned pixel size, which is
      ; BINSZ in log-lambda units.
      
      ; this probably should be fixed elsewhere but limit where the fit range
      tempxmax=tempwset.xmax
      tempwset.xmax=(size(tempdispersion,/dimens))[0]-1
      correct_dlam, tempdispersion, 0, tempwset, dlam=binsz, /inverse
      tempwset.xmax=tempxmax
   
      ;----------

      dims = size(tempflux, /dimens)
      npix = dims[0]
      nfib = dims[1]

      ;----------
      ; Make a map of the size of each pixel in delta-(log10-Angstroms),
      ; and re-normalize the flux to electrons/(dloglam)

      correct_dlam, tempflux, tempivar, tempwset, dlam=binsz
      correct_dlam, tempsky, 0, tempwset, dlam=binsz

      ;----------
      ; Determine if this is a blue or red spectrum

      icam = (where(cameras EQ camnames))[0]
      if (icam EQ -1) then $
       message, 'Unknown camera ' + cameras
      nexpvec[icam] = nexpvec[icam] + 1
      exptimevec[icam] = exptimevec[icam] + sxpar(hdr, 'EXPTIME')

      ;----------
      ; Apply spectro-photometric calibration

      expnum = sxpar(hdr, 'EXPOSURE')
      calibfile = djs_filepath(string(camnames[icam], expnum, $
       format='("spFluxcalib-", a2, "-", i8.8, ".fits")'), $
       root_dir=combinedir)
      calibfile = (findfile(calibfile+'*'))[0]

      if (keyword_set(calibfile)) then begin
         calibfac = mrdfits(calibfile, 0, calibhdr, /silent)
      endif else begin
         splog, 'WARNING: Reading default flux-calib vectors for camera=' $
          + camnames[icam]
         calibfac = fcalib_default(camnames[icam], tempwave, exptimevec[icam])
      endelse
      minval = 0.05 * mean(calibfac)
      divideflat, tempflux, invvar=tempivar, calibfac, minval=minval
      divideflat, tempsky, calibfac, minval=minval
      temppixmask = temppixmask $
       OR ((calibfac LE minval OR keyword_set(calibfile) EQ 0) $
       * pixelmask_bits('BADFLUXFACTOR'))

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
      minval = 0.05 / mean(aterm)
      nrownative=(size(tempflux,/dimens))[0]
      divideflat, tempflux, invvar=tempivar, invertcorr[0:nrownative-1,*], minval=minval
      tempflux = tempflux + bterm
      divideflat, tempsky, invertcorr[0:nrownative-1,*], minval=minval
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
         ; Construct the image arrays
         flux = make_array(npixmax,nobj*nfiles,type=size(tempflux,/type))
         fluxivar = make_array(npixmax,nobj*nfiles,type=size(tempivar,/type))
         wave = make_array(npixmax,nobj*nfiles,type=size(tempwave,/type))
         dispersion = make_array(npixmax,nobj*nfiles,type=size(tempdisp,/type))
         pixelmask = make_array(npixmax,nobj*nfiles,type=size(temppixmask,/type))
         skyflux = make_array(npixmax,nobj*nfiles,type=size(tempsky,/type))
         ximg = make_array(npixmax,nobj*nfiles,type=size(tempximg,/type))
         superflat = make_array(npixmax,nobj*nfiles,type=size(tempsuperflat,/type))

         ; Append as vectors...
         camerasvec = cameras
         label = makelabel(hdr)
         filenum = lonarr(nfib) + ifile
         plugmap = tempplug
      endif else begin
         ; Append as vectors...
         camerasvec = [camerasvec, cameras]
         label = [label, makelabel(hdr)]
         filenum = [filenum, lonarr(nfib) + ifile]
         plugmap = [plugmap, tempplug]
      endelse

      flux[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempflux
      fluxivar[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempivar
      wave[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempwave
      ; Pad the wavelengths with reasonable values
      if (npixarr[ifile] LT npixmax) then begin
         dwave = tempwave[npixarr[ifile]-1,*] - tempwave[npixarr[ifile]-2,*]
         addwave = tempwave[npixarr[ifile]-1,*] $
          ## (1+lonarr(npixmax-npixarr[ifile])) $
          + dwave ## (1+lindgen(npixmax-npixarr[ifile]))
         wave[npixarr[ifile]:npixmax-1,nobj*ifile:nobj*(ifile+1)-1] = addwave
      endif
      dispersion[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempdispersion
      pixelmask[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = temppixmask
      skyflux[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempsky
      ximg[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempximg
      superflat[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempsuperflat

      splog, prename=''
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
   if (nminfile LT 1) then begin
      splog, 'ABORT: At least 1 file needed for each camera'
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

;   nfiber = max(plugmap.fiberid)
   nfiber = 2 * nfib

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
          plugmap[indx].mag, format = '(a, i5.4, a, a, f6.2, 5f6.2)'

         finalplugmap[ifiber] = plugmap[indx]

         ; Identify all objects with the same XFOCAL,YFOCAL plate position, and
         ; combine all these objects into a single spectrum.
         ; If all pluggings are identical, then this will always be
         ; the same fiber ID.
         ; Also, insist that the object type is not 'NA', which would
         ; occur for unplugged fibers. <--- Disable this for BOSS ???

         indx = where(abs(plugmap.xfocal - plugmap[indx].xfocal) LT 0.0001 $
          AND abs(plugmap.yfocal - plugmap[indx].yfocal) LT 0.0001)
;          AND strtrim(plugmap.objtype,2) NE 'NA')
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
   ;----------
   ; Modify the 1st file's header to use for the combined plate header.

   bighdr = *hdrarr[0]

   ;---------------------------------------------------------------------------
   ; FLUX DISTORTION IMAGE
   ;---------------------------------------------------------------------------

   ; Compute the flux distortion image
   corrimg = flux_distortion(finalflux, finalivar, finalandmask, finalormask, $
    plugmap=finalplugmap, loglam=finalwave, plotfile=distortpsfile, hdr=bighdr)

   igood = where(finalivar GT 0)
   thismin = min(corrimg[igood], max=thismax)
   cratio = thismin / thismax
   if (cratio LT 1./100) then begin
      splog, 'WARNING: Flux distortion image dynamic range = ', 1./cratio, $
       ' (DISABLE)'
      corrimg[*] = 1.
   endif else begin
      splog, 'Flux distortion image dynamic range = ', 1./cratio
   endelse

   ; Plot S/N and throughput **before** this distortion-correction.
   splog, prelog='Initial'
   platesn, finalflux, finalivar, finalandmask, finalplugmap, finalwave, $
    hdr=bighdr, plotfile=djs_filepath(plotsnfile+'.orig', root_dir=combinedir)
   splog, prelog=''

   ; Apply this flux-distortion to the final, co-added fluxes.
   invcorrimg = 1. / corrimg
   minicorrval = 0.05 / mean(corrimg)
   divideflat, finalflux, invvar=finalivar, invcorrimg, minval=minicorrval
   divideflat, bestsky, invcorrimg, minval=minicorrval
   finalandmask = finalandmask $
    OR (invcorrimg LE minicorrval) * pixelmask_bits('BADFLUXFACTOR')
   finalormask = finalormask $
    OR (invcorrimg LE minicorrval) * pixelmask_bits('BADFLUXFACTOR')

   ; Plot S/N and throughput **after** this distortion-correction.
   ; (This over-writes header cards written in the first call.)
   splog, prelog='Final'
   platesn, finalflux, finalivar, finalandmask, finalplugmap, finalwave, $
    hdr=bighdr, plotfile=djs_filepath(plotsnfile, root_dir=combinedir)
   splog, prelog=''

   ;---------------------------------------------------------------------------
   ; Write the corrected spCFrame files.
   ; All the fluxes + their errors are calibrated.
   ; The wavelengths + dispersions are converted from trace sets to 2D images.
   ; The pixel mask has the COMBINEREJ bit set.
   ;---------------------------------------------------------------------------

   for ifile=0, nfiles-1 do begin
      thisfile = fileandpath(filenames[ifile], path=thispath)
      thisfile = djs_filepath(repstr(thisfile,'spFrame','spCFrame'), $
       root_dir=thispath)
      splog, 'Writing file #', ifile, ': ', thisfile, prename=filenames[ifile]
      indx = where(filenum EQ ifile, nthis)

      hdr = *hdrarr[ifile]
      sxaddpar, hdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'

      ; Apply the flux-distortion image to each individual frame, by
      ; interpolating off the full wavelength-scale distortion image
      ; onto the wavelength mapping of each individual exposure+CCD.
      for i=0L, nthis-1 do begin
         thisflux1 = flux[*,indx[i]]
         thisivar1 = fluxivar[*,indx[i]]
         thissky1 = skyflux[*,indx[i]]
         j = plugmap[indx[i]].fiberid - 1
         thisicorr = interpol(invcorrimg[*,j], finalwave, wave[*,indx[i]])
         divideflat, thisflux1, invvar=thisivar1, thisicorr, minval=minicorrval
         flux[*,indx[i]] = thisflux1
         fluxivar[*,indx[i]] = thisivar1

         divideflat, thissky1, thisicorr, minval=minicorrval
         skyflux[*,indx[i]] = thissky1
      endfor

      mwrfits, flux[*,indx], thisfile, hdr, /create
      
      fxaddpar, exthdr, 'EXTNAME', 'IVAR', ' Inverse variance'
      mwrfits, fluxivar[*,indx], thisfile, exthdr
      
      fxaddpar, exthdr, 'EXTNAME', 'MASK', ' Pixel mask'
      mwrfits, pixelmask[*,indx], thisfile, exthdr
      
      fxaddpar, exthdr, 'EXTNAME', 'WAVELENGTH', ' Wavelength solution'
      mwrfits, wave[*,indx], thisfile, exthdr
      
      fxaddpar, exthdr, 'EXTNAME', 'WAVEDISP', ' Wavelength dispersion'
      mwrfits, dispersion[*,indx], thisfile, exthdr
      
      ;; need a different header for plugmap structure
      ;; undefine it first so that mwrfits doesn't duplicate comments
      ;; on successive writes
      undefine, plughdr
      fxaddpar, plughdr, 'EXTNAME', 'PLUGMAP', ' Plugmap structure'
      mwrfits, plugmap[indx], thisfile, plughdr
      
      fxaddpar, exthdr, 'EXTNAME', 'SKY', ' Subtracted sky flux'
      mwrfits, skyflux[*,indx], thisfile, exthdr
      
      fxaddpar, exthdr, 'EXTNAME', 'X', ' Trace X locations on CCD'
      mwrfits, ximg[*,indx], thisfile, exthdr
      
      fxaddpar, exthdr, 'EXTNAME', 'SUPERFLAT', ' Superflat'
      mwrfits, superflat[*,indx], thisfile, exthdr
      
   endfor
   splog, prename=''

   ;----------
   ; Clear memory

   wave = 0
   flux = 0
   fluxivar = 0
   temppixmask = 0
   dispersion = 0
   skyflux = 0
   superflat = 0

   ;---------------------------------------------------------------------------
   ; Create the output header
   ;---------------------------------------------------------------------------

   ;----------
   ; Print roll call of bad fibers and bad pixels.

   fiber_rollcall, finalandmask, finalwave

   ;----------
   ; Remove header cards that were specific to this first exposure
   ; (where we got the header).

   ncoeff = sxpar(bighdr, 'NWORDER')
   for i=2, ncoeff-1 do sxdelpar, bighdr, 'COEFF'+strtrim(string(i),2)

   sxdelpar, bighdr, ['SPA', 'IPA', 'IPARATE']
   sxdelpar, bighdr, 'EXPOSURE'
   sxdelpar, bighdr, 'REQTIME'
   sxdelpar, bighdr, 'QUALITY'
   sxdelpar, bighdr, 'FILENAME'
   sxdelpar, bighdr, 'SEQID'
   sxdelpar, bighdr, 'DARKTIME'
   sxdelpar, bighdr, 'CAMERAS'
   sxdelpar, bighdr, 'PLUGMAPO'
   for i=1, 4 do sxdelpar, bighdr, 'GAIN'+strtrim(string(i),2)
   for i=1, 4 do sxdelpar, bighdr, 'RDNOISE'+strtrim(string(i),2)
   sxdelpar, bighdr, ['CAMCOL', 'CAMROW']
   sxdelpar, bighdr, ['AMPLL', 'AMPLR', 'AMPUL', 'AMPUR']
   sxdelpar, bighdr, ['FFS', 'FF', 'NE', 'HGCD']
   sxdelpar, bighdr, ['SPEC1', 'SPEC2']
   sxdelpar, bighdr, 'NBLEAD'
   sxdelpar, bighdr, 'PIXFLAT'
   sxdelpar, bighdr, 'PIXBIAS'
   sxdelpar, bighdr, 'FLATFILE'
   sxdelpar, bighdr, 'ARCFILE'
   sxdelpar, bighdr, 'OBJFILE'
   sxdelpar, bighdr, 'FRAMESN2'
   sxdelpar, bighdr, 'DEREDSN2'

   ;----------
   ; Average together some of the fields from the individual headers.

   cardname = [ 'AZ', 'ALT', 'TAI', 'WTIME', 'AIRTEMP', 'DEWPOINT', $
    'DEWDEP', 'DUSTA', 'DUSTB', 'DUSTC', 'DUSTD', 'GUSTS', 'HUMIDITY', $
    'HUMIDOUT', 'PRESSURE', 'WINDD', 'WINDS', 'TEMP01', 'TEMP02', $
    'TEMP03', 'TEMP04', 'HELIO_RV', 'SEEING20', 'SEEING50', 'SEEING80', $
    'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
    'WSIGMA', 'XSIGMA' ]
   sxcombinepar, hdrarr, cardname, bighdr, func='average'

   sxcombinepar, hdrarr, 'TAI-BEG', bighdr, func='min'
   sxcombinepar, hdrarr, 'TAI-END', bighdr, func='max'

   sxcombinepar, hdrarr, 'XCHI2', bighdr, func='max', outcard='XCHI2MAX', $
    after='XCHI2'
   sxcombinepar, hdrarr, 'XCHI2', bighdr, func='min', outcard='XCHI2MIN', $
    after='XCHI2'

   sxcombinepar, hdrarr, 'SKYCHI2', bighdr, func='max', outcard='SCHI2MAX', $
    after='SKYCHI2'
   sxcombinepar, hdrarr, 'SKYCHI2', bighdr, func='min', outcard='SCHI2MIN', $
    after='SKYCHI2'

   sxcombinepar, hdrarr, 'WSIGMA', bighdr, func='max', outcard='WSIGMAX', $
    after='WSIGMA'
   sxcombinepar, hdrarr, 'WSIGMA', bighdr, func='min', outcard='WSIGMIN', $
    after='WSIGMA'

   sxcombinepar, hdrarr, 'XSIGMA', bighdr, func='max', outcard='XSIGMAX', $
    after='XSIGMA'
   sxcombinepar, hdrarr, 'XSIGMA', bighdr, func='min', outcard='XSIGMIN', $
    after='XSIGMA'

   ; Add the NGUIDE keywords for all headers of one flavor of CAMERAS
   ; (e.g., for all the 'b1' exposures if the first frame is 'b1'.)
   cardname = 'NGUIDE'
   sxcombinepar, hdrarr[0], cardname, bighdr, func='total'
   cameras0 = sxpar(*(hdrarr[0]), 'CAMERAS')
   for ihdr=1, n_elements(hdrarr)-1 do begin
      if (sxpar(*(hdrarr[ihdr]), 'CAMERAS') EQ cameras0) then $
       sxcombinepar, hdrarr[ihdr], cardname, bighdr, func='total'
   endfor

   ;----------
   ; Use the MJD passed as a keyword, which will typically be for the most
   ; observation, and be consistent with the output file names

   if (keyword_set(mjd)) then $
    sxaddpar, bighdr, 'MJD', mjd

   ; Get the list of MJD's used for these reductions, then convert to a string
   mjdlist = mjdlist[uniq(mjdlist, sort(mjdlist))]
   mjdlist = strtrim(strcompress(string(mjdlist,format='(99a)')),2)
   sxaddpar, bighdr, 'MJDLIST', mjdlist, after='MJD'

   ;----------
   ; Add new header cards

   sxaddpar, bighdr, 'VERSCOMB', idlspec2d_version(), $
    ' Version of idlspec2d for combining multiple spectra', after='VERS2D'
   sxaddpar, bighdr, 'NEXP', nfiles, $
    ' Number of exposures in this file', before='EXPTIME'
   for ifile=0,nfiles-1 do $
    sxaddpar, bighdr, string('EXPID',ifile+1, format='(a5,i2.2)'), label[ifile], $
     ' ID string for exposure '+strtrim(ifile+1,2), before='EXPTIME'
   if (keyword_set(bestexpnum)) then $
    sxaddpar, bighdr, 'BESTEXP', bestexpnum, before='EXPID01'

   sxaddpar, bighdr, 'EXPTIME', min(exptimevec), $
    ' Minimum of exposure times for all cameras'
   for icam=0, ncam-1 do $
    sxaddpar, bighdr, 'NEXP_'+camnames[icam], nexpvec[icam], $
     ' '+camnames[icam]+' camera number of exposures', before='EXPTIME'
   for icam=0, ncam-1 do $
    sxaddpar, bighdr, 'EXPT_'+camnames[icam], exptimevec[icam], $
     ' '+camnames[icam]+' camera exposure time (seconds)', before='EXPTIME'
   sxaddpar, bighdr, 'SPCOADD', systime(), $
    ' SPCOADD finished', after='EXPTIME'

   sxaddpar, bighdr, 'NWORDER', 2, ' Linear-log10 coefficients'
   sxaddpar, bighdr, 'NWORDER', 2, ' Linear-log10 coefficients'
   sxaddpar, bighdr, 'WFITTYPE', 'LOG-LINEAR', ' Linear-log10 dispersion'
   sxaddpar, bighdr, 'COEFF0', wavemin, $
    ' Central wavelength (log10) of first pixel'
   sxaddpar, bighdr, 'COEFF1', binsz, ' Log10 dispersion per pixel'

   sxaddpar, bighdr, 'NAXIS1', n_elements(bestflux)
   sxaddpar, bighdr, 'NAXIS2', nfiber

   spawn, 'uname -n', uname
   sxaddpar, bighdr, 'UNAME', uname[0]

   ;----------
   ; Check for smear exposure used and place info in header

;   smearused = total((finalandmask AND pixelmask_bits('SMEARIMAGE')) NE 0) $
;    GT 0 ? 'T' : 'F'
;   sxaddpar, bighdr, 'SMEARUSE', smearused, ' Smear image used?'

   ;----------
   ; Compute the fraction of bad pixels in total, and on each spectrograph.
   ; Bad pixels are any with SKYMASK(INVVAR)=0, excluding those where
   ; the NODATA bit is set in the pixel mask.

   ifib1 = where(finalplugmap.spectrographid EQ 1, nfib1)
   ifib2 = where(finalplugmap.spectrographid EQ 2, nfib2)
   qbadpix = skymask(finalivar, finalandmask, finalormask) EQ 0 $
    AND (finalandmask AND pixelmask_bits('NODATA')) EQ 0
   if (nfib1 GT 0) then $
    fbadpix1 = total(qbadpix[*,ifib1]) / (nfib1 * nfinalpix) $
   else $
    fbadpix1 = 0
   if (nfib2 GT 0) then $
    fbadpix2 = total(qbadpix[*,ifib2]) / (nfib2 * nfinalpix) $
   else $
    fbadpix2 = 0
   if (nfib1 GT 0 AND nfib2 GT 0) then $
    fbadpix = total(qbadpix[*,[ifib1,ifib2]]) / ((nfib1+nfib2) * nfinalpix) $
   else if (nfib1 GT 0) then $
    fbadpix = fbadpix1 $
   else if (nfib2 GT 0) then $
    fbadpix = fbadpix1 $
   else $
    fbadpix = 0

   sxaddpar, bighdr, 'FBADPIX', fbadpix, ' Fraction of bad pixels'
   sxaddpar, bighdr, 'FBADPIX1', fbadpix1, ' Fraction of bad pixels on spectro-1'
   sxaddpar, bighdr, 'FBADPIX2', fbadpix2, ' Fraction of bad pixels on spectro-2'

   ;----------
   ; Add keywords for IRAF-compatability

   add_iraf_keywords, bighdr, wavemin, binsz

   mkhdr, hdrfloat, finalivar, /image, /extend
   add_iraf_keywords, hdrfloat, wavemin, binsz

   mkhdr, hdrlong, finalandmask, /image, /extend
   add_iraf_keywords, hdrlong, wavemin, binsz

   ;---------------------------------------------------------------------------
   ; Write combined output file
   ;---------------------------------------------------------------------------

   ; First write the file with the flux distortion vectors
   mwrfits, corrimg, distortfitsfile, bighdr, /create

   fulloutname = djs_filepath(outputname, root_dir=combinedir)

   ; HDU #0 is flux
   sxaddpar, bighdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
   mwrfits, finalflux, fulloutname, bighdr, /create

   ; HDU #1 is inverse variance
   sxaddpar, hdrfloat, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
   sxaddpar, hdrfloat, 'EXTNAME', 'IVAR', ' Inverse variance'
   mwrfits, finalivar, fulloutname, hdrfloat

   ; HDU #2 is AND-pixelmask
   sxaddpar, hdrlong, 'EXTNAME', 'ANDMASK', ' AND Mask'
   mwrfits, finalandmask, fulloutname, hdrlong

   ; HDU #3 is OR-pixelmask
   sxaddpar, hdrlong, 'EXTNAME', 'ORMASK', ' OR Mask'
   mwrfits, finalormask, fulloutname, hdrlong

   ; HDU #4 is dispersion map
   sxaddpar, hdrfloat, 'BUNIT', 'pixels'
   sxaddpar, hdrfloat, 'EXTNAME', 'WAVEDISP', ' Wavelength dispersion'
   mwrfits, finaldispersion, fulloutname, hdrfloat

   ; HDU #5 is plugmap
   sxaddpar, hdrplug, 'EXTNAME', 'PLUGMAP', ' Plugmap structure'
   mwrfits, finalplugmap, fulloutname, hdrplug

   ; HDU #6 is the sky
   sxaddpar, hdrsky, 'EXTNAME', 'SKY', ' Subtracted sky flux'
   mwrfits, finalsky, fulloutname, hdrsky

   return
end
;------------------------------------------------------------------------------
