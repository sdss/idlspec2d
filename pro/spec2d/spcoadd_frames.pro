;+
; NAME:
;   spcoadd_frames
;
; PURPOSE:
;   Combine several reduced frames of the same objects
;
; CALLING SEQUENCE:
;   spcoadd_frames, filenames, outputname, $
;    fcalibprefix=, [ binsz=, zeropoint=, nord=, wavemin=, $
;    bkptbin=, window=, maxsep=, adderr= ]
;
; INPUTS:
;   filenames      - Name(s) of files to combine (written by SPREDUCE)
;   outputname     - Output file name
;
; REQUIRED KEYWORDS:
;   fcalibprefix   - Prefix for flux-calibration files.
;
; OPTIONAL KEYWORDS:
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
;   divideflat
;   djs_diff_angle()
;   idlspec2d_version()
;   modfits
;   mrdfits()
;   pixelmask_bits()
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

   plate = strtrim(string(sxpar(hdr, 'PLATEID')),2)
   camera = strtrim(sxpar(hdr, 'CAMERAS'),2)
   mjd = strtrim(string(sxpar(hdr, 'MJD')),2)
   seqid =  strtrim(string(sxpar(hdr, 'SEQID')),2)
   expos =  strtrim(string(sxpar(hdr, 'EXPOSURE')),2)

   return, plate+'-'+camera+'-'+mjd+'-'+seqid+'-'+expos
end

;------------------------------------------------------------------------------
pro spcoadd_frames, filenames, outputname, fcalibprefix=fcalibprefix, $
 binsz=binsz, zeropoint=zeropoint, nord=nord, wavemin=wavemin, $
 bkptbin=bkptbin, window=window, maxsep=maxsep, adderr=adderr

   if (NOT keyword_set(binsz)) then binsz = 1.0d-4 $
    else binsz = double(binsz)
   if (NOT keyword_set(zeropoint)) then zeropoint = 3.5D
   if (n_elements(window) EQ 0) then window = 100

   ;----------
   ; Sort filenames such that this list contains first the blue then the red
; Why should this matter to sort them ???

   nfiles = n_elements(filenames)
   if (nfiles EQ 0) then return

   filenames = filenames[sort(filenames)]

   ;---------------------------------------------------------------------------

   camnames = ['b1', 'b2', 'r1', 'r2']
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
      tempflux = mrdfits(filenames[ifile], 0, hdr)
      tempivar = mrdfits(filenames[ifile], 1)
      temppixmask = mrdfits(filenames[ifile], 2)
      tempwset = mrdfits(filenames[ifile], 3)
      tempdispset = mrdfits(filenames[ifile], 4)
      tempplug = mrdfits(filenames[ifile], 5, structyp='PLUGMAPOBJ')
      if (NOT keyword_set(tempflux)) then $
       message, 'Error reading file ' + filenames[ifile]

      ;----------
      ; Add an additional error term equal to ADDERR of the flux.

      if (keyword_set(adderr)) then begin
         gmask = tempivar NE 0 ; =1 for good points
         tempivar = 1.0 / ( 1.0/(tempivar + (1-gmask)) $
          + (adderr * (tempflux>0))^2 ) * gmask
      endif

      ;----------
      ; Read header info

      cameras = strtrim(sxpar(hdr, 'CAMERAS'),2)
      expstr = string(sxpar(hdr,'EXPOSURE'), format='(i8.8)')

      ;----------
      ; Solve for wavelength and lambda-dispersion at each pixel in the image

      traceset2xy, tempwset, junk, tempwave
      traceset2xy, tempdispset, junk, tempdispersion

      dims = size(tempflux, /dimens)
      npix = dims[0]
      nfib = dims[1]

      ;----------
      ; Determine if this is a blue or red spectrum

      icam = (where(cameras EQ camnames))[0]
      if (icam EQ -1) then $
       message, 'Unknown camera ' + cameras
      exptimevec[icam] = exptimevec[icam] + sxpar(hdr, 'EXPTIME')

      ;----------
      ; Apply spectro-photometric calibration

      fcalibfile = fcalibprefix + '-' + camnames[icam] + '.fits'
      calibhdr = headfits(fcalibfile)
      cwavemin = sxpar(calibhdr, 'WAVEMIN')
      cwavemax = sxpar(calibhdr, 'WAVEMAX')
      calibset = mrdfits(fcalibfile, 1)
      calibfac = bspline_valu(tempwave, calibset)
 
      ; Set to bad any pixels whose wavelength is outside the known
      ; flux-calibration region.
      ibad = where(tempwave LT alog10(cwavemin) OR tempwave GT alog10(cwavemax))
      if (ibad[0] NE -1) then calibfac[ibad] = 0

      divideflat, tempflux, tempivar, calibfac, minval=0.05*mean(calibfac)

      ;----------
      ; Apply flux-correction factor between spectro-photometric exposure
      ; and this exposure.

      spectroid = strmid(cameras,1,1)
      corrset = mrdfits('spFluxcorr-'+expstr+'-'+spectroid+'.fits', 1)
      traceset2xy, corrset, tempwave, corrimg

      divideflat, tempflux, tempivar, 1.0/corrimg, minval=0.05*mean(1.0/corrimg)

      ;----------
      ; Apodize the errors

      if (keyword_set(window)) then begin
         swin = window < npix
         tempivar[0:swin-1] = $
          tempivar[0:swin-1] * findgen(swin) / window
         tempivar[npix-swin:npix-1] = $
          tempivar[npix-swin:npix-1] * (swin-1-findgen(swin)) / window
      endif

      ;----------
      ; Concatenate data from all images

      if (ifile EQ 0) then begin
         flux = tempflux
         fluxivar = tempivar
         wave = tempwave
         dispersion = tempdispersion
         pixelmask = temppixmask

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

         ; Append as vectors...
         camerasvec = [camerasvec, cameras]
         label = [label, makelabel(hdr)]
         filenum = [filenum, lonarr(nfib) + ifile]
         plugmap = [plugmap, tempplug]
      endelse

   endfor

   ;----------
   ; Remove the COMBINEREJ bit from the input pixel masks, since this
   ; is the procedure that should set that bit.

   bitval = pixelmask_bits('COMBINEREJ')
   indx = where(pixelmask AND bitval)
   if (indx[0] NE -1) then pixelmask[indx] = pixelmask[indx] - bitval

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

;stop
;set_plot,'x'
;ifib=295
;colorv=['default','red','green','blue','magenta','yellow','cyan','default']
;lam = 10^wave
;indx = where(plugmap.fiberid EQ ifib)
;splot, lam[*,indx[0]], flux[*,indx[0]], color=colorv[0]
;for jnum=1, nfiles-1 do $
; soplot, lam[*,indx[jnum]], flux[*,indx[jnum]], color=colorv[jnum]

   ;---------------------------------------------------------------------------
   ; Construct output data structures, including the wavelength scale
   ;---------------------------------------------------------------------------

   totalpix = (size(flux, /dimens))[0]

   nonzero = where(fluxivar GT 0.0)
   minfullwave = min(wave[nonzero])
   maxfullwave = max(wave[nonzero])

   ; Get max and min wavelength from good pixels

   if (NOT keyword_set(wavemin)) then begin
      spotmin = fix((minfullwave - zeropoint)/binsz) + 1
      spotmax = fix((maxfullwave - zeropoint)/binsz)
      wavemin = spotmin * binsz + zeropoint
      wavemax = spotmax * binsz + zeropoint
   endif else begin
      spotmin = 0
      if (NOT keyword_set(wavemax)) then begin
        spotmax = fix((maxfullwave - wavemin)/binsz)
        wavemax = spotmax * binsz + wavemin
      endif else spotmax = fix((wavemax - wavemin)/binsz)
   endelse

   nfinalpix = spotmax - spotmin + 1
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
          newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
          andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
          nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep

         finalflux[*,ifiber] = bestflux
         finalivar[*,ifiber] = bestivar
         finalandmask[*,ifiber] = bestandmask
         finalormask[*,ifiber] = bestormask
         finaldispersion[*,ifiber] = bestdispersion

         ; The following adds the COMBINEREJ bit to the input pixel masks
         pixelmask[*,indx] = temppixmask
      endif else begin
         splog, 'Fiber', ifiber+1, ' NO DATA'
         finalandmask[*,ifiber] = pixelmask_bits('NODATA')
         finalormask[*,ifiber] = pixelmask_bits('NODATA')
      endelse
   endfor

   ;---------------------------------------------------------------------------
   ; Create the output header
   ;---------------------------------------------------------------------------

   ncoeff = sxpar(hdr, 'NWORDER')
   for i=2, ncoeff-1 do sxdelpar, hdr, 'COEFF'+strtrim(string(i),2)

   ; Now get rid of exposure, and add list of exposures

   sxdelpar, hdr, 'EXPOSURE'
   sxdelpar, hdr, 'SEQID'

   sxaddpar, objhdr, 'VERSCOMB', idlspec2d_version(), $
    'Version of idlspec2d for combining multiple spectra', after='VERS2D'
   sxaddpar, hdr, 'NEXP', nfiles, $
    'Number of exposures in this file', before='EXPTIME'
   for ifile=0,nfiles-1 do $
    sxaddpar, hdr, 'EXPID'+strtrim(string(ifile),2), label[ifile], $
     'ID string for exposure '+strtrim(string(ifile),2), before='EXPTIME'

   sxaddpar, hdr, 'EXPTIME', min(exptimevec), $
    ' Minimum of exposure times for all cameras'
   for icam=0, ncam-1 do $
    sxaddpar, hdr, 'EXPT_'+camnames[icam], exptimevec[icam], $
     ' '+camnames[icam]+' camera exposure time (seconds)', before='EXPTIME'
   sxaddpar, hdr, 'SPCOADD', systime(), $
    ' SPCOADD finished', after='EXPTIME'

   sxaddpar, hdr, 'NWORDER', 2, 'Linear-log10 coefficients'
   sxaddpar, hdr, 'WFITTYPE', 'LOG-LINEAR', 'Linear-log10 dispersion'
   sxaddpar, hdr, 'COEFF0', wavemin, $
    'Central wavelength (log10) of first pixel'
   sxaddpar, hdr, 'COEFF1', binsz, 'Log10 dispersion per pixel'

   sxaddpar, hdr, 'NAXIS1', n_elements(bestflux)
   sxaddpar, hdr, 'NAXIS2', nfiber

   ;----------
   ; Add keywords for anyone dumb enough to use IRAF

   sxaddpar, hdr, 'WAT0_001', 'system=linear'
   sxaddpar, hdr, 'WAT1_001', $
    'wtype=linear label=Wavelength units=Angstroms'
   sxaddpar, hdr, 'CRVAL1', wavemin, 'Iraf zero point'
   sxaddpar, hdr, 'CD1_1', binsz, 'Iraf dispersion'
   sxaddpar, hdr, 'CRPIX1', 1, 'Iraf starting pixel'
   sxaddpar, hdr, 'CTYPE1', 'LINEAR'
   sxaddpar, hdr, 'WCSDIM', nfiber
   sxaddpar, hdr, 'DC-FLAG', 1, 'Log-linear flag'

   ;---------------------------------------------------------------------------
   ; Write combined output file
   ;---------------------------------------------------------------------------

   ; 1st HDU is flux
   mwrfits, finalflux, outputname, hdr, /create

   ; 2nd HDU is inverse variance
   mwrfits, finalivar, outputname

   ; 3rd HDU is AND-pixelmask
   mwrfits, finalandmask, outputname

   ; 4th HDU is OR-pixelmask
   mwrfits, finalormask, outputname

   ; 5th HDU is dispersion map
   mwrfits, finaldispersion, outputname

   ; 6th HDU is plugmap
   mwrfits, finalplugmap, outputname

   ;---------------------------------------------------------------------------
   ; Write the modified pixel masks to the input files
   ;---------------------------------------------------------------------------

   for ifile=0, nfiles-1 do begin
      splog, 'Modifying file #', ifile, ': ', filenames[ifile]
      indx = where(filenum EQ ifile)
      modfits, filenames[ifile], pixelmask[*,indx], exten_no=2
   endfor

   return
end
;------------------------------------------------------------------------------
