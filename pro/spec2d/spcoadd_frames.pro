;+
; NAME:
;   spcoadd_frames
;
; PURPOSE:
;   Combine several reduced frames of the same objects
;
; CALLING SEQUENCE:
;   spcoadd_frames, filenames, outputname, $
;    binsz, zeropoint, [ nord=, wavemin=, $
;    bkptbin=, window=, maxsep= ]
;
; INPUTS:
;   filenames      - Name(s) of files to combine (written by SPREDUCE)
;   outputname     - Output file name
;
; REQUIRED KEYWORDS:
;
; OPTIONAL KEYWORDS:
;   binsz          - Bin size (in log-10 wavelength) in output spectra;
;                    default to 10d-4, which corresponds to 69.02977415 km/s.
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
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine also outputs original 2048 spectra with mask pixels
;   replaced with their b-spline values. ??? NOPE - This code removed!!!
;
;   This routine can combine data from multiple (different) plug maps.
;   Objects are matched based upon their positions agreeing to 2 arc sec.
;
;   All input files must have the same number of pixels per spectrum,
;   i.e. 2048 wavelength samplings, although those wavelengths can
;   be different.
;
; EXAMPLES:
;
; BUGS:
;   Should only apodize starting with the first/last good pixel of a spectrum.
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   djs_laxisgen()
;   djs_maskinterp()
;   djs_median()
;   idlspec2d_version()
;   mrdfits()
;   pixelmask_bits()
;   sxaddpar
;   sxdelpar
;   sxpar()
;   writefits
;
; INTERNAL SUPPORT PROCEDURES:
;   combine1fiber
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
pro combine1fiber, fullwave, fullspec, fullivar, fullpixelmask, fulldisp, $
 finalwave, bestflux, bestivar, andmask, ormask, bestdisp, $
 nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep

; ??? Return fullcombmask for modifying the masks in the original input files.

   if (NOT keyword_set(nord)) then nord = 3
   if (NOT keyword_set(maxsep)) then maxsep = 2.0 * binsz
   if (NOT keyword_set(bkptbin)) then bkptbin = 1.2 * binsz

   specnum = djs_laxisgen( size(fullwave,/dimens), iaxis=1)

   fullcombmask = bytarr(n_elements(fullspec))

   nfinalpix = N_elements(finalwave)
   bestflux = fltarr(nfinalpix)
   bestivar = bestflux*0.0
   bestdisp = bestflux*0.0

   nonzero = where(fullivar GT 0.0, ngood)

   if (ngood EQ 0) then begin

      splog, 'No good points'
      andmask = pixelmask_bits('NODATA')
      ormask = pixelmask_bits('NODATA')
      return

   endif else begin

      ; Now let's break sorted wavelengths into groups where
      ; pixel separations are larger than maxsep

      isort = nonzero[sort(fullwave[nonzero])]
      wavesort = fullwave[isort]

      padwave = [min(wavesort) - 2.0*maxsep, wavesort, $
       max(wavesort) + 2.0*maxsep]

      ig1 = where(padwave[1:ngood] - padwave[0:ngood-1] GT maxsep, nstart)
      ig2 = where(padwave[2:ngood+1] - padwave[1:ngood] GT maxsep, nend)
      if (nstart NE nend) then $
       message, 'ABORT: Grouping tricks did not work!'

      for igrp=0, nstart-1 do begin

         ss = isort[ig1[igrp] : ig2[igrp]]
         bkpt = 0

         fullbkpt = slatec_splinefit(fullwave[ss], fullspec[ss], coeff, $
          nord=nord, eachgroup=1, maxiter=20, upper=3.0, lower=3.0, $
          bkspace=bkptbin, bkpt=bkpt, invvar=fullivar[ss], mask=bmask, /silent)

         inside = where(finalwave GE min(bkpt) $
          AND finalwave LE max(bkpt), numinside)

         if (numinside EQ 0) then begin
            splog,'WARNING: No wavelengths inside breakpoints'
         endif else if (total(abs(coeff)) EQ 0.0) then begin
            splog,'WARNING: All B-spline coefficients have been set to zero!'
         endif else begin         

            bestflux[inside] = slatec_bvalu(finalwave[inside],fullbkpt,coeff)
fwave = float(finalwave[inside])
print,10^[min(fwave),max(fwave)]

            splog, 'Masked ', fix(total(1-bmask)), ' of', $
             n_elements(bmask), ' pixels'

            ;-----------------------------------------------------------------
            ;  Here replace original flux values of masked pixels with b-spline
            ;  evaluations.

            ireplace = where(bmask EQ 0)

            if (ireplace[0] NE -1) then begin 
               fullspec[ss[ireplace]] = $
                slatec_bvalu(fullwave[ss[ireplace]],fullbkpt,coeff)

               fullpixelmask[ss[ireplace]] = fullpixelmask[ss[ireplace]] OR $
                pixelmask_bits('COMBINEREJ')
            endif

         endelse
         fullcombmask[ss] = bmask

      endfor

      ;---------------------------------------------------------------------
      ; Combine inverse variance and pixel masks.

      andmask = lonarr(nfinalpix) - 1 ; Start with all bits set in AND-mask.
      ormask = lonarr(nfinalpix)

      for j=0, max(specnum) do begin
         these = where(specnum EQ j)
         if (these[0] NE -1) then begin

            inbetween = where(finalwave GE min(fullwave[these]) AND $
                              finalwave LE max(fullwave[these]))
            if (inbetween[0] NE -1) then begin

               ; Conserve inverse variance by doing a linear interpolation
               ; on that quantity.

               result = interpol(fullivar[these] * fullcombmask[these], $
                fullwave[these], finalwave[inbetween])

               ; Grow the fullcombmask below to reject any new sampling
               ; containing even a partial masked pixel.

               smask = interpol(float(fullcombmask[these]), $
                fullwave[these], finalwave[inbetween])
               ibad = where(smask LT 1.0)
               if (ibad[0] NE -1) then result[ibad] = 0

               bestivar[inbetween] = bestivar[inbetween] + result

               lowside = fix((fullwave[these]-finalwave[0])/binsz)
               highside = lowside + 1
               andmask[lowside]  = andmask[lowside] AND fullpixelmask[these]
               andmask[highside] = andmask[highside] AND fullpixelmask[these]
               ormask[lowside]   = ormask[lowside] OR fullpixelmask[these]
               ormask[highside]  = ormask[highside] OR fullpixelmask[these]

               ; Combine the dispersions in the dumbest way possible

               bestdisp[inbetween] = interpol(fulldisp[these], $
                fullwave[these], finalwave[inbetween])

            endif

         endif
      endfor
      splog, 'Medians:', string(format='(f7.2)', $
       djs_median(fullspec[these]))

   endelse

   ;----------
   ; Replace NaN's in combined spectra; this should really never happen

   inff = where(finite(bestflux) EQ 0 OR finite(bestivar) EQ 0)
   if (inff[0] NE -1) then begin
      splog, 'WARNING: NaNs in combined spectra ', N_elements(inff)
      bestflux[inff] = 0.0
      bestivar[inff] = 0.0
   endif

   ;----------
   ; Interpolate over masked pixels, just for aesthetic purposes

   bestflux = djs_maskinterp(bestflux, bestivar EQ 0, /const)

   ;----------
   ; Set the NODATA mask bit wherever there is no good data

   ibad = where(bestivar EQ 0)
   if (ibad[0] NE -1) then begin
      andmask[ibad] = andmask[ibad] AND pixelmask_bits('NODATA')
      ormask[ibad] = ormask[ibad] AND pixelmask_bits('NODATA')
   endif

   ;----------
   ; Replace values of -1 in the AND mask with 0's

   andmask = andmask * (andmask NE -1)

   return
end
;------------------------------------------------------------------------------

pro spcoadd_frames, filenames, outputname, $
 binsz, zeropoint, nord=nord, wavemin=wavemin, $
 bkptbin=bkptbin, window=window, maxsep=maxsep

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

      tempflux = mrdfits(filenames[ifile], 0, hdr)
      tempivar = mrdfits(filenames[ifile], 1)
      temppixmask = mrdfits(filenames[ifile], 2)
      tempwset = mrdfits(filenames[ifile], 3)
      tempdispset = mrdfits(filenames[ifile], 4)
      tempplug = mrdfits(filenames[ifile], 5, structyp='PLUGMAPOBJ')
      if (NOT keyword_set(tempflux)) then $
       message, 'Error reading file ' + filenames[ifile]

      ;----------
      ; Solve for wavelength and lambda-dispersion at each pixel in the image

      traceset2xy, tempwset, junk, tempwave
      traceset2xy, tempdispset, junk, tempdispersion

      dims = size(tempflux, /dimens)
      npix = dims[0]

      ;----------
      ; Determine if this is a blue or red spectrum

      cameras = strtrim(sxpar(hdr, 'CAMERAS'),2)
      icam = where(cameras EQ camnames)
      if (icam[0] EQ -1) then $
       message, 'Unknown camera ' + cameras
      exptimevec[icam] = exptimevec[icam] + sxpar(hdr, 'EXPTIME')

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
         plugmap = [plugmap, tempplug]
      endelse

   endfor

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
   finaldispersion = lonarr(nfinalpix, nfiber)
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
         combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
          pixelmask[*,indx], dispersion[*,indx], $
          finalwave, bestflux, bestivar, bestandmask, bestormask, $
          bestdispersion, nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep

         finalflux[*,ifiber] = bestflux
         finalivar[*,ifiber] = bestivar
         finalandmask[*,ifiber] = bestandmask
         finalormask[*,ifiber] = bestormask
         finaldispersion[*,ifiber] = bestdispersion
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
   ; Write output file

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

   return
end
;------------------------------------------------------------------------------
