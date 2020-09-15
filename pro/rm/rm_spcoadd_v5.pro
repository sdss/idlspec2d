;+
; NAME:
;   rm_spcoadd_v5
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
; Simply scale the overlaping wavelength region of blue/red to a common scale
pro rm_fluxscale_ccd, waveb, fluxb, ivarb, waver, fluxr, ivarr, red=red,blue=blue

   ; common wavelength range
   lam_range = [5900., 6250.]
   d_loglam = 1d-4 
   npix = alog10(lam_range[1]/lam_range[0])/d_loglam
   loglam_vec = alog10(5900.) + d_loglam*findgen(npix+1)

   ind1 = where(waveb ge lam_range[0] and waveb le lam_range[1] $
        and ivarb gt 0, npix1)
   ind2 = where(waver ge lam_range[0] and waver le lam_range[1] $
        and ivarr gt 0, npix2)
   if npix1 gt 100 and npix2 gt 100 then begin
       flux_interp1 = interpol(fluxb[ind1],waveb[ind1],10.D^loglam_vec)
       flux_interp2 = interpol(fluxr[ind2],waver[ind2],10.D^loglam_vec)
       scale = median(flux_interp1/flux_interp2)
       if keyword_set(red) and scale gt 0 then begin
          fluxb = fluxb/scale & ivarb = ivarb*scale^2
       endif 

       if keyword_set(blue) and scale gt 0 then begin
          fluxr = fluxr*scale & ivarr = ivarr/scale^2
       endif

       if not keyword_set(red) and not keyword_set(blue) and scale gt 0 then begin
          fluxb = fluxb*0.5*(1.+1./scale) & ivarb = ivarb/(0.5*(1.+1./scale))^2
          fluxr = fluxr*0.5*(1.+scale) & ivarb = ivarr/(0.5*(1.+scale))^2
       endif

   endif

end

;------------------------------------------------------------------------------
pro rm_spcoadd_v5, spframes, outputname, $
 mjd=mjd, binsz=binsz, zeropoint=zeropoint, nord=nord, $
 wavemin=wavemin, wavemax=wavemax, $
 bkptbin=bkptbin, window=window, maxsep=maxsep, adderr=adderr, $
 docams=camnames, plotsnfile=plotsnfile, combinedir=combinedir, $
 bestexpnum=bestexpnum,nofcorr=nofcorr,nodist=nodist, $
 plates=plates, legacy=legacy, single_spectra=single_spectra

   if (NOT keyword_set(binsz)) then binsz = 1.0d-4 $
    else binsz = double(binsz)
   if (NOT keyword_set(zeropoint)) then zeropoint = 3.5D
   if (n_elements(window) EQ 0) then window = 100
   if (NOT keyword_set(combinedir)) then combinedir=''

   string1 = repstr(outputname,'spField','spFluxdistort')
   string2 = repstr(string1,'.fits','.ps')
   distortfitsfile = djs_filepath(string1, root_dir=combinedir)
   distortpsfile = djs_filepath(string2, root_dir=combinedir)

   ;----------
   ; Sort filenames such that this list contains first the blue then the red

   nfiles = n_elements(spframes)
   if (nfiles EQ 0) then return

   filenames = spframes[sort(spframes)]

   ;---------------------------------------------------------------------------
   if keyword_set(legacy) then begin
      if NOT keyword_set(camnames) then camnames = ['b1', 'r1', 'b2', 'r1']
   endif else begin
      if NOT keyword_set(camnames) then camnames = ['b1', 'r1']
   endelse
   ncam = N_elements(camnames)
   nexpvec = lonarr(ncam)
   exptimevec = fltarr(ncam)
   single_spectra=1
   ;---------------------------------------------------------------------------
   ; Loop through each 2D output and read in the data
   ;---------------------------------------------------------------------------

   ; Start by determining the size of all the files
   npixarr = lonarr(nfiles)
   plugmap_rm=create_struct('CONFIGURATION','','RA0',0.D,'DEC0',0.D,'TAI',0.D,'MJD',0.0,'AIRMASS',0.D,'DATE','')
   if keyword_set(legacy) then begin
     nexp_tmp2 = nfiles/4 ;Get data for each exposure
   endif else begin
     nexp_tmp2 = nfiles/2
   endelse
   rm_plugmap = replicate(plugmap_rm, nexp_tmp2)
   for ifile=0, nfiles-1 do begin
      spframe_read, filenames[ifile], hdr=objhdr
      npixarr[ifile] = sxpar(objhdr,'NAXIS1')
      if ifile LT nexp_tmp2 then begin
        if keyword_set(legacy) or keyword_set(plates) then begin
          rm_plugmap[ifile].configuration=sxpar(objhdr,'PLATEID')
        endif else begin
          rm_plugmap[ifile].configuration=sxpar(objhdr,'CONFIID')
        endelse
        rm_plugmap[ifile].ra0=sxpar(objhdr,'RA')
        rm_plugmap[ifile].dec0=sxpar(objhdr,'DEC')
        rm_plugmap[ifile].tai=sxpar(objhdr,'TAI-BEG')
        rm_plugmap[ifile].mjd=sxpar(objhdr,'MJD')
        rm_plugmap[ifile].airmass=sxpar(objhdr,'AIRMASS')
        mjdt=sxpar(objhdr,'TAI-BEG')/(24.D*3600.D)
        mjd2datelist,mjdt,datelist=date
        rm_plugmap[ifile].date=date
      endif
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

      if not keyword_set(nofcorr) then begin
         nexp_tmp = nfiles/4
         corrfile = djs_filepath(string(camnames[icam], expnum, $
          format='("spFluxcorr-", a2, "-", i8.8, ".fits")'), $
          root_dir=combinedir)
         thisfile = (findfile(corrfile+'*'))[0]
         if (NOT keyword_set(thisfile)) then $
          message,' Could not find flux-corr file ' + corrfile

         aterm = mrdfits(thisfile, 0, corrhdr, /silent)
         bterm = mrdfits(thisfile, 1)
         ; Only apply the fcorr vector to objects with med_SN>3
         tempsn = median(tempflux*sqrt(tempivar),dim=1)
         ind_notapply = where(tempsn lt 1.)
         if ind_notapply[0] ne -1 then aterm[*, ind_notapply] = 1.

         invertcorr = 1. / aterm
         minval = 0.05 / mean(aterm)
         nrownative=(size(tempflux,/dimens))[0]
         divideflat, tempflux, invvar=tempivar, invertcorr[0:nrownative-1,*], minval=minval
         tempflux = tempflux + bterm
         divideflat, tempsky, invertcorr[0:nrownative-1,*], minval=minval
         temppixmask = temppixmask $
          OR (invertcorr LE minval) * pixelmask_bits('BADFLUXFACTOR')

      endif
    
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
         for it = 0, n_elements(tempplug.fiberid)-1 do begin
          if it eq 0 then begin 
            expnumvec1 = expnum
          endif else begin
            expnumvec1 = [expnumvec1, expnum] 
          endelse
         endfor
         expnumvec = expnumvec1
      endif else begin
         ; Append as vectors...
         camerasvec = [camerasvec, cameras]
         label = [label, makelabel(hdr)]
         filenum = [filenum, lonarr(nfib) + ifile]
         plugmap = [plugmap, tempplug]
         for it = 0, n_elements(tempplug.fiberid)-1 do begin
           if it eq 0 then begin
             expnumvec1 = expnum
           endif else begin
             expnumvec1 = [expnumvec1, expnum]
           endelse
         endfor
         expnumvec = [expnumvec, expnumvec1]
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
   pixelmask_rm=pixelmask
   expnumf=expnumvec[uniq(expnumvec, sort(expnumvec))]
   ;print, expnumvec
   ;print, expnumf
   ;indx = (where((plugmap.fiberid EQ 0+1) $
   ;  AND (expnumvec EQ expnumf[0])));[0]
   ;print, indx
   ;print, expnumvec[indx]
   ;indx = (where((plugmap.fiberid EQ 0+1)));[0]
   ;print, indx
   ;print, expnumvec[indx]
   ;exit
   ;----------
   ; Scale the blue and red flux to the same flux level
   ; note that the filenames are sorted as b1[nexp],r1[nexp]
   if keyword_set(legacy) then begin
     nexp_tmp = nfiles/4
     for iexp=0, nexp_tmp - 1 do begin
       for iobj=0L,499L do begin
         ; for b1 and r1
         ifile = iexp
         waveb = wave[0:npixarr[ifile]-1,nobj*ifile+iobj]
         fluxb = flux[0:npixarr[ifile]-1,nobj*ifile+iobj]
         ivarb = fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj]
         ifile = iexp + nexp_tmp*2
         waver = wave[0:npixarr[ifile]-1,nobj*ifile+iobj]
         fluxr = flux[0:npixarr[ifile]-1,nobj*ifile+iobj]
         ivarr = fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj]
         rm_fluxscale_ccd,waveb,fluxb,ivarb,waver,fluxr,ivarr,/red
         ifile = iexp
         flux[0:npixarr[ifile]-1,nobj*ifile+iobj] = fluxb
         fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj] = ivarb
         ifile = iexp + nexp_tmp*2
         flux[0:npixarr[ifile]-1,nobj*ifile+iobj] = fluxr
         fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj] = ivarr
         ; for b2 and r2
         ifile = iexp + nexp_tmp*1
         waveb = wave[0:npixarr[ifile]-1,nobj*ifile+iobj]
         fluxb = flux[0:npixarr[ifile]-1,nobj*ifile+iobj]
         ivarb = fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj]
         ifile = iexp + nexp_tmp*3
         waver = wave[0:npixarr[ifile]-1,nobj*ifile+iobj]
         fluxr = flux[0:npixarr[ifile]-1,nobj*ifile+iobj]
         ivarr = fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj]
         rm_fluxscale_ccd,waveb,fluxb,ivarb,waver,fluxr,ivarr,/red
         ifile = iexp + nexp_tmp*1
         flux[0:npixarr[ifile]-1,nobj*ifile+iobj] = fluxb
         fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj] = ivarb
         ifile = iexp + nexp_tmp*3
         flux[0:npixarr[ifile]-1,nobj*ifile+iobj] = fluxr
         fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj] = ivarr
       endfor
     endfor
   endif else begin
     nexp_tmp = nfiles/2;#4
     for iexp=0, nexp_tmp - 1 do begin
         for iobj=0L,499L do begin
            ; for b1 and r1
            ifile = iexp
            waveb = wave[0:npixarr[ifile]-1,nobj*ifile+iobj]
            fluxb = flux[0:npixarr[ifile]-1,nobj*ifile+iobj]
            ivarb = fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj]
            ifile = iexp + nexp_tmp;*4
            waver = wave[0:npixarr[ifile]-1,nobj*ifile+iobj]
            fluxr = flux[0:npixarr[ifile]-1,nobj*ifile+iobj]
            ivarr = fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj]
            rm_fluxscale_ccd,waveb,fluxb,ivarb,waver,fluxr,ivarr,/red
            ifile = iexp
            flux[0:npixarr[ifile]-1,nobj*ifile+iobj] = fluxb
            fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj] = ivarb
            ifile = iexp + nexp_tmp;*4
            flux[0:npixarr[ifile]-1,nobj*ifile+iobj] = fluxr
            fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj] = ivarr
        endfor
     endfor
   endelse
   ;-----------

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

   nfiber = max(plugmap.fiberid)
   ;nfiber = 2 * nfib

   ;finalflux = fltarr(nfinalpix, nfiber)
   ;finalivar = fltarr(nfinalpix, nfiber)
   ;finalandmask = lonarr(nfinalpix, nfiber)
   ;finalormask = lonarr(nfinalpix, nfiber)
   ;finaldispersion = fltarr(nfinalpix, nfiber)
   ;finalsky = fltarr(nfinalpix, nfiber)
   ;finalplugmap = replicate(plugmap[0], nfiber)
   
   finalflux_rm = fltarr(nfinalpix, nfiber, nexp_tmp)
   finalivar_rm = fltarr(nfinalpix, nfiber, nexp_tmp)
   finalandmask_rm = lonarr(nfinalpix, nfiber, nexp_tmp)
   finalormask_rm = lonarr(nfinalpix, nfiber, nexp_tmp)
   finaldispersion_rm = fltarr(nfinalpix, nfiber, nexp_tmp)
   finalsky_rm = fltarr(nfinalpix, nfiber, nexp_tmp)
   finalplugmap_rm = replicate(plugmap[0], nfiber, nexp_tmp)
   finalra_rm = fltarr(nfiber, nexp_tmp)
   finaldec_rm = fltarr(nfiber, nexp_tmp)
   mjds_rm = lonarr(nfiber, nexp_tmp)
   config_rm = lonarr(nfiber, nexp_tmp)
   

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
      for iexp=0, nexp_tmp - 1 do begin
       ; Find the first occurance of fiber number IFIBER+1
       indx = (where((plugmap.fiberid EQ ifiber+1) $
         AND (expnumvec EQ expnumf[iexp])));[0]
       if (indx[0] NE -1) then begin
         splog, 'Coadd read & blue exposure ',expnumf[iexp]
         splog, 'Fiber', ifiber+1, ' ', plugmap[indx[0]].objtype, $
           plugmap[indx[0]].mag, format = '(a, i5.4, a, a, f6.2, 5f6.2)'
         finalplugmap_rm[ifiber,iexp] = plugmap[indx[0]]
         ; Identify all objects with the same XFOCAL,YFOCAL plate position, and
         ; combine all these objects into a single spectrum.
         ; If all pluggings are identical, then this will always be
         ; the same fiber ID.
         ; Also, insist that the object type is not 'NA', which would
         ; occur for unplugged fibers. <--- Disable this for BOSS ???
         ;indx = where(abs(plugmap.xfocal - plugmap[indx].xfocal) LT 0.0001 $
         ;  AND abs(plugmap.yfocal - plugmap[indx].yfocal) LT 0.0001)
         ;          AND strtrim(plugmap.objtype,2) NE 'NA')
       endif
       if (indx[0] NE -1) then begin
        
         temppixmask = pixelmask_rm[*,indx]
         combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
           finalmask=temppixmask, indisp=dispersion[*,indx], $
           skyflux=skyflux[*,indx], $
           newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
           andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
           newsky=bestsky, $
           nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
           maxiter=50, upper=3.0, lower=3.0, maxrej=1
         finalflux_rm[*,ifiber,iexp] = bestflux
         finalivar_rm[*,ifiber,iexp] = bestivar
         finalandmask_rm[*,ifiber,iexp] = bestandmask
         finalormask_rm[*,ifiber,iexp] = bestormask
         finaldispersion_rm[*,ifiber,iexp] = bestdispersion
         finalsky_rm[*,ifiber,iexp] = bestsky
         ; The following adds the COMBINEREJ bit to the input pixel masks
         pixelmask_rm[*,indx] = temppixmask
         ratemp=plugmap[indx].ra
         finalra_rm[ifiber,iexp]=ratemp[0]
         dectemp=plugmap[indx].dec
         finaldec_rm[ifiber,iexp]=dectemp[0]
         mjds_rm[ifiber,iexp]=rm_plugmap[iexp].mjd
         ; use expuse number instad of configuration number for legacy
         config_rm[ifiber,iexp]=rm_plugmap[iexp].configuration
       endif else begin
         splog, 'Fiber', ifiber+1, ' NO DATA'
         finalandmask_rm[*,ifiber,iexp] = pixelmask_bits('NODATA')
         finalormask_rm[*,ifiber,iexp] = pixelmask_bits('NODATA')
       endelse
      endfor
   ;   ; Find the first occurance of fiber number IFIBER+1
   ;   indx = (where(plugmap.fiberid EQ ifiber+1));[0]
   ;   if (indx[0] NE -1) then begin
   ;      splog, 'Coadd all the exposures'
   ;      splog, 'Fiber', ifiber+1, ' ', plugmap[indx[0]].objtype, $
   ;       plugmap[indx[0]].mag, format = '(a, i5.4, a, a, f6.2, 5f6.2)'
   ;      finalplugmap[ifiber] = plugmap[indx[0]]
   ;      ;Check this part, this will no longer the case for the BHM
   ;      ; Identify all objects with the same XFOCAL,YFOCAL plate position, and
   ;      ; combine all these objects into a single spectrum.
   ;      ; If all pluggings are identical, then this will always be
   ;      ; the same fiber ID.
   ;      ; Also, insist that the object type is not 'NA', which would
   ;      ; occur for unplugged fibers. <--- Disable this for BOSS ???
   ;      ;indx = where(abs(plugmap.xfocal - plugmap[indx].xfocal) LT 0.0001 $
   ;      ; AND abs(plugmap.yfocal - plugmap[indx].yfocal) LT 0.0001)
   ;       AND strtrim(plugmap.objtype,2) NE 'NA')
   ;   endif
   ;   if (indx[0] NE -1) then begin
   ;      temppixmask = pixelmask[*,indx]
   ;      combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
   ;       finalmask=temppixmask, indisp=dispersion[*,indx], $
   ;       skyflux=skyflux[*,indx], $
   ;       newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
   ;       andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
   ;       newsky=bestsky, $
   ;       nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
   ;       maxiter=50, upper=3.0, lower=3.0, maxrej=1
   ;      finalflux[*,ifiber] = bestflux
   ;      finalivar[*,ifiber] = bestivar
   ;      finalandmask[*,ifiber] = bestandmask
   ;      finalormask[*,ifiber] = bestormask
   ;      finaldispersion[*,ifiber] = bestdispersion
   ;      finalsky[*,ifiber] = bestsky
   ;      ; The following adds the COMBINEREJ bit to the input pixel masks
   ;      pixelmask[*,indx] = temppixmask
   ;   endif else begin
   ;      splog, 'Fiber', ifiber+1, ' NO DATA'
   ;      finalandmask[*,ifiber] = pixelmask_bits('NODATA')
   ;      finalormask[*,ifiber] = pixelmask_bits('NODATA')
   ;   endelse
   endfor
   
   ;Set a list of targets from its coordinates, this block code considers
   ;the posibility to observe the same target at a diferent fiber in a
   ;diferent FPS configuartion
   brake=0
   indx0=0
   ra_rm=plugmap.ra
   dec_rm=plugmap.dec
   ra_tp=ra_rm;plugmap.ra
   dec_tp=dec_rm;plugmap.dec
   ;indx_tar=list() ; I would prefer to use the list entity but 
   ;idl version 7.7 doesnt have include yet the list definition
   while brake eq 0 do begin
     indx1=where(ra_rm eq ra_tp[0])
     indx1=indx1[0]
     nt1=where(abs(ra_rm - ra_rm[indx1])*3600 LE 0.5 $
       AND abs(dec_rm - dec_rm[indx1])*3600 LE 0.5)
     nt2=where(abs(ra_tp - ra_tp[0])*3600 LE 0.5 $
       AND abs(dec_tp - dec_tp[0])*3600 LE 0.5)
     if (nt1[0] NE -1) then begin
       if indx0 eq 0 then begin
         ;indx_tar.add,nt1
         indx_tar=[nt1]
       endif else begin
         indx_tar=[indx_tar,'-10',nt1]
       endelse
       if n_elements(ra_tp) gt n_elements(nt2) then begin
         remove,nt2,ra_tp,dec_tp
         indx1=where(ra_rm eq ra_tp[0])
         indx1=indx1[0]
       endif else begin
         brake=1
       endelse
       indx0+=1
     endif
   endwhile
   indx_tar=['-10',indx_tar,'-10']
   nt=where(indx_tar eq -10,ntarget)
   ntarget=ntarget-1
   finalflux = fltarr(nfinalpix, ntarget)
   finalivar = fltarr(nfinalpix, ntarget)
   finalandmask = lonarr(nfinalpix, ntarget)
   finalormask = lonarr(nfinalpix, ntarget)
   finaldispersion = fltarr(nfinalpix, ntarget)
   finalsky = fltarr(nfinalpix, ntarget)
   finalplugmap = replicate(plugmap[0], ntarget)
   mjds = lonarr(ntarget)
   final_ra = dblarr(ntarget)
   final_dec = dblarr(ntarget)
   indx_target=intarr(ntarget)
   nexp_target=intarr(ntarget)
   indx_target_s=replicate(create_struct('target_index',0),ntarget)
   nexp_target_s=replicate(create_struct('nexp',0),ntarget)
   struct_assign, {fiberid: 0L}, finalplugmap ; Zero out all elements in this
   ; FINALPLUGMAP structure.
   for itarget=0, ntarget-1 do begin
      indx=indx_tar[nt[itarget]+1:nt[itarget+1]-1]
      if (indx[0] NE -1) then begin
         splog, 'Coadd all the exposures with the same coordinates'
         splog, 'Target', itarget+1, ' ', plugmap[indx[0]].objtype, $
          plugmap[indx[0]].mag, format = '(a, i5.4, a, a, f6.2, 5f6.2)'
         finalplugmap[itarget] = plugmap[indx[0]]
         mjds[itarget]=mjds_rm[indx[0]]
         final_ra[itarget]=plugmap[indx[0]].ra
         final_dec[itarget]=plugmap[indx[0]].dec
      ;      ;Check this part, this will no longer the case for the BHM
      ;      ; Identify all objects with the same XFOCAL,YFOCAL plate position, and
      ;      ; combine all these objects into a single spectrum.
      ;      ; If all pluggings are identical, then this will always be
      ;      ; the same fiber ID.
      ;      ; Also, insist that the object type is not 'NA', which would
      ;      ; occur for unplugged fibers. <--- Disable this for BOSS ???
      ;      ;indx = where(abs(plugmap.xfocal - plugmap[indx].xfocal) LT 0.0001 $
      ;      ; AND abs(plugmap.yfocal - plugmap[indx].yfocal) LT 0.0001)
      ;      ; AND strtrim(plugmap.objtype,2) NE 'NA')
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
         finalflux[*,itarget] = bestflux
         finalivar[*,itarget] = bestivar
         finalandmask[*,itarget] = bestandmask
         finalormask[*,itarget] = bestormask
         finaldispersion[*,itarget] = bestdispersion
         finalsky[*,itarget] = bestsky
         ; The following adds the COMBINEREJ bit to the input pixel masks
         pixelmask[*,indx] = temppixmask
         indx_target[itarget]=itarget+1
         ;if keyword_set(legacy) then begin
         ;  nexp_target[itarget]=n_elements(indx)/2
         ;endif else begin
           nexp_target[itarget]=n_elements(indx)/2
         ;endelse
      endif else begin
         splog, 'Target', itarget+1, ' NO DATA'
         finalandmask[*,itarget] = pixelmask_bits('NODATA')
         finalormask[*,itarget] = pixelmask_bits('NODATA')
         indx_target[itarget]=itarget
         nexp_target[itarget]=0
      endelse
   endfor
   ;print, indx_tar[nt[0]+1:nt[1]-1] mod nfiber
   ;print, indx_tar[nt[0]+1:nt[1]-1]/nfiber mod nexp_tmp
   indx_target_s.target_index=indx_target
   nexp_target_s.nexp=nexp_target
   finalplugmap=struct_addtags(finalplugmap,indx_target_s)
   finalplugmap=struct_addtags(finalplugmap,nexp_target_s)
   ;----------
   ; Modify the 1st file's header to use for the combined plate header.

   bighdr = *hdrarr[0]

   ;---------------------------------------------------------------------------
   ; FLUX DISTORTION IMAGE
   ;---------------------------------------------------------------------------

   ; Compute the flux distortion image
   if not keyword_set(nodist) then begin
      invcorrimg_rm = fltarr(nfinalpix, nfiber, nexp_tmp)
      minicorrval_rm = fltarr(nexp_tmp)
      splog, 'Compute the flux distortion image for each exposure'
      for iexp=0, nexp_tmp - 1 do begin
        splog, 'EXPOSURE number ', expnumf[iexp]
        corrimg = flux_distortion(finalflux_rm[*,*,iexp], finalivar_rm[*,*,iexp], $
          finalandmask_rm[*,*,iexp], finalormask_rm[*,*,iexp], $
          plugmap=finalplugmap_rm[*,iexp], loglam=finalwave, plotfile=distortpsfile, hdr=bighdr, $
          legacy=legacy)
        igood = where(finalivar_rm[*,*,iexp] GT 0)
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
        platesn, finalflux_rm[*,*,iexp], finalivar_rm[*,*,iexp], $
          finalandmask_rm[*,*,iexp], finalplugmap_rm[*,iexp], finalwave, $
          hdr=bighdr, legacy=legacy, plotfile=djs_filepath(repstr(plotsnfile+'.orig','X',string(iexp,format='(i2.2)')), root_dir=combinedir)
        splog, prelog=''
        ; Apply this flux-distortion to the final, co-added fluxes.
        invcorrimg = 1. / corrimg
        minicorrval = 0.05 / mean(corrimg)
        invcorrimg_rm[*,*,iexp]=invcorrimg
        minicorrval_rm[iexp]=minicorrval
        final_flux=finalflux_rm[*,*,iexp]
        final_ivar=finalivar_rm[*,*,iexp]
        final_sky=finalsky_rm[*,*,iexp]
        finaland_mask=finalandmask_rm[*,*,iexp]
        finalor_mask=finalormask_rm[*,*,iexp]
        divideflat, final_flux, invvar=final_ivar, invcorrimg, minval=minicorrval
        divideflat, final_sky, invcorrimg, minval=minicorrval
        finalandmask_t = finaland_mask $
          OR (invcorrimg LE minicorrval) * pixelmask_bits('BADFLUXFACTOR')
        finalormask_t = finalor_mask $
          OR (invcorrimg LE minicorrval) * pixelmask_bits('BADFLUXFACTOR')
        ; Plot S/N and throughput **after** this distortion-correction.
        ; (This over-writes header cards written in the first call.)
        splog, prelog='Final'
        finalflux_rm[*,*,iexp]=final_flux
        finalivar_rm[*,*,iexp]=final_ivar
        finalsky_rm[*,*,iexp]=final_sky
        finalandmask_rm[*,*,iexp]=finaland_mask
        finalormask_rm[*,*,iexp]=finalor_mask
        platesn, finalflux_rm[*,*,iexp], finalivar_rm[*,*,iexp], $
          finalandmask_rm[*,*,iexp], finalplugmap_rm[*,iexp], finalwave, $
          hdr=bighdr, legacy=legacy, plotfile=djs_filepath(repstr(plotsnfile,'X',string(iexp,format='(i2.2)')), root_dir=combinedir)
        splog, prelog=''
      endfor
      splog, 'Compute the flux distortion image for all exposures'  
      corrimg = flux_distortion(finalflux, finalivar, finalandmask, finalormask, $
       plugmap=finalplugmap, loglam=finalwave, plotfile=distortpsfile, hdr=bighdr, $
       legacy=legacy)
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
       hdr=bighdr, legacy=legacy, plotfile=djs_filepath(repstr(plotsnfile+'.orig','-X',''), root_dir=combinedir)
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
       hdr=bighdr, legacy=legacy, plotfile=djs_filepath(repstr(plotsnfile,'-X',''), root_dir=combinedir)
      splog, prelog=''

   endif
   ;---------------------------------------------------------------------------
   ; Write the corrected spCFrame files.
   ; All the fluxes + their errors are calibrated.
   ; The wavelengths + dispersions are converted from trace sets to 2D images.
   ; The pixel mask has the COMBINEREJ bit set.
   ;---------------------------------------------------------------------------

   for ifile=0, nfiles-1 do begin
      thisfile = fileandpath(filenames[ifile], path=thispath)
      ; MODIFIED BY YUE SHEN TO WRITE THE SPCFRAME FILES IN THE SAME COMBINEDIR
      ; 01-17-2014
      thispath = combinedir

      thisfile = djs_filepath(repstr(thisfile,'spFrame','spCFrame'), $
       root_dir=thispath)
      splog, 'Writing file #', ifile, ': ', thisfile, prename=filenames[ifile]
      indx = where(filenum EQ ifile, nthis)

      hdr = *hdrarr[ifile]
      sxaddpar, hdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'

      ; Apply the flux-distortion image to each individual frame, by
      ; interpolating off the full wavelength-scale distortion image
      ; onto the wavelength mapping of each individual exposure+CCD.
      if not keyword_set(nodist) then begin
         ;print, nthis, ifile mod 4
         for i=0L, nthis-1 do begin
            thisflux1 = flux[*,indx[i]]
            thisivar1 = fluxivar[*,indx[i]]
            thissky1 = skyflux[*,indx[i]]
            j = plugmap[indx[i]].fiberid - 1
            tt=ifile mod nexp_tmp
            thisicorr = interpol(invcorrimg_rm[*,j,tt], finalwave, wave[*,indx[i]])
            divideflat, thisflux1, invvar=thisivar1, thisicorr, minval=minicorrval
            flux[*,indx[i]] = thisflux1
            fluxivar[*,indx[i]] = thisivar1

            divideflat, thissky1, thisicorr, minval=minicorrval
            skyflux[*,indx[i]] = thissky1
         endfor

      endif

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
   sxaddpar, bighdr, 'FIELDID', strmid(outputname,8,5)
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
   ; Average together some of the fields from the individual headers. fieldid

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
   bighdr_rm=bighdr
   sxaddpar, bighdr_rm, 'NAXIS3', nexp_tmp, ''

   ;----------
   ; Add keywords for IRAF-compatability

   add_iraf_keywords, bighdr, wavemin, binsz

   mkhdr, hdrfloat, finalivar, /image, /extend
   add_iraf_keywords, hdrfloat, wavemin, binsz
   
   mkhdr, hdrfloat_rm, finalivar_rm, /image, /extend
   add_iraf_keywords, hdrfloat_rm, wavemin, binsz

   mkhdr, hdrlong, finalandmask, /image, /extend
   add_iraf_keywords, hdrlong, wavemin, binsz
   
   mkhdr, hdrlong_rm, finalandmask_rm, /image, /extend
   add_iraf_keywords, hdrlong_rm, wavemin, binsz

   ;---------------------------------------------------------------------------
   ; Write combined output file
   ;---------------------------------------------------------------------------

   ; First write the file with the flux distortion vectors
   if not keyword_set(nodist) then $
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
   
   ;; HDU #7 is flux
   ;sxaddpar, bighdr_rm, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
   ;mwrfits, finalflux_rm, fulloutname, bighdr_rm

   ;; HDU #8 is inverse variance
   ;sxaddpar, hdrfloat_rm, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
   ;sxaddpar, hdrfloat_rm, 'EXTNAME', 'IVAR', ' Inverse variance'
   ;mwrfits, finalivar_rm, fulloutname, hdrfloat_rm

   ;; HDU #9 is AND-pixelmask
   ;sxaddpar, hdrlong_rm, 'EXTNAME', 'ANDMASK', ' AND Mask'
   ;mwrfits, finalandmask, fulloutname, hdrlong_rm

   ;; HDU #10 is OR-pixelmask
   ;sxaddpar, hdrlong_rm, 'EXTNAME', 'ORMASK', ' OR Mask'
   ;mwrfits, finalormask_rm, fulloutname, hdrlong_rm

   ;; HDU #11 is dispersion map
   ;sxaddpar, hdrfloat_rm, 'BUNIT', 'pixels'
   ;sxaddpar, hdrfloat_rm, 'EXTNAME', 'WAVEDISP', ' Wavelength dispersion'
   ;mwrfits, finaldispersion_rm, fulloutname, hdrfloat_rm

   ;; HDU #12 is the sky
   ;sxaddpar, hdrsky_rm, 'EXTNAME', 'SKY', ' Subtracted sky flux'
   ;mwrfits, finalsky_rm, fulloutname, hdrsky_rm
   
   ;; HDU #13 is rm plugmap
   ;sxaddpar, hdrplug, 'EXTNAME', 'PLUGMAP', ' Plugmap structure'
   ;mwrfits, rm_plugmap, fulloutname, hdrplug
   ;delvar, hdrplug
   if keyword_set(single_spectra) then begin
   ;writing each individual coadd spectrum on the field
   sxdelpar, bighdr, 'NAXIS2'
   spawn,'mkdir -p '+combinedir+'coadd'
   for itarget=0, ntarget-1 do begin
    
     finalvalues=replicate(create_struct('flux',0.0),n_elements(finalwave))
     values_t=replicate(create_struct('loglam',0.0),n_elements(finalwave))
     finalvalues=struct_addtags(finalvalues,values_t)
     values_t=replicate(create_struct('ivar',0.0),n_elements(finalwave))
     finalvalues=struct_addtags(finalvalues,values_t)
     values_t=replicate(create_struct('and_mask',0.0),n_elements(finalwave))
     finalvalues=struct_addtags(finalvalues,values_t)
     values_t=replicate(create_struct('or_mask',0.0),n_elements(finalwave))
     finalvalues=struct_addtags(finalvalues,values_t)
     values_t=replicate(create_struct('wdisp',0.0),n_elements(finalwave))
     finalvalues=struct_addtags(finalvalues,values_t)
     values_t=replicate(create_struct('sky',0.0),n_elements(finalwave))
     finalvalues=struct_addtags(finalvalues,values_t)
     finalvalues.flux=finalflux[*,itarget]
     finalvalues.loglam=finalwave
     finalvalues.ivar=finalivar[*,itarget]
     finalvalues.and_mask=finalandmask[*,itarget]
     finalvalues.or_mask=finalormask[*,itarget]
     finalvalues.wdisp=finaldispersion[*,itarget]
     finalvalues.sky=finalsky[*,itarget]
     if keyword_set(legacy) or keyword_set(plates) then begin
        if keyword_set(legacy) then begin
           targid_tar=string(finalplugmap[itarget].fiberid,format='(i4.4)');string(itarget,format='(i3.3)');strtrim(strcompress(string(itarget,format='(99a)')),2)
        endif else begin
           targid_tar=string(finalplugmap[itarget].targetid,format='(i10.10)')
        endelse
     endif else begin   
        targid_tar=finalplugmap[itarget].targetid
     endelse
     sxaddpar, bighdr, 'PLUG_RA', final_ra[itarget], $
       ' RA of Target'
     sxaddpar, bighdr, 'PLUG_DEC', final_dec[itarget], $
       ' DEC of Target'
     if keyword_set(legacy) or keyword_set(plates) then begin
       if keyword_set(mjd) then begin
         thismjd=mjd
       endif else begin
         thismjd=mjds[itarget]
       endelse
     endif else begin
       thismjd=mjds[itarget]
     endelse
     thismjd=strtrim(strcompress(string(thismjd,format='(99a)')),2)
     coadddir=combinedir+'coadd/'+thismjd
     spawn,'mkdir -p '+coadddir
     ;coaddname = repstr(repstr(outputname,'spField','spSpec'),'.fits', $ 
     ; '-'+string(itarget,format='(i3.3)')+'.fits')
     coaddname = repstr(repstr(outputname,'spField','spSpec'),'.fits', $
      '-'+targid_tar+'.fits') 
     fulloutname_coadd = djs_filepath(coaddname, root_dir=coadddir)
     ; HDU # 0 header
     mwrfits, junk_d, fulloutname_coadd, bighdr, /create
     ; HDU # 1 header
     if itarget eq 0 then begin
       sxaddpar, coadd_val, 'EXTNAME', 'COADD', ' Coadded spectrum'
     endif
     mwrfits, finalvalues, fulloutname_coadd, coadd_val
     sxdelpar, coadd_val, 'COMMENT'
     ;sxdelpar, coadd_val, 'EXTNAME'
     ;delvar,coadd_val     
     ;indx=indx_tar[nt[itarget]+1:nt[itarget+1]-1]
     ;if (indx[0] NE -1) then begin
     ;endif
     
     ;; HDU #0 is flux
     ;sxaddpar, bighdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
     ;mwrfits, finalflux[*,itarget], fulloutname_coadd, bighdr, /create

     ;; HDU #1 is inverse variance
     ;sxaddpar, hdrfloat, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
     ;sxaddpar, hdrfloat, 'EXTNAME', 'IVAR', ' Inverse variance'
     ;mwrfits, finalivar[*,itarget], fulloutname_coadd, hdrfloat

     ;; HDU #2 is AND-pixelmask
     ;sxaddpar, hdrlong, 'EXTNAME', 'ANDMASK', ' AND Mask'
     ;mwrfits, finalandmask[*,itarget], fulloutname_coadd, hdrlong

     ;; HDU #3 is OR-pixelmask
     ;sxaddpar, hdrlong, 'EXTNAME', 'ORMASK', ' OR Mask'
     ;mwrfits, finalormask[*,itarget], fulloutname_coadd, hdrlong

     ;; HDU #4 is dispersion map
     ;sxaddpar, hdrfloat, 'BUNIT', 'pixels'
     ;sxaddpar, hdrfloat, 'EXTNAME', 'WAVEDISP', ' Wavelength dispersion'
     ;mwrfits, finaldispersion[*,itarget], fulloutname_coadd, hdrfloat

     ; HDU #5 is plugmap
     if itarget eq 0 then begin
       sxaddpar, hdrplug, 'EXTNAME', 'PLUGMAP', ' Plugmap structure'
     endif
     mwrfits, finalplugmap[itarget], fulloutname_coadd, hdrplug
     sxdelpar, hdrplug, 'COMMENT'
     ;delvar,hdrplug;
     ;; HDU #6 is the sky
     ;sxaddpar, hdrsky, 'EXTNAME', 'SKY', ' Subtracted sky flux'
     ;mwrfits, finalsky[*,itarget], fulloutname_coadd, hdrsky
   ;endfor
   ;
   ;writing each individual single exposure spectrum on the field
   ;spawn,'mkdir -p '+combinedir+'/single'
   for ifiber=0, nfiber-1 do begin
     for iexp=0, nexp_tmp - 1 do begin
       if keyword_set(legacy) or keyword_set(plates) then begin
         if keyword_set(legacy) then begin
            targid_rm=string(finalplugmap_rm[ifiber,iexp].fiberid,format='(i4.4)');targid_tar
         endif else begin
            targid_rm=string(finalplugmap_rm[ifiber,iexp].targetid,format='(i10.10)');targid_tar
         endelse
       endif else begin
         targid_rm=finalplugmap_rm[ifiber,iexp].targetid
       endelse
       if targid_rm eq targid_tar then begin
         finalvalues_rm=replicate(create_struct('flux',0.0),n_elements(finalwave))
         values_t=replicate(create_struct('loglam',0.0),n_elements(finalwave))
         finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
         values_t=replicate(create_struct('ivar',0.0),n_elements(finalwave))
         finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
         values_t=replicate(create_struct('and_mask',0.0),n_elements(finalwave))
         finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
         values_t=replicate(create_struct('or_mask',0.0),n_elements(finalwave))
         finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
         values_t=replicate(create_struct('wdisp',0.0),n_elements(finalwave))
         finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
         values_t=replicate(create_struct('sky',0.0),n_elements(finalwave))
         finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
         finalvalues_rm.flux=finalflux_rm[*,ifiber,iexp]
         finalvalues_rm.loglam=finalwave
         finalvalues_rm.ivar=finalivar_rm[*,ifiber,iexp]
         finalvalues_rm.and_mask=finalandmask_rm[*,ifiber,iexp]
         finalvalues_rm.or_mask=finalormask_rm[*,ifiber,iexp]
         finalvalues_rm.wdisp=finaldispersion_rm[*,ifiber,iexp]
         finalvalues_rm.sky=finalsky_rm[*,ifiber,iexp]
         ; HDU # N header
         thisconf=config_rm[ifiber,iexp]
         thisconf=string(thisconf,format='(i6.6)')
         sxaddpar, indv_val, 'EXTNAME', 'CONFIG_'+thisconf, ' Single exposure spectrum'
         mwrfits, finalvalues_rm, fulloutname_coadd, indv_val
         sxdelpar, indv_val, 'COMMENT'
         ;sxdelpar, indv_val, 'EXTNAME'
         ;delvar,indv_val
         
       endif
       
       ;sxaddpar, bighdr, 'RA', finalra_rm[ifiber,iexp], $
       ;  ' RA of the Target'
       ;sxaddpar, bighdr, 'DEC', finaldec_rm[ifiber,iexp], $
       ;  ' DEC of the Target'
       ;thismjd=mjds_rm[ifiber,iexp]
       ;thismjd=strtrim(strcompress(string(thismjd,format='(99a)')),2)
       ;thisconf=config_rm[ifiber,iexp]
       ;thisconf=string(thisconf,format='(i6.6)')
       ;singledir=combinedir+'/single/'+thismjd
       ;spawn,'mkdir -p '+singledir
       ;singlename = repstr(repstr(repstr(outputname,'spField','spConfig'),'.fits', $
       ; '-'+string(ifiber,format='(i3.3)')+'-'+string(iexp,format='(i2.2)')+'.fits'),strmid(outputname,8,4),thisconf)
       ;fulloutname_single = djs_filepath(singlename, root_dir=singledir)
       ;print, fulloutname_single
       ;; HDU #0 is flux
       ;sxaddpar, bighdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
       ;mwrfits, finalflux_rm[*,ifiber,iexp], fulloutname_single, bighdr, /create
       ;
       ;; HDU #1 is inverse variance
       ;sxaddpar, hdrfloat, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
       ;sxaddpar, hdrfloat, 'EXTNAME', 'IVAR', ' Inverse variance'
       ;mwrfits, finalivar_rm[*,ifiber,iexp], fulloutname_single, hdrfloat
       ;
       ;; HDU #2 is AND-pixelmask
       ;sxaddpar, hdrlong, 'EXTNAME', 'ANDMASK', ' AND Mask'
       ;mwrfits, finalandmask_rm[*,ifiber,iexp], fulloutname_single, hdrlong
       ;
       ;; HDU #3 is OR-pixelmask
       ;sxaddpar, hdrlong, 'EXTNAME', 'ORMASK', ' OR Mask'
       ;mwrfits, finalormask_rm[*,ifiber,iexp], fulloutname_single, hdrlong
       ;
       ;; HDU #4 is dispersion map
       ;sxaddpar, hdrfloat, 'BUNIT', 'pixels'
       ;sxaddpar, hdrfloat, 'EXTNAME', 'WAVEDISP', ' Wavelength dispersion'
       ;mwrfits, finaldispersion_rm[*,ifiber,iexp], fulloutname_single, hdrfloat
       ;
       ;; HDU #5 is plugmap
       ;sxaddpar, hdrplug, 'EXTNAME', 'PLUGMAP', ' Plugmap structure'
       ;mwrfits, finalplugmap_rm[ifiber,iexp], fulloutname_single, hdrplug
       ;
       ;; HDU #6 is the sky
       ;sxaddpar, hdrsky, 'EXTNAME', 'SKY', ' Subtracted sky flux'
       ;mwrfits, finalsky_rm[*,ifiber,iexp], fulloutname_single, hdrsky
       ;
     endfor
   endfor
   endfor
   endif
   return
end
;------------------------------------------------------------------------------
