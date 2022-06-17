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
;   radec_coadd    - Coadd using ra-dec matching rather then catalogID matching
;   no_reject      - Turns off rejection in the coadding of exposures
;   onestep_coadd  - Legacy algorithm for coadd. Coadding blue+red and all exposures
;                    at the the same time.
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
function unique_plmap_values, tag
   tag = STRSPLIT(tag, /EXTRACT)
   tag=tag[UNIQ(tag, SORT(tag))]
   return,STRJOIN(tag,' ')
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
 plates=plates, legacy=legacy, single_spectra=single_spectra, $
 radec_coadd=radec_coadd, no_reject=no_reject, onestep_coadd=onestep_coadd

    @specFileHdr_cards.idl


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
      radec_coadd = 1 ; no catalogid for legacy
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
   plugmap_rm=create_struct('DESIGN','','CONFIGURATION','','RA0',0.D,'DEC0',0.D,$
                            'TAI',0.D,'MJD',0.0,'AIRMASS',0.D,'DATE','','EXPTIME',0.0,$
                            'SEEING20',0.D,'SEEING50',0.D,'SEEING80',0.D)
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
          rm_plugmap[ifile].configuration=sxpar(objhdr,'CONFID')
          rm_plugmap[ifile].DESIGN=sxpar(objhdr,'DESIGNID')
        endelse
        rm_plugmap[ifile].ra0=sxpar(objhdr,'RA')
        rm_plugmap[ifile].dec0=sxpar(objhdr,'DEC')
        rm_plugmap[ifile].tai=sxpar(objhdr,'TAI-BEG')
        rm_plugmap[ifile].mjd=sxpar(objhdr,'MJD')
        rm_plugmap[ifile].airmass=sxpar(objhdr,'AIRMASS')
        rm_plugmap[ifile].SEEING20=sxpar(objhdr,'SEEING20')
        rm_plugmap[ifile].SEEING50=sxpar(objhdr,'SEEING50')
        rm_plugmap[ifile].SEEING80=sxpar(objhdr,'SEEING80')
        mjdt=sxpar(objhdr,'TAI-BEG')/(24.D*3600.D)
        mjd2datelist,mjdt,datelist=date
        rm_plugmap[ifile].date=date
        rm_plugmap[ifile].exptime=sxpar(objhdr,'EXPTIME')
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
       hdr=hdr, adderr=adderr, reslset=tempreslset

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

      thisdesign = sxpar(hdr,'DESIGNID')
      thisconfig = sxpar(hdr,'CONFID')
      if (NOT keyword_set(designlist)) then designlist = thisdesign $
       else designlist = [designlist, thisdesign]
      if (NOT keyword_set(configlist)) then configlist = thisconfig $
       else configlist = [configlist, thisconfig]

      ;----------
      ; Solve for wavelength and lambda-dispersion at each pixel in the image

      ;need to reset junk since the array lengths change
      junk=0
      traceset2xy, tempwset, junk, tempwave
      traceset2xy, tempdispset, junk, tempdispersion
      traceset2xy, tempreslset, junk, tempresolution
      ;----------
      ; Here is the correct conversion from pixels to log-lambda dispersion.
      ; We are converting from the dispersion in units of spFrame pixel sizes
      ; to the dispersion in units of the new rebinned pixel size, which is
      ; BINSZ in log-lambda units.
      
      ; this probably should be fixed elsewhere but limit where the fit range
      tempxmax=tempwset.xmax
      tempwset.xmax=(size(tempdispersion,/dimens))[0]-1
      correct_dlam, tempdispersion, 0, tempwset, dlam=binsz, /inverse
      ;correct_dlam, tempresolution, 0, tempwset, dlam=binsz, /inverse
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
         resolution = make_array(npixmax,nobj*nfiles,type=size(tempresl,/type))
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
      resolution[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempresolution
      pixelmask[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = temppixmask
      skyflux[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempsky
      ximg[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempximg
      superflat[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempsuperflat

      splog, prename=''
   endfor
   pixelmask_rm=pixelmask
   expnumf=expnumvec[uniq(expnumvec, sort(expnumvec))]

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
     nexp_tmp = nfiles/2
     for iexp=0, nexp_tmp - 1 do begin
         for iobj=0L,499L do begin
            ; for b1 and r1
            ifile = iexp
            waveb = wave[0:npixarr[ifile]-1,nobj*ifile+iobj]
            fluxb = flux[0:npixarr[ifile]-1,nobj*ifile+iobj]
            ivarb = fluxivar[0:npixarr[ifile]-1,nobj*ifile+iobj]
            ifile = iexp + nexp_tmp
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
   tempresolution = 0
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
   
   finalflux_rm = fltarr(nfinalpix, nfiber, nexp_tmp)
   finalivar_rm = fltarr(nfinalpix, nfiber, nexp_tmp)
   finalandmask_rm = lonarr(nfinalpix, nfiber, nexp_tmp)
   finalormask_rm = lonarr(nfinalpix, nfiber, nexp_tmp)
   finaldispersion_rm = fltarr(nfinalpix, nfiber, nexp_tmp)
   finalresolution_rm = fltarr(nfinalpix, nfiber, nexp_tmp)
   finalsky_rm = fltarr(nfinalpix, nfiber, nexp_tmp)
   finalplugmap_rm = replicate(plugmap[0], nfiber, nexp_tmp)


   if keyword_set(onestep_coadd) then begin
       finalra_rm = fltarr(nfiber, nexp_tmp)
       finaldec_rm = fltarr(nfiber, nexp_tmp)
       fiberid_rm = lonarr(nfiber,nexp_tmp)
       firstcarton_rm=strarr(nfiber, nexp_tmp)
       carton2TarPK_rm=LON64ARR(nfiber, nexp_tmp)
       Assigned_rm=lonarr(nfiber, nexp_tmp)
       on_target_rm=lonarr(nfiber, nexp_tmp)
       valid_rm=lonarr(nfiber, nexp_tmp)
       xfocal_rm=fltarr(nfiber,nexp_tmp)
       yfocal_rm=fltarr(nfiber,nexp_tmp)

       mjds_rm = lonarr(nfiber, nexp_tmp)
       mjds_rm_summ = lonarr(nfiber, nexp_tmp)
       config_rm = lonarr(nfiber, nexp_tmp)
       tai_rm = fltarr(nfiber, nexp_tmp)
       exptime_rm = fltarr(nfiber, nexp_tmp)
       moon_target_rm=strarr(nfiber, nexp_tmp)
       moon_phasef_rm=strarr(nfiber, nexp_tmp)
       snr2listG=strarr(nfiber, nexp_tmp)
       snr2listR=strarr(nfiber, nexp_tmp)
       snr2listI=strarr(nfiber, nexp_tmp)
       configs=strarr(nfiber, nexp_tmp)
       mjdlist_fib=strarr(nfiber, nexp_tmp)
       designs=strarr(nfiber, nexp_tmp)
       airmass_rm=strarr(nfiber, nexp_tmp)
       seeing20_rm=strarr(nfiber, nexp_tmp)
       seeing50_rm=strarr(nfiber, nexp_tmp)
       seeing80_rm=strarr(nfiber, nexp_tmp)
       weights_rm=strarr(nfiber, nexp_tmp)
       expid_rm = lonarr(nfiber, nexp_tmp)
    endif else begin
       mjds_rm_summ = lonarr(nfiber, nexp_tmp)
       finalra_rm = []
       finaldec_rm = []
       fiberid_rm = []
       firstcarton_rm=[]
       carton2TarPK_rm=[]
       Assigned_rm=[]
       on_target_rm=[]
       valid_rm=[]
       xfocal_rm=[]
       yfocal_rm=[]

       mjds_rm = []
       config_rm = []
       tai_rm = []
       exptime_rm = []
       moon_target_rm=[]
       moon_phasef_rm=[]
       snr2listG=strarr(nfiber, nexp_tmp)
       snr2listR=strarr(nfiber, nexp_tmp)
       snr2listI=strarr(nfiber, nexp_tmp)

       configs=[]
       mjdlist_fib=[]
       designs=[]
       airmass_rm=[]
       seeing20_rm=[]
       seeing50_rm=[]
       seeing80_rm=[]
       weights_rm=[]
       expid_rm=[]
    endelse
   ;----------
   
   
   combinedwave    = dblarr(nfinalpix,nobj*nexp_tmp)
   combinedflux    = fltarr(nfinalpix,nobj*nexp_tmp)
   combinedivar    = fltarr(nfinalpix,nobj*nexp_tmp)
   combineddisp    = fltarr(nfinalpix,nobj*nexp_tmp)
   combinedskyflux = fltarr(nfinalpix,nobj*nexp_tmp)
   combinedresl    = fltarr(nfinalpix,nobj*nexp_tmp)
   combinedandmask = lonarr(nfinalpix,nobj*nexp_tmp)
   combinedormask  = lonarr(nfinalpix,nobj*nexp_tmp)
   comb_plugmap    = []

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
   i=0
   for ifiber=0, nfiber-1 do begin
      for iexp=0, nexp_tmp - 1 do begin
       ; Find the first occurance of fiber number IFIBER+1
       indx = (where((plugmap.fiberid EQ ifiber+1) $
         AND (expnumvec EQ expnumf[iexp])));[0]
       if (indx[0] NE -1) then begin
         splog, 'Coadd red & blue exposure ',expnumf[iexp]
         splog, 'Fiber', ifiber+1, ' ', plugmap[indx[0]].objtype, $
           plugmap[indx[0]].mag, format = '(a, i5.4, a, a, f6.2, 5f6.2)'
         finalplugmap_rm[ifiber,iexp] = plugmap[indx[0]]
       endif
       if (indx[0] NE -1) then begin
        
         temppixmask = pixelmask_rm[*,indx]
         combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
           finalmask=temppixmask, indisp=dispersion[*,indx], $
           skyflux=skyflux[*,indx], $
           newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
           andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
           newsky=bestsky, inresl=resolution[*,indx], newresl=bestresolution, $
           nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
           maxiter=50, upper=3.0, lower=3.0, maxrej=1
         finalflux_rm[*,ifiber,iexp] = bestflux
         finalivar_rm[*,ifiber,iexp] = bestivar
         finalandmask_rm[*,ifiber,iexp] = bestandmask
         finalormask_rm[*,ifiber,iexp] = bestormask
         finaldispersion_rm[*,ifiber,iexp] = bestdispersion
         finalresolution_rm[*,ifiber,iexp] = bestresolution
         finalsky_rm[*,ifiber,iexp] = bestsky
         
         
         ; The Following are used as input for the target centric coadd
         combinedwave[*,i]    = finalwave
         combinedflux[*,i]    = bestflux
         combinedivar[*,i]    = bestivar
         combineddisp[*,i]    = bestdispersion
         combinedskyflux[*,i] = bestsky
         combinedresl[*,i]    = bestresolution
         combinedandmask[*,i] = bestandmask
         combinedormask[*,i]  = bestormask
         comb_plugmap         = [comb_plugmap, plugmap[indx[0]]]
         i++
         
         ; The following adds the COMBINEREJ bit to the input pixel masks
         pixelmask_rm[*,indx] = temppixmask
         pixelmask[*,indx] = temppixmask
      

        tai_t=rm_plugmap[iexp].tai+double(rm_plugmap[iexp].exptime/2.0)
        jdtemp=tai_t/(24.D*3600.D)
        jdtemp=jdtemp+2400000.5
        mphase,jdtemp,mfrac
        moonpos,jdtemp,ra_moon,dec_moon
        ratemp=plugmap[indx[0]].ra
        dectemp=plugmap[indx[0]].dec
        moon_dist = djs_diff_angle(ra_moon, dec_moon, ratemp, dectemp)

         if keyword_set(onestep_coadd) then begin
            finalra_rm[ifiber,iexp]=ratemp
            finaldec_rm[ifiber,iexp]=dectemp
            mjds_rm[ifiber,iexp]=rm_plugmap[iexp].mjd
            mjds_rm_summ[ifiber,iexp]=rm_plugmap[iexp].mjd
            fiberid_rm[ifiber, iexp] = plugmap[indx[0]].fiberid
            firstcarton_rm[ifiber, iexp] = plugmap[indx[0]].firstcarton
            carton2TarPK_rm[ifiber, iexp] = plugmap[indx[0]].carton_to_target_pk
            Assigned_rm[ifiber, iexp] = plugmap[indx[0]].assigned
            on_target_rm[ifiber, iexp] = plugmap[indx[0]].on_target
            valid_rm[ifiber, iexp] = plugmap[indx[0]].valid
            xfocal_rm[ifiber, iexp] = plugmap[indx[0]].xfocal
            yfocal_rm[ifiber, iexp] = plugmap[indx[0]].yfocal
            tai_rm[ifiber,iexp]=rm_plugmap[iexp].tai+double(rm_plugmap[iexp].exptime/2.0)
            exptime_rm[ifiber,iexp]=rm_plugmap[iexp].exptime
            ; use expuse number instad of configuration number for legacy
            config_rm[ifiber,iexp]=rm_plugmap[iexp].configuration
            designs[ifiber,iexp] =rm_plugmap[iexp].design
            airmass_rm[ifiber,iexp] = rm_plugmap[iexp].airmass
            seeing20_rm[ifiber,iexp] =rm_plugmap[iexp].SEEING20
            seeing50_rm[ifiber,iexp] =rm_plugmap[iexp].SEEING50
            seeing80_rm[ifiber,iexp] =rm_plugmap[iexp].SEEING80
            expid_rm[ifiber,iexp] = iexp
            moon_target_rm[ifiber,iexp] =moon_dist
            moon_phasef_rm[ifiber,iexp] =mfrac
         endif else begin
            finalra_rm=[finalra_rm,ratemp]
            finaldec_rm=[finaldec_rm,dectemp]
            mjds_rm=[mjds_rm,rm_plugmap[iexp].mjd]
            mjds_rm_summ[ifiber,iexp]=rm_plugmap[iexp].mjd
            fiberid_rm = [fiberid_rm,plugmap[indx[0]].fiberid]
            firstcarton_rm = [firstcarton_rm,plugmap[indx[0]].firstcarton]
            carton2TarPK_rm = [carton2TarPK_rm,plugmap[indx[0]].carton_to_target_pk]
            Assigned_rm = [Assigned_rm,plugmap[indx[0]].assigned]
            on_target_rm = [on_target_rm,plugmap[indx[0]].on_target]
            valid_rm = [valid_rm,plugmap[indx[0]].valid]
            xfocal_rm = [xfocal_rm,plugmap[indx[0]].xfocal]
            yfocal_rm = [yfocal_rm,plugmap[indx[0]].yfocal]
            tai_rm = [tai_rm,rm_plugmap[iexp].tai+double(rm_plugmap[iexp].exptime/2.0)]
            exptime_rm = [exptime_rm,rm_plugmap[iexp].exptime]
            ; use expuse number instad of configuration number for legacy
            config_rm = [config_rm,rm_plugmap[iexp].configuration]
            designs = [designs,rm_plugmap[iexp].design]
            airmass_rm=[airmass_rm,rm_plugmap[iexp].airmass]
            seeing20_rm=[seeing20_rm,rm_plugmap[iexp].SEEING20]
            seeing50_rm=[seeing50_rm,rm_plugmap[iexp].SEEING50]
            seeing80_rm=[seeing80_rm,rm_plugmap[iexp].SEEING80]
            expid_rm=[expid_rm,iexp]
            
            moon_target_rm=[moon_target_rm,moon_dist]
            moon_phasef_rm=[moon_phasef_rm,mfrac]
         endelse
       endif else begin
         splog, 'Fiber', ifiber+1, ' NO DATA'
         finalandmask_rm[*,ifiber,iexp] = pixelmask_bits('NODATA')
         finalormask_rm[*,ifiber,iexp] = pixelmask_bits('NODATA')
       endelse
      endfor
   endfor
   
   
   
   ;---------------------------------------------------------------------------
   ; FLUX DISTORTION IMAGE
   ;---------------------------------------------------------------------------
   bighdr = *hdrarr[0] ;1st file's header to use for the combined plate header

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

      endfor
   endif
   
   if not keyword_set(onestep_coadd) then begin
        in_plugmap=plugmap
        plugmap=comb_plugmap
   endif
   if keyword_set(radec_coadd) then begin
       if (not keyword_set(legacy)) and (not keyword_set(plates)) then $
           splog, 'WARNING: Coadds May include mix of valid/on_target flags'
       ;Set a list of targets from its coordinates, this block code considers
       ;the posibility to observe the same target at a diferent fiber in a
       ;diferent FPS configuartion
       brake=0
       indx0=0
       ra_rm=plugmap.ra
       dec_rm=plugmap.dec
       ra_tp=ra_rm
       dec_tp=dec_rm
       while brake eq 0 do begin
         indx1=where((ra_rm eq ra_tp[0]) and (dec_rm eq dec_tp[0]))
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
             indx1=where((ra_rm eq ra_tp[0]) and (dec_rm eq dec_tp[0]))
             indx1=indx1[0]
           endif else begin
             brake=1
           endelse
           indx0+=1
         endif
       endwhile
   endif else begin
       ;Set a list of targets from its catalogid, this block code considers
       ;the posibility to observe the same target at a diferent fiber in a
       ;diferent FPS configuartion, but will have the same catalogid
       ;(within a catalog cross match - TODO: coadd across catalog cross matchs)
       catid_rm = plugmap.catalogid
       catid_rm = catid_rm[UNIQ(catid_rm, SORT(catid_rm))]
       undefine, indx_tar
       foreach catid, catid_rm do begin
         if keyword_set(plates) then begin
            indx1=where(plugmap.catalogid eq catid, ct)
         endif else begin
            indx1=where((plugmap.catalogid eq catid) AND (plugmap.Assigned eq 1) $
                         AND (plugmap.valid eq 1) AND (plugmap.On_target eq 1), ct)
            ; below is primarily to produce outputs of "unplugged" targets
            ; so we can include extra skys fibers for outputs
            if ct eq 0 then indx1 = where((plugmap.catalogid eq catid), ct)
         endelse
         if ct ne 0 then begin
            if not keyword_set(indx_tar) then indx_tar=[indx1] $
            else indx_tar=[indx_tar,'-10',indx1]
         endif
       endforeach
   endelse
   
   indx_tar=['-10',indx_tar,'-10']
   nt=where(indx_tar eq -10,ntarget)
   ntarget=ntarget-1
   finalflux = fltarr(nfinalpix, ntarget)
   finalivar = fltarr(nfinalpix, ntarget)
   finalandmask = lonarr(nfinalpix, ntarget)
   finalormask = lonarr(nfinalpix, ntarget)
   finaldispersion = fltarr(nfinalpix, ntarget)
   finalresolution = fltarr(nfinalpix, ntarget)
   finalsky = fltarr(nfinalpix, ntarget)
   finalplugmap = replicate(plugmap[0], ntarget)
   mjds = lonarr(ntarget)
   final_ra = dblarr(ntarget)
   final_dec = dblarr(ntarget)
   indx_target=intarr(ntarget)
   nexp_target=intarr(ntarget)
   mjdsfinal = dblarr(ntarget)
   mjd_t = dblarr(ntarget)
   snr_t = dblarr(ntarget)
   exp_disp_med = dblarr(ntarget)
   weight = dblarr(nexp_tmp)
   snr2G_target = strarr(ntarget)
   snr2R_target = strarr(ntarget)
   snr2I_target = strarr(ntarget)

   mjdlist_target = strarr(ntarget)
   designs_target = strarr(ntarget)
   configs_target = strarr(ntarget)

   moon_target = strarr(ntarget)
   moon_phasef = strarr(ntarget)
   airmass_target = strarr(ntarget)
   seeing20_target = strarr(ntarget)
   seeing50_target = strarr(ntarget)
   seeing80_target = strarr(ntarget)
   tai_target = strarr(ntarget)
   fiber_target = strarr(ntarget)
   RA_target = strarr(ntarget)
   DEC_target = strarr(ntarget)
   Firstcarton_target = strarr(ntarget)
   carton2TarPK_target = strarr(ntarget)
   Assigned_target = strarr(ntarget)
   on_target_target = strarr(ntarget)
   valid_target = strarr(ntarget)
   xfocal_target = strarr(ntarget)
   yfocal_target = strarr(ntarget)
   exptime_target = dblarr(ntarget)
   weights_target = strarr(ntarget)

   fiber_target_s=replicate(create_struct('FIBERID_LIST',' '),ntarget)
   RA_target_s=replicate(create_struct('RA_LIST',' '),ntarget)
   DEC_target_s=replicate(create_struct('DEC_LIST',' '),ntarget)
   FIRSTCARTON_target_s=replicate(create_struct('FIRSTCARTON_LIST',' '),ntarget)
   cartoon2TarPK_target_s=replicate(create_struct('CARTON_TO_TARGET_PK_LIST',' '),ntarget)
   Assigned_target_s=replicate(create_struct('ASSIGNED_LIST',' '),ntarget)
   on_target_target_s=replicate(create_struct('ON_TARGET_LIST',' '),ntarget)
   valid_target_s=replicate(create_struct('VALID_LIST',' '),ntarget)
   exp_disp_med_s=replicate(create_struct('EXP_DISP_MED',0.D),ntarget)
   xfocal_target_s=replicate(create_struct('XFOCAL_LIST',' '),ntarget)
   yfocal_target_s=replicate(create_struct('YFOCAL_LIST',' '),ntarget)
   indx_target_s=replicate(create_struct('target_index',0),ntarget)
   nexp_target_s=replicate(create_struct('nexp',0),ntarget)
   exptime_target_s=replicate(create_struct('exptime',0),ntarget)
   mjdf_target_s=replicate(create_struct('MJD_FINAL',0.D),ntarget)
   moon_target_s=replicate(create_struct('MOON_DIST',' '),ntarget)
   moon_phasef_s=replicate(create_struct('MOON_PHASE',' '),ntarget)
   airmass_s=replicate(create_struct('AIRMASS',0.D, 'AIRMASS_LIST',' '),ntarget)
   seeing20_s=replicate(create_struct('SEEING20',0.D,'SEEING20_LIST',' '),ntarget)
   seeing50_s=replicate(create_struct('SEEING50',0.D,'SEEING50_LIST',' '),ntarget)
   seeing80_s=replicate(create_struct('SEEING80',0.D,'SEEING80_LIST',' '),ntarget)
   tai_target_s=replicate(create_struct('TAI_LIST',' '),ntarget)
   snr2G_target_s=replicate(create_struct('FIELDSNR2G_LIST',' '),ntarget)
   snr2R_target_s=replicate(create_struct('FIELDSNR2R_LIST',' '),ntarget)
   snr2I_target_s=replicate(create_struct('FIELDSNR2I_LIST',' '),ntarget)
   mjdlist_target_s=replicate(create_struct('MJDLIST',' '),ntarget)
   designs_target_s=replicate(create_struct('DESIGNS',' '),ntarget)
   configs_target_s=replicate(create_struct('CONFIGS',' '),ntarget)
   struct_assign, {fiberid: 0L}, finalplugmap ; Zero out all elements in this
   ; FINALPLUGMAP structure.
   for itarget=0, ntarget-1 do begin
      indx=indx_tar[nt[itarget]+1:nt[itarget+1]-1]
      if (indx[0] NE -1) then begin
         if keyword_set(radec_coadd) then $
           splog, 'Coadd all the exposures with the same coordinates ('+strtrim(itarget+1,2)+'/'+strtrim(ntarget,2)+')' $
         else splog, 'Coadd all the exposures with the same CatalogID ('+strtrim(itarget+1,2)+'/'+strtrim(ntarget,2)+')'
         splog, 'Target', itarget+1, ' ', plugmap[indx[0]].objtype, $
          plugmap[indx[0]].mag, format = '(a, i5.4, a, a, f6.2, 5f6.2)'
         finalplugmap[itarget] = plugmap[indx[0]]
         mjds[itarget]=mjds_rm[indx[0]]
         final_ra[itarget]=plugmap[indx[0]].ra
         final_dec[itarget]=plugmap[indx[0]].dec


         if keyword_set(no_reject) then begin
           bestandmask= combinedandmask[*,indx]
           bestormask = combinedormask[*,indx]
           temppixmask = combinedandmask[*,indx]
           rm_combine1fiber, combinedwave[*,indx], combinedflux[*,indx], $
                     combinedivar[*,indx], finalmask=temppixmask, $
                     indisp=combineddisp[*,indx], skyflux=combinedskyflux[*,indx],$
                     inormask=combinedormask[*,indx], inandmask=combinedandmask[*,indx],$
                     inresl=combinedresl[*,indx], $
                     andmask=bestandmask, ormask=bestormask, $
                     newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
                     newdisp=bestdispersion, newsky=bestsky, newresl=bestresolution, $
                     nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
                     maxiter=0, upper=3d6, lower=3d6, maxrej=1
                     
                     
           ; only calculate and report median exposure dispersion if valid target
           if (not keyword_set(legacy)) and (not keyword_set(plates)) then begin
            test=where((plugmap[indx].Assigned eq 1) $
                         AND (plugmap[indx].valid eq 1) $
                         AND (plugmap[indx].On_target eq 1), ctv)
           endif else ctv=1

           if ctv ne 0 then begin
             if n_elements(indx) gt 1 then begin
               exp_disp_med[itarget] = abs(median(STDDEV(combinedflux,DIMENSION=2,/double)/bestflux,$
                                              /EVEN,/double))
             endif else exp_disp_med[itarget] = 0.0
           endif else exp_disp_med[itarget]= -1.d
         endif else begin
           if not keyword_set(onestep_coadd) then begin
             bestandmask= combinedandmask[*,indx]
             bestormask = combinedormask[*,indx]
             temppixmask = combinedandmask[*,indx]
             rm_combine1fiber, combinedwave[*,indx], combinedflux[*,indx], $
                     combinedivar[*,indx], finalmask=temppixmask, $
                     indisp=combineddisp[*,indx], skyflux=combinedskyflux[*,indx],$
                     inormask=combinedormask[*,indx], inandmask=combinedandmask[*,indx],$
                     inresl=combinedresl[*,indx], $
                     andmask=bestandmask, ormask=bestormask, $
                     newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
                     newdisp=bestdispersion, newsky=bestsky, newresl=bestresolution, $
                     nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
                     maxiter=50, upper=3.0, lower=3.0, maxrej=1
                     
             ; only calculate and report median exposure dispersion if valid target
             if (not keyword_set(legacy)) and (not keyword_set(plates)) then begin
               test=where((plugmap[indx].Assigned eq 1) $
                           AND (plugmap[indx].valid eq 1) $
                           AND (plugmap[indx].On_target eq 1), ctv)
             endif else ctv=1

             if ctv ne 0 then begin
               if n_elements(indx) gt 1 then begin
                exp_disp_med[itarget] = abs(median(STDDEV(combinedflux,DIMENSION=2,/double)/bestflux,$
                                                /EVEN,/double))
               endif else exp_disp_med[itarget] = 0.0
             endif else exp_disp_med[itarget]= -1.d

           endif else begin
             ; This is for legacy to match older coadd algorithm
             temppixmask = pixelmask[*,indx]
             combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
                     finalmask=temppixmask, indisp=dispersion[*,indx], $
                     skyflux=skyflux[*,indx], $
                     newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
                     andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
                     newsky=bestsky, inresl=resolution[*,indx], newresl=bestresolution, $
                     nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
                     maxiter=50, upper=3.0, lower=3.0, maxrej=1
             exp_disp_med[itarget]= -1.
             ; The following adds the COMBINEREJ bit to the input pixel masks
             pixelmask[*,indx] = temppixmask
           endelse
         endelse
         finalflux[*,itarget] = bestflux
         finalivar[*,itarget] = bestivar
         finalandmask[*,itarget] = bestandmask
         finalormask[*,itarget] = bestormask
         finaldispersion[*,itarget] = bestdispersion
         finalresolution[*,itarget] = bestresolution
         finalsky[*,itarget] = bestsky
         
         indx_target[itarget]=itarget+1
         nexp_target[itarget]=n_elements(indx)
      endif else begin
         splog, 'Target', itarget+1, ' NO DATA'
         finalandmask[*,itarget] = pixelmask_bits('NODATA')
         finalormask[*,itarget] = pixelmask_bits('NODATA')
         indx_target[itarget]=itarget
         nexp_target[itarget]=0
         exp_disp_med[itarget]= -1.
      endelse
   endfor
   indx_target_s.target_index=indx_target
   nexp_target_s.nexp=nexp_target
   finalplugmap=struct_addtags(finalplugmap,indx_target_s)
   finalplugmap=struct_addtags(finalplugmap,nexp_target_s)

   
   ;---------------------------------------------------------------------------
   ; FLUX DISTORTION IMAGE
   ;---------------------------------------------------------------------------
   if keyword_set(onestep_coadd) then begin
     ; compute the Flux distortion of the coadd of exposures
      splog, 'Compute the flux distortion image for all exposures'
      corrimg = flux_distortion(finalflux, finalivar, finalandmask, finalormask, $
       plugmap=finalplugmap, loglam=finalwave, plotfile=distortpsfile, hdr=bighdr, $
       legacy=legacy)
      igood = where(finalivar GT 0)
      thismin = min(corrimg[igood], max=thismax)
      cratio = thismin / thismax
      if (cratio LT 1./100) then begin
         splog, 'WARNING: Flux distortion image dynamic range = ', 1./cratio, ' (DISABLE)'
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
   endif

   ;----------
   ; Modify the 1st file's header to use for the combined plate header.


   tai_flag=0
   ;mjd_t=0.0
   ;snr_t=0.0
   
   master_snr2=dblarr(2,3)
   master_snr2_dered=dblarr(2,3)



   for iexp=0, nexp_tmp - 1 do begin
        platesn, finalflux_rm[*,*,iexp], finalivar_rm[*,*,iexp], $
          finalandmask_rm[*,*,iexp], finalplugmap_rm[*,iexp], finalwave, hdr=bighdr, $
          legacy=legacy, plotfile=djs_filepath(repstr(plotsnfile,'X',string(iexp,format='(i2.2)')), root_dir=combinedir), $
          coeffs=coeffs, snplate=snplate, specsnlimit=specsnlimit, dered_snplate=dered_snplate
        splog, prelog=''
        bands = ['G','R','I']
        if keyword_set(legacy) then sp_n=2 else sp_n=1
        for ispec=1, sp_n do begin
           for bb=0, n_elements(bands)-1 do begin
              key1 = 'SNC0'+ strtrim(ispec,2)+strupcase(bands[bb])+string(iexp,format='(i2.2)')
              comment = string(format='(a,i2,a,i2.2,a)', $
               ' SN fit coeff for spec', ispec, ', exp ', iexp ,' at '+strupcase(bands[bb])+'-band')
              sxaddpar, bighdr, key1, coeffs[bb,(ispec-1)*2], comment, before='LOWREJ'
              key1 = 'SNC1'+strtrim(ispec,2)+strupcase(bands[bb])+string(iexp,format='(i2.2)')
              comment = string(format='(a,i2,a,i2.2,a)', $
               ' SN fit coeff for spec', ispec, ', exp ', iexp ,' at '+strupcase(bands[bb])+'-band')
              sxaddpar, bighdr, key1, coeffs[bb,(ispec-1)*2+1], $
               comment, before='LOWREJ'
              key1 = 'SN2_'+strtrim(ispec,2)+strupcase(bands[bb])+string(iexp,format='(i2.2)')
              comment = string(format='(a,i2,a,i2.2,a,f5.2,a)', $
               ' (S/N)^2 for spec ', ispec, ', exp ', iexp ,' at mag ', specsnlimit[bb].snmag,' at '+strupcase(bands[bb])+'-band' )
              sxaddpar, bighdr, key1, snplate[ispec-1,bb], comment, before='LOWREJ'
           endfor
        endfor

        master_snr2=master_snr2+snplate
        master_snr2_dered=master_snr2_dered+dered_snplate


;        snr_t=snplate[0,2]+snr_t;use SNR2 in i-band
;        tai_t=rm_plugmap[iexp].tai+rm_plugmap[iexp].exptime/2.0
;        mjd_t=tai_t/(24.D*3600.D)*snplate[0,2]+mjd_t
;        jdtemp=tai_t/(24.D*3600.D)
;        jdtemp=jdtemp+2400000.5
;        mphase,jdtemp,mfrac
;        moonpos,jdtemp,ra_moon,dec_moon
;        ra_t=finalplugmap_rm[*,iexp].ra
;        dec_t=finalplugmap_rm[*,iexp].dec
        for ifib=0, nfiber-1 do begin
;            moon_dist = djs_diff_angle(ra_moon, dec_moon, ra_t[ifib], dec_t[ifib])
;            if keyword_set(onestep_coadd) then begin
;                moon_target_rm[ifib,iexp]=strtrim(strcompress(string(string(moon_dist,format='(f0.1)'),format='(999a)')),2)
;                moon_phasef_rm[ifib,iexp]=strtrim(strcompress(string(string(mfrac,format='(f0.2)'),format='(999a)')),2)
;                airmass_rm[ifib,iexp]=strtrim(strcompress(string(rm_plugmap[iexp].airmass,format='(999a)')),2)
;                seeing20_rm[ifib,iexp]=strtrim(strcompress(string(rm_plugmap[iexp].SEEING20,format='(999a)')),2)
;                seeing50_rm[ifib,iexp]=strtrim(strcompress(string(rm_plugmap[iexp].SEEING50,format='(999a)')),2)
;                seeing80_rm[ifib,iexp]=strtrim(strcompress(string(rm_plugmap[iexp].SEEING80,format='(999a)')),2)
                snr2listG[ifib,iexp]=strtrim(strcompress(string(string(snplate[0,0],format='(f0.2)'),format='(999a)')),2)
                snr2listR[ifib,iexp]=strtrim(strcompress(string(string(snplate[0,1],format='(f0.2)'),format='(999a)')),2)
                snr2listI[ifib,iexp]=strtrim(strcompress(string(string(snplate[0,2],format='(f0.2)'),format='(999a)')),2)
;                mjdlist_fib[ifib,iexp]=strtrim(strcompress(string(string(rm_plugmap[iexp].mjd,format='(i5)'),format='(999a)')),2)
;                configs[ifib,iexp]=strtrim(strcompress(string(rm_plugmap[iexp].configuration,format='(999a)')),2)
;                designs[ifib,iexp]=strtrim(strcompress(string(rm_plugmap[iexp].DESIGN,format='(999a)')),2)
;                weights_rm[ifib,iexp] =strtrim(strcompress(string(snplate[0,2],format='(999a)')),2)
;            endif else begin
;                moon_target_rm=[moon_target_rm,strtrim(strcompress(string(string(moon_dist,format='(f0.1)'),format='(999a)')),2)]
;                moon_phasef_rm=[moon_phasef_rm,strtrim(strcompress(string(string(mfrac,format='(f0.2)'),format='(999a)')),2)]
;                airmass_rm=[airmass_rm,strtrim(strcompress(string(rm_plugmap[iexp].airmass,format='(999a)')),2)]
;                seeing20_rm=[seeing20_rm,strtrim(strcompress(string(rm_plugmap[iexp].SEEING20,format='(999a)')),2)]
;                seeing50_rm=[seeing50_rm,strtrim(strcompress(string(rm_plugmap[iexp].SEEING50,format='(999a)')),2)]
;                seeing80_rm=[seeing80_rm,strtrim(strcompress(string(rm_plugmap[iexp].SEEING80,format='(999a)')),2)]
;                snr2listG=[snr2listG,strtrim(strcompress(string(string(snplate[0,0],format='(f0.2)'),format='(999a)')),2)]
;                snr2listR=[snr2listR,strtrim(strcompress(string(string(snplate[0,1],format='(f0.2)'),format='(999a)')),2)]
;                snr2listI=[snr2listI,strtrim(strcompress(string(string(snplate[0,2],format='(f0.2)'),format='(999a)')),2)]
;                mjdlist_fib=[mjdlist_fib,strtrim(strcompress(string(string(rm_plugmap[iexp].mjd,format='(i5)'),format='(999a)')),2)]
;                configs=[configs,strtrim(strcompress(string(rm_plugmap[iexp].configuration,format='(999a)')),2)]
;                designs=[designs,strtrim(strcompress(string(rm_plugmap[iexp].DESIGN,format='(999a)')),2)]
;                weights_rm =[weights_rm,strtrim(strcompress(string(snplate[0,2],format='(999a)')),2)]
;            endelse
        endfor
        if (NOT keyword_set(tailist)) then tailist = tai_t $
        else tailist = [tailist, tai_t]
        
        tai_flag=1
        weight[iexp]=snplate[0,2]
   endfor

   ; Plot S/N and throughput **after** this distortion-correction.
   ; (This over-writes header cards written in the first call.)
   splog, prelog='Final'
   platesn, finalflux, finalivar, finalandmask, finalplugmap, finalwave, $
   hdr=bighdr, legacy=legacy, plotfile=djs_filepath(repstr(plotsnfile,'-X',''), root_dir=combinedir), coeffs=coeffs
   splog, prelog=''
   bands = ['G','R','I']
   if keyword_set(legacy) then sp_n=2 else sp_n=1
   for ispec=1, sp_n do begin
        for bb=0, n_elements(bands)-1 do begin
            key1 = 'SNC0'+ strtrim(ispec,2)+strupcase(bands[bb])
            comment = string(format='(a,i2,a)', $
            ' Total SN fit coeff for spec', ispec,' at '+strupcase(bands[bb])+'-band' )
            sxaddpar, bighdr, key1, coeffs[bb,(ispec-1)*2], comment, before='LOWREJ'
            key1 = 'SNC1'+strtrim(ispec,2)+strupcase(bands[bb])
            comment = string(format='(a,i2,a)', $
            ' Total SN fit coeff for spec', ispec,' at '+strupcase(bands[bb])+'-band')
            sxaddpar, bighdr, key1, coeffs[bb,(ispec-1)*2+1], $
            comment, before='LOWREJ'
        endfor
   endfor

   bands = ['G','R','I']
   if keyword_set(legacy) then sp_n=2 else sp_n=1
   for ispec=1, sp_n do begin
       for bb=0, n_elements(bands)-1 do begin

           ; Standard (S/N)^2
           key1 = 'SPEC'+ strtrim(ispec,2)+'_'+strupcase(bands[bb])
           snr2= sxpar(bighdr, key1)
           comment = string(format='(a,i2,a,f5.2)', $
               ' (S/N)^2 for spec ', ispec, ' at mag ', specsnlimit[bb].snmag)
           sxaddpar, bighdr, key1, master_snr2[ispec-1,bb], comment;, before='NSTD'
           key2 = 'FSPEC'+ strtrim(ispec,2)+'_'+strupcase(bands[bb])
           comment = string(format='(a,i2,a,f5.2)', $
               'Fit (S/N)^2 for spec ', ispec, ' at mag ', specsnlimit[bb].snmag)
           sxaddpar, bighdr, key2, snr2, comment, after=key1

           ; Extinction corrected (S/N)^2
           key1 = 'SN2EXT'+strtrim(ispec,2)+strupcase(bands[bb])
           snr2= sxpar(bighdr, key1)
           comment = ' Extinction corrected (S/N)^2'
           sxaddpar, bighdr, key1, master_snr2_dered[ispec-1,bb], $
                comment;, before='NSTD'
           key2 = 'FSN2EX'+ strtrim(ispec,2)+strupcase(bands[bb])
           comment = ' Extinction corrected Fit (S/N)^2'
           sxaddpar, bighdr, key2, snr2, comment, after=key1
       endfor
   endfor




   for itarget=0, ntarget-1 do begin
         indx=indx_tar[nt[itarget]+1:nt[itarget+1]-1]
         if (indx[0] NE -1) then begin
            moon_target[itarget]=strtrim(strcompress(string(string(moon_target_rm[indx[0]],format='(f0.1)'),format='(999a)')),2)
            moon_phasef[itarget]=strtrim(strcompress(string(string(moon_phasef_rm[indx[0]],format='(f0.2)'),format='(999a)')),2)
            fiber_target[itarget]=strtrim(strcompress(string(fiberid_rm[indx[0]],format='(999a)')),2)
            RA_target[itarget]=strtrim(strcompress(string(string(finalra_rm[indx[0]],format='(f0.6)'),format='(999a)')),2)
            DEC_target[itarget]=strtrim(strcompress(string(string(finaldec_rm[indx[0]],format='(f0.6)'),format='(999a)')),2)
            Firstcarton_target[itarget]=strtrim(strcompress(string(firstcarton_rm[indx[0]],format='(999a)')),2)
            carton2TarPK_target[itarget]=strtrim(strcompress(string(carton2TarPK_rm[indx[0]],format='(999a)')),2)
            Assigned_target[itarget]=strtrim(strcompress(string(Assigned_rm[indx[0]],format='(999a)')),2)
            on_target_target[itarget]=strtrim(strcompress(string(string(on_target_rm[indx[0]],format='(i15)'),format='(999a)')),2)
            valid_target[itarget]=strtrim(strcompress(string(valid_rm[indx[0]],format='(999a)')),2)
            xfocal_target[itarget]=strtrim(strcompress(string(string(xfocal_rm[indx[0]],format='(f0.3)'),format='(999a)')),2)
            yfocal_target[itarget]=strtrim(strcompress(string(string(yfocal_rm[indx[0]],format='(f0.3)'),format='(999a)')),2)
            tai_target[itarget]=strtrim(strcompress(string(string(tai_rm[indx[0]],format='(i15)'),format='(999a)')),2)
            tai_t=tai_rm[indx[0]]
            if nexp_tmp eq 1 then begin
                match = [fiberid_rm[indx[0]]-1]
                mjd_t[itarget]=tai_t/(24.D*3600.D)*snr2listI[match]
                snr_t[itarget]=snr2listI[match] ;use SNR2 in i-band
                snr2G_target[itarget]=strtrim(strcompress(string(string(snr2listG[match],format='(f0.2)'),format='(999a)')),2)
                snr2R_target[itarget]=strtrim(strcompress(string(string(snr2listR[match],format='(f0.2)'),format='(999a)')),2)
                snr2I_target[itarget]=strtrim(strcompress(string(string(snr2listI[match],format='(f0.2)'),format='(999a)')),2)
            endif else begin
                match=[fiberid_rm[indx[0]]-1,expid_rm[indx[0]]]
                mjd_t[itarget]=tai_t/(24.D*3600.D)*snr2listI[match[0],match[1]]
                snr_t[itarget]=snr2listI[match[0],match[1]] ;use SNR2 in i-band
                snr2G_target[itarget]=strtrim(strcompress(string(string(snr2listG[match[0],match[1]],format='(f0.2)'),format='(999a)')),2)
                snr2R_target[itarget]=strtrim(strcompress(string(string(snr2listR[match[0],match[1]],format='(f0.2)'),format='(999a)')),2)
                snr2I_target[itarget]=strtrim(strcompress(string(string(snr2listI[match[0],match[1]],format='(f0.2)'),format='(999a)')),2)
                weights_target[itarget]=strtrim(strcompress(string(string(snr2listI[match[0],match[1]],format='(f0.10)'),format='(999a)')),2)
            endelse
            mjdlist_target[itarget]=strtrim(strcompress(string(string(mjds_rm[indx[0]],format='(i5)'),format='(999a)')),2)
            designs_target[itarget]=strtrim(strcompress(string(designs[indx[0]],format='(999a)')),2)
            configs_target[itarget]=strtrim(strcompress(string(config_rm[indx[0]],format='(999a)')),2)
            airmass_target[itarget]=strtrim(strcompress(string(airmass_rm[indx[0]],format='(999a)')),2)
            seeing20_target[itarget]=strtrim(strcompress(string(seeing20_rm[indx[0]],format='(999a)')),2)
            seeing50_target[itarget]=strtrim(strcompress(string(seeing50_rm[indx[0]],format='(999a)')),2)
            seeing80_target[itarget]=strtrim(strcompress(string(seeing80_rm[indx[0]],format='(999a)')),2)
            exptime_target[itarget]=exptime_rm[indx[0]]
            if n_elements(indx) gt 1 then begin
               if not keyword_set(onestep_coadd) then nindx=n_elements(indx) else nindx=n_elements(indx)/2
               for iexp=1, nindx-1 do begin
                  moon_target[itarget]=moon_target[itarget]+' '+strtrim(strcompress(string(string(moon_target_rm[indx[iexp]],format='(f0.1)'),format='(999a)')),2)
                  moon_phasef[itarget]=moon_phasef[itarget]+' '+strtrim(strcompress(string(string(moon_phasef_rm[indx[iexp]],format='(f0.2)'),format='(999a)')),2)
                  fiber_target[itarget]=fiber_target[itarget]+' '+strtrim(strcompress(string(fiberid_rm[indx[iexp]],format='(999a)')),2)
                  RA_target[itarget]=RA_target[itarget]+' '+strtrim(strcompress(string(string(finalra_rm[indx[iexp]],format='(f0.6)'),format='(999a)')),2)
                  DEC_target[itarget]=DEC_target[itarget]+' '+strtrim(strcompress(string(string(finaldec_rm[indx[iexp]],format='(f0.6)'),format='(999a)')),2)
                  Firstcarton_target[itarget]=Firstcarton_target[itarget]+' '+strtrim(strcompress(string(firstcarton_rm[indx[iexp]],format='(999a)')),2)
                  carton2TarPK_target[itarget]=carton2TarPK_target[itarget]+' '+strtrim(strcompress(string(carton2TarPK_rm[indx[iexp]],format='(999a)')),2)
                  Assigned_target[itarget]=Assigned_target[itarget]+' '+strtrim(strcompress(string(string(Assigned_rm[indx[iexp]], format='(i15)'),format='(999a)')),2)
                  on_target_target[itarget]=on_target_target[itarget]+' '+strtrim(strcompress(string(string(on_target_rm[indx[iexp]], format='(i15)'),format='(999a)')),2)
                  valid_target[itarget]=valid_target[itarget]+' '+strtrim(strcompress(string(string(valid_rm[indx[iexp]], format='(i15)'),format='(999a)')),2)
                  xfocal_target[itarget]=xfocal_target[itarget]+' '+strtrim(strcompress(string(string(xfocal_rm[indx[iexp]],format='(f0.3)'),format='(999a)')),2)
                  yfocal_target[itarget]=yfocal_target[itarget]+' '+strtrim(strcompress(string(string(yfocal_rm[indx[iexp]],format='(f0.3)'),format='(999a)')),2)
                  tai_target[itarget]=tai_target[itarget]+' '+strtrim(strcompress(string(string(tai_rm[indx[iexp]],format='(i15)'),format='(999a)')),2)
                  tai_t=tai_rm[indx[iexp]]
                  if nexp_tmp eq 1 then begin
                      match = [fiberid_rm[indx[iexp]]-1]
                      mjd_t[itarget]=mjd_t[itarget]+tai_t/(24.D*3600.D)*snr2listI[match]
                      snr_t[itarget]=snr_t[itarget]+snr2listI[match] ;use SNR2 in i-band
                      snr2G_target[itarget]=snr2G_target[itarget]+' '+strtrim(strcompress(string(string(snr2listG[match],format='(f0.2)'),format='(999a)')),2)
                      snr2R_target[itarget]=snr2R_target[itarget]+' '+strtrim(strcompress(string(string(snr2listR[match],format='(f0.2)'),format='(999a)')),2)
                      snr2I_target[itarget]=snr2I_target[itarget]+' '+strtrim(strcompress(string(string(snr2listI[match],format='(f0.2)'),format='(999a)')),2)
                      weights_target[itarget]=weights_target[itarget]+' '+strtrim(strcompress(string(string(snr2listI[match],format='(f0.10)'),format='(999a)')),2)
                  endif else begin
                      match=[fiberid_rm[indx[iexp]]-1,expid_rm[indx[iexp]]]
                      mjd_t[itarget]=mjd_t[itarget]+tai_t/(24.D*3600.D)*snr2listI[match[0],match[1]]
                      snr_t[itarget]=snr_t[itarget]+snr2listI[match[0],match[1]] ;use SNR2 in i-band
                      snr2G_target[itarget]=snr2G_target[itarget]+' '+strtrim(strcompress(string(string(snr2listG[match[0],match[1]],format='(f0.2)'),format='(999a)')),2)
                      snr2R_target[itarget]=snr2R_target[itarget]+' '+strtrim(strcompress(string(string(snr2listR[match[0],match[1]],format='(f0.2)'),format='(999a)')),2)
                      snr2I_target[itarget]=snr2I_target[itarget]+' '+strtrim(strcompress(string(string(snr2listI[match[0],match[1]],format='(f0.2)'),format='(999a)')),2)
                      weights_target[itarget]=weights_target[itarget]+' '+strtrim(strcompress(string(string(snr2listI[match[0],match[1]],format='(f0.10)'),format='(999a)')),2)
                  endelse
                  mjdlist_target[itarget]=mjdlist_target[itarget]+' '+strtrim(strcompress(string(string(mjds_rm[indx[0]],format='(i5)'),format='(999a)')),2)
                  designs_target[itarget]=designs_target[itarget]+' '+strtrim(strcompress(string(designs[indx[iexp]],format='(999a)')),2)
                  configs_target[itarget]=configs_target[itarget]+' '+strtrim(strcompress(string(config_rm[indx[iexp]],format='(999a)')),2)
                  airmass_target[itarget]=airmass_target[itarget]+' '+strtrim(strcompress(string(airmass_rm[indx[iexp]],format='(999a)')),2)
                  seeing20_target[itarget]=seeing20_target[itarget]+' '+strtrim(strcompress(string(seeing20_rm[indx[iexp]],format='(999a)')),2)
                  seeing50_target[itarget]=seeing50_target[itarget]+' '+strtrim(strcompress(string(seeing50_rm[indx[iexp]],format='(999a)')),2)
                  seeing80_target[itarget]=seeing80_target[itarget]+' '+strtrim(strcompress(string(seeing80_rm[indx[iexp]],format='(999a)')),2)
                  exptime_target[itarget]=exptime_target[itarget]+exptime_rm[indx[iexp]]
               endfor
            endif
            Firstcarton_target[itarget] = unique_plmap_values(Firstcarton_target[itarget])
            carton2TarPK_target[itarget] = unique_plmap_values(carton2TarPK_target[itarget])
         endif
   endfor
   mjd_t=mjd_t/snr_t
   mjdsfinal[*]=mjd_t

   airmass_target_f  = dblarr(ntarget)
   seeing20_target_f = dblarr(ntarget)
   seeing50_target_f = dblarr(ntarget)
   seeing80_target_f = dblarr(ntarget)

   for itar=0, ntarget-1 do begin
        weights_target_f_tmp   = double((strsplit(weights_target[itar], /extract)))
        airmass_target_f_tmp   = double((strsplit(airmass_target[itar], /extract)))
        seeing20_target_f_tmp  = double((strsplit(seeing20_target[itar],/extract)))
        seeing50_target_f_tmp  = double((strsplit(seeing50_target[itar],/extract)))
        seeing80_target_f_tmp  = double((strsplit(seeing80_target[itar],/extract)))
      
        airmass_target_f[itar]  = total(airmass_target_f_tmp *weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
        seeing20_target_f[itar] = total(seeing20_target_f_tmp*weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
        seeing50_target_f[itar] = total(seeing50_target_f_tmp*weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
        seeing80_target_f[itar] = total(seeing80_target_f_tmp*weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
   endfor
   
   mjdf_target_s.mjd_final=mjdsfinal
   finalplugmap=struct_addtags(finalplugmap,mjdf_target_s)
   moon_target_s.moon_dist=moon_target
   finalplugmap=struct_addtags(finalplugmap,moon_target_s)
   moon_phasef_s.moon_phase=moon_phasef
   finalplugmap=struct_addtags(finalplugmap,moon_phasef_s)
   fiber_target_s.fiberid_list=fiber_target
   finalplugmap=struct_addtags(finalplugmap,fiber_target_s)

   RA_target_s.RA_list=RA_target
   finalplugmap=struct_addtags(finalplugmap,RA_target_s)
   DEC_target_s.DEC_list=DEC_target
   finalplugmap=struct_addtags(finalplugmap,DEC_target_s)
   exptime_target_s.exptime=exptime_target
   finalplugmap=struct_addtags(finalplugmap,exptime_target_s)

   FIRSTCARTON_target_s.FIRSTCARTON_LIST=Firstcarton_target
   finalplugmap=struct_addtags(finalplugmap,FIRSTCARTON_target_s)
   cartoon2TarPK_target_s.CARTON_TO_TARGET_PK_LIST=carton2TarPK_target
   finalplugmap=struct_addtags(finalplugmap,cartoon2TarPK_target_s)
   Assigned_target_s.ASSIGNED_LIST=Assigned_target
   finalplugmap=struct_addtags(finalplugmap,Assigned_target_s)
   on_target_target_s.ON_TARGET_LIST=on_target_target
   finalplugmap=struct_addtags(finalplugmap,on_target_target_s)
   valid_target_s.VALID_LIST=on_target_target
   finalplugmap=struct_addtags(finalplugmap,valid_target_s)
   exp_disp_med_s.EXP_DISP_MED=exp_disp_med
   finalplugmap=struct_addtags(finalplugmap,exp_disp_med_s)
   xfocal_target_s.XFOCAL_LIST=xfocal_target
   finalplugmap=struct_addtags(finalplugmap,xfocal_target_s)
   yfocal_target_s.YFOCAL_LIST=yfocal_target
   finalplugmap=struct_addtags(finalplugmap,yfocal_target_s)
   
   tai_target_s.tai_list=tai_target
   finalplugmap=struct_addtags(finalplugmap,tai_target_s)
   snr2G_target_s.fieldsnr2g_list=snr2G_target
   finalplugmap=struct_addtags(finalplugmap,snr2G_target_s)
   snr2R_target_s.fieldsnr2r_list=snr2R_target
   finalplugmap=struct_addtags(finalplugmap,snr2R_target_s)
   snr2I_target_s.fieldsnr2i_list=snr2I_target
   finalplugmap=struct_addtags(finalplugmap,snr2I_target_s)

   mjdlist_target_s.mjdlist=mjdlist_target
   finalplugmap=struct_addtags(finalplugmap,mjdlist_target_s)
   designs_target_s.designs=designs_target
   finalplugmap=struct_addtags(finalplugmap,designs_target_s)
   configs_target_s.configs=configs_target
   finalplugmap=struct_addtags(finalplugmap,configs_target_s)
   
   airmass_s.AIRMASS=airmass_target_f
   airmass_s.AIRMASS_LIST=airmass_target
   finalplugmap=struct_addtags(finalplugmap,airmass_s)
   seeing20_s.SEEING20=seeing20_target_f
   seeing20_s.SEEING20_LIST=seeing20_target
   finalplugmap=struct_addtags(finalplugmap,seeing20_s)
   seeing50_s.SEEING50=seeing50_target_f
   seeing50_s.SEEING50_LIST=seeing50_target
   finalplugmap=struct_addtags(finalplugmap,seeing50_s)
   seeing80_s.SEEING80=seeing80_target_f
   seeing80_s.SEEING80_LIST=seeing80_target
   finalplugmap=struct_addtags(finalplugmap,seeing80_s)
   
   finalplugmap=struct_addtags(finalplugmap,replicate(create_struct('FIBER_RA',' '),ntarget))
   finalplugmap=struct_addtags(finalplugmap,replicate(create_struct('FIBER_DEC',' '),ntarget))
   finalplugmap.FIBER_RA=finalplugmap.RA
   finalplugmap.FIBER_DEC=finalplugmap.DEC
   ;---------------------------------------------------------------------------
   ; Write the corrected spCFrame files.
   ; All the fluxes + their errors are calibrated.
   ; The wavelengths + dispersions are converted from trace sets to 2D images.
   ; The pixel mask has the COMBINEREJ bit set.
   ;---------------------------------------------------------------------------

   if not keyword_set(onestep_coadd) then begin
           plugmap=comb_plugmap
           plugmap=in_plugmap
   endif

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
      
      fxaddpar, exthdr, 'EXTNAME', 'SPRESL', ' Spectral resolution'
      mwrfits, resolution[*,indx], thisfile, exthdr
      
   endfor
   splog, prename=''
   ;----------
   ; Clear memory

   wave = 0
   flux = 0
   fluxivar = 0
   temppixmask = 0
   dispersion = 0
   resolution = 0
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
   sxaddpar, bighdr, 'FIELDID', (strsplit(outputname,'-',/EXTRACT))[1]
   ;sxaddpar, bighdr, 'FIELDID', strmid(outputname,8,5)
   ncoeff = sxpar(bighdr, 'NWORDER')
   for i=2, ncoeff-1 do sxdelpar, bighdr, 'COEFF'+strtrim(string(i),2)

   sxdelpar, bighdr, ['SPA', 'IPA', 'IPARATE']
   sxdelpar, bighdr, 'CONFID'
   sxdelpar, bighdr, 'DESIGNID'
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

   cardname = [ 'AZ', 'ALT', 'AIRMASS', 'TAI', 'WTIME', 'AIRTEMP', 'DEWPOINT', $
    'DEWDEP', 'DUSTA', 'DUSTB', 'DUSTC', 'DUSTD', 'GUSTS', 'HUMIDITY', $
    'HUMIDOUT', 'PRESSURE', 'WINDD', 'WINDS', 'TEMP01', 'TEMP02', $
    'TEMP03', 'TEMP04', 'HELIO_RV', 'SEEING20', 'SEEING50', 'SEEING80', $
    'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
    'WSIGMA', 'XSIGMA' ]
   sxdelpar, bighdr,cardname
   sxdelpar, bighdr,'TAI-BEG'
   sxdelpar, bighdr,'TAI-END'
   sxdelpar, bighdr,'XCHI2'
   sxdelpar, bighdr,'SKYCHI2'
   sxdelpar, bighdr,'WSIGMA'
   sxdelpar, bighdr,'XSIGMA'
   sxdelpar, bighdr,'NGUIDE'

   cameras0 = sxpar(*(hdrarr[0]), 'CAMERAS')
   for ihdr=1, n_elements(hdrarr)-1 do begin
      if (sxpar(*(hdrarr[ihdr]), 'CAMERAS') EQ cameras0) then $
       ;sxcombinepar, hdrarr[ihdr], cardname, bighdr, func='total'
       sxdelpar, bighdr,cardname
   endfor

   ;----------
   ; Use the MJD passed as a keyword, which will typically be for the most
   ; observation, and be consistent with the output file names

   if (keyword_set(mjd)) then $
    sxaddpar, bighdr, 'MJD', mjd

   ; Get the list of MJD's used for these reductions, then convert to a string
   mjdlist = mjdlist[uniq(mjdlist, sort(mjdlist))]
   mjdlist = strtrim(strcompress(string(mjdlist,format='(999a)')),2)

   ; Get the list of Designs used for these reductions, then convert to a string
   designlist = strtrim(strcompress(string(designlist,format='(999a)')),2)
 
   ; Get the list of configurations used for these reductions, then convert to a string
   configlist = strtrim(strcompress(string(configlist,format='(999a)')),2)

   if keyword_set(tai_flag) then begin
    indtai=uniq(tailist, sort(tailist))
    tailist = tailist[indtai]
    snr2listG = snr2listG[indtai]
    snr2listR = snr2listR[indtai]
    snr2listI = snr2listI[indtai]
    tailist = strtrim(strcompress(string(string(tailist,format='(i15)'),format='(999a)')),2)
    snr2listG = strtrim(strcompress(string(snr2listG,format='(999a)')),2)
    snr2listR = strtrim(strcompress(string(snr2listR,format='(999a)')),2)
    snr2listI = strtrim(strcompress(string(snr2listI,format='(999a)')),2)
   endif
   ;----------
   ; Add new header cards

   sxaddpar, bighdr, 'VERSCOMB', idlspec2d_version(), $
    ' Version of idlspec2d for combining multiple spectra', after='VERS2D'
   sxaddpar, bighdr, 'NEXP', nfiles, $
    ' Number of exposures in this file', before='EXPTIME'
   for ifile=0,nfiles-1 do $
    sxaddpar, bighdr, string('EXPID',ifile+1, format='(a5,i2.2)'), label[ifile], $
     ' ID string for exposure '+strtrim(ifile+1,2), before='EXPTIME'
   if (keyword_set(bestexpnum)) then sxaddpar, bighdr, 'BESTEXP', bestexpnum, before='EXPID01'

   sxaddpar, bighdr, 'EXPTIME', min(exptimevec), $
    ' Minimum of exposure times for all cameras'
   for icam=0, ncam-1 do $
    sxaddpar, bighdr, 'NEXP_'+camnames[icam], nexpvec[icam], $
     ' '+camnames[icam]+' camera number of exposures', before='EXPTIME'
   for icam=0, ncam-1 do $
    sxaddpar, bighdr, 'EXPT_'+camnames[icam], exptimevec[icam], $
     ' '+camnames[icam]+' camera exposure time (seconds)', before='EXPTIME'
   sxaddpar, bighdr, 'SPCOADD', systime(), ' SPCOADD finished', after='EXPTIME'

   sxaddpar, bighdr, 'NWORDER', 2, ' Linear-log10 coefficients'
   sxaddpar, bighdr, 'NWORDER', 2, ' Linear-log10 coefficients'
   sxaddpar, bighdr, 'WFITTYPE', 'LOG-LINEAR', ' Linear-log10 dispersion'
   sxaddpar, bighdr, 'COEFF0', wavemin, ' Central wavelength (log10) of first pixel'
   sxaddpar, bighdr, 'COEFF1', binsz, ' Log10 dispersion per pixel'

   sxaddpar, bighdr, 'NAXIS1', n_elements(bestflux)
   sxaddpar, bighdr, 'NAXIS2', nfiber

   spawn, 'uname -n', uname
   sxaddpar, bighdr, 'UNAME', uname[0]

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
   ; Clean plugmap
   tags_to_delete= ['POSITIONERID','HOLEID', 'XWOK', 'YWOK', 'ZWOK', $
                    $;'XFOCAL', 'YFOCAL', 'ZFOCAL', $
                    'ALPHA', 'BETA', 'FIBERID']
   foreach tag, tags_to_delete do begin
      if tag_exist(finalplugmap,tag) then $
          finalplugmap = struct_trimtags(finalplugmap,except_tags=[tag])
   endforeach

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
   
   ; HDU #7 is the resolution map
   sxaddpar, hdrfloat, 'BUNIT', 'angstroms'
   sxaddpar, hdrfloat, 'EXTNAME', 'SPECRESL', ' Spectral resolution'
   mwrfits, finalresolution, fulloutname, hdrfloat
   
   if keyword_set(single_spectra) then begin
    ;writing each individual coadd spectrum on the field
    sxdelpar, bighdr, 'NAXIS2'
    spawn,'mkdir -p '+combinedir+'coadd'
    for itarget=0, ntarget-1 do begin
       added_exp=[]

       finalvalues=replicate(create_struct('flux',0.0),n_elements(finalwave))
       values_t=replicate(create_struct('loglam',0.0),n_elements(finalwave))
       finalvalues=struct_addtags(finalvalues,values_t)
       values_t=replicate(create_struct('ivar',0.0),n_elements(finalwave))
       finalvalues=struct_addtags(finalvalues,values_t)
       values_t=replicate(create_struct('and_mask',long(0)),n_elements(finalwave))
       finalvalues=struct_addtags(finalvalues,values_t)
       values_t=replicate(create_struct('or_mask',long(0)),n_elements(finalwave))
       finalvalues=struct_addtags(finalvalues,values_t)
       values_t=replicate(create_struct('wdisp',0.0),n_elements(finalwave))
       finalvalues=struct_addtags(finalvalues,values_t)
       values_t=replicate(create_struct('sky',0.0),n_elements(finalwave))
       finalvalues=struct_addtags(finalvalues,values_t)
       values_t=replicate(create_struct('wresl',0.0),n_elements(finalwave))
       finalvalues=struct_addtags(finalvalues,values_t)
       finalvalues.flux=finalflux[*,itarget]
       finalvalues.loglam=finalwave
       finalvalues.ivar=finalivar[*,itarget]
       finalvalues.and_mask=finalandmask[*,itarget]
       finalvalues.or_mask=finalormask[*,itarget]
       finalvalues.wdisp=finaldispersion[*,itarget]
       finalvalues.sky=finalsky[*,itarget]
       finalvalues.wresl=finalresolution[*,itarget]
       if keyword_set(legacy) or keyword_set(plates) then begin
          if keyword_set(legacy) then begin
             targid_tar=string(finalplugmap[itarget].fiberid,format='(i4.4)')
          endif else begin
             if finalplugmap[itarget].catalogid eq 0 then begin
                targid_tar=string(finalplugmap[itarget].fiberid,format='(i11.11)')
             endif else begin
                  if finalplugmap[itarget].program.contains('offset', /FOLD_CASE ) then begin
                      targid_tar=strtrim(string(finalplugmap[itarget].catalogid),1)
                  endif else begin
                      targid_tar=string(finalplugmap[itarget].catalogid,format='(i11.11)')
                  endelse
             endelse
          endelse
       endif else begin   
          targid_tar=finalplugmap[itarget].catalogid
       endelse
       sxaddpar, bighdr, 'PLUG_RA', final_ra[itarget], ' RA of Target'
       sxaddpar, bighdr, 'PLUG_DEC', final_dec[itarget], ' DEC of Target'
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
       coaddname = repstr(repstr(outputname,'spField','spSpec'),'.fits', $
        '-'+strtrim(targid_tar,2)+'.fits') 
       fulloutname_coadd = djs_filepath(coaddname, root_dir=coadddir)
       ; HDU # 0 header
       mwrfits, junk_d, fulloutname_coadd, bighdr, /create
       ; HDU # 1 header
       if itarget eq 0 then begin
         sxaddpar, coadd_val, 'EXTNAME', 'COADD', ' Coadded spectrum'
       endif
       mwrfits, finalvalues, fulloutname_coadd, coadd_val
       sxdelpar, coadd_val, 'COMMENT'
  
       ; HDU #2 is plugmap
       if itarget eq 0 then begin
         sxaddpar, hdrplug, 'EXTNAME', 'PLUGMAP', ' Plugmap structure'
       endif
       mwrfits, finalplugmap[itarget], fulloutname_coadd, hdrplug
       sxdelpar, hdrplug, 'COMMENT'
  
       for ifiber=0, nfiber-1 do begin
        for iexp=0, nexp_tmp - 1 do begin
         if keyword_set(legacy) or keyword_set(plates) then begin
           if keyword_set(legacy) then begin
              targid_rm=string(finalplugmap_rm[ifiber,iexp].fiberid,format='(i4.4)');targid_tar
           endif else begin
              if finalplugmap_rm[ifiber,iexp].catalogid eq 0 then begin
                  targid_rm=string(finalplugmap_rm[ifiber,iexp].fiberid,format='(i11.11)');targid_tar
              endif else begin
                  if finalplugmap[itarget].program.contains('offset', /FOLD_CASE ) then begin
                      targid_rm=strtrim(string(finalplugmap_rm[ifiber,iexp].catalogid),1)
                  endif else begin
                      targid_rm=string(finalplugmap_rm[ifiber,iexp].catalogid,format='(i11.11)');targid_tar
                  endelse
              endelse
           endelse
         endif else begin
           targid_rm=finalplugmap_rm[ifiber,iexp].catalogid
         endelse
         if targid_rm eq targid_tar then begin
           added_exp=[added_exp,iexp]
           finalvalues_rm=replicate(create_struct('flux',0.0),n_elements(finalwave))
           values_t=replicate(create_struct('loglam',0.0),n_elements(finalwave))
           finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
           values_t=replicate(create_struct('ivar',0.0),n_elements(finalwave))
           finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
           values_t=replicate(create_struct('and_mask',long(0)),n_elements(finalwave))
           finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
           values_t=replicate(create_struct('or_mask',long(0)),n_elements(finalwave))
           finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
           values_t=replicate(create_struct('wdisp',0.0),n_elements(finalwave))
           finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
           values_t=replicate(create_struct('sky',0.0),n_elements(finalwave))
           finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
           values_t=replicate(create_struct('wresl',0.0),n_elements(finalwave))
           finalvalues_rm=struct_addtags(finalvalues_rm,values_t)
           finalvalues_rm.flux=finalflux_rm[*,ifiber,iexp]
           finalvalues_rm.loglam=finalwave
           finalvalues_rm.ivar=finalivar_rm[*,ifiber,iexp]
           finalvalues_rm.and_mask=finalandmask_rm[*,ifiber,iexp]
           finalvalues_rm.or_mask=finalormask_rm[*,ifiber,iexp]
           finalvalues_rm.wdisp=finaldispersion_rm[*,ifiber,iexp]
           finalvalues_rm.sky=finalsky_rm[*,ifiber,iexp]
           finalvalues_rm.wresl=finalresolution_rm[*,ifiber,iexp]
           ; HDU # N header
           thisconf=mjds_rm_summ[ifiber,iexp];config_rm[ifiber,iexp]
           thisconf=string(thisconf,format='(i5.5)')+'-'+string(iexp,format='(i2.2)')
           
           indv_val=*hdrarr[iexp]
           sxdelpar,indv_val, ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND','CAMERAS']
           sxaddpar, indv_val, 'EXTNAME', 'MJD_EXP_'+thisconf, ' Single exposure spectrum'
           mwrfits, finalvalues_rm, fulloutname_coadd, indv_val
           sxdelpar, indv_val, 'COMMENT'
           
         endif
         
        endfor
      endfor

      if n_elements(added_exp) ne 0 then begin
        iused_hdrarr=hdrarr[added_exp]
        used_weight=weight[added_exp]
        cardnames_avg = ['AZ', 'ALT', 'AIRMASS', 'WTIME', 'AIRTEMP', 'DEWPOINT', $
                         'DEWDEP', 'DUSTA', 'DUSTB', 'DUSTC', 'DUSTD', 'GUSTS', 'GUSTD', $
                         'WINDD25M', 'WINDS25M', $
                         'HUMIDITY', 'HUMIDOUT', 'PRESSURE', 'WINDD', 'WINDS', 'TEMP01', 'TEMP02', $
                         'TEMP03', 'TEMP04', 'HELIO_RV', 'SEEING20', 'SEEING50', 'SEEING80', $
                         'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
                         'WSIGMA', 'XSIGMA' , 'CCDTEMP', 'LN2TEMP']
    
        h = headfits(fulloutname_coadd) ;Read primary header
    
        if n_elements(added_exp) eq 1 then begin
          foreach card, cardname do begin
            val= SXPAR(*iused_hdrarr[0], card, COMMENT = cmt)
            sxaddpar, h, card, val, cmt, /SaveComment
          endforeach
        endif else begin
           foreach cardname, cardnames_avg do begin
               sxcombinepar_v2, iused_hdrarr, cardname, h, Comment=key_match_dict[cardname], func='average', weights=used_weight, /SaveComment
           endforeach
           sxcombinepar_v2, iused_hdrarr, 'TAI-BEG', h, Comment=key_match_dict['TAIBEG'], func='min'
           sxcombinepar_v2, iused_hdrarr, 'TAI-END', h, Comment=key_match_dict['TAIEND'], func='max'
           
           sxcombinepar_v2, iused_hdrarr, 'XCHI2', h, Comment=key_match_dict['XCHI2MAX'], func='max', outcard='XCHI2MAX', after='XCHI2'
           sxcombinepar_v2, iused_hdrarr, 'XCHI2', h, Comment=key_match_dict['XCHI2MIN'], func='min', outcard='XCHI2MIN', after='XCHI2'
           sxcombinepar_v2, iused_hdrarr, 'SKYCHI2', h, Comment=key_match_dict['SCHI2MAX'], func='max', outcard='SCHI2MAX', after='SKYCHI2'
           sxcombinepar_v2, iused_hdrarr, 'SKYCHI2', h, Comment=key_match_dict['SCHI2MIN'], func='min', outcard='SCHI2MIN', after='SKYCHI2'
           sxcombinepar_v2, iused_hdrarr, 'WSIGMA', h, Comment=key_match_dict['WSIGMAX'], func='max', outcard='WSIGMAX', after='WSIGMA'
           sxcombinepar_v2, iused_hdrarr, 'WSIGMA', h, Comment=key_match_dict['WSIGMIN'], func='min', outcard='WSIGMIN', after='WSIGMA'
           sxcombinepar_v2, iused_hdrarr, 'XSIGMA', h, Comment=key_match_dict['XSIGMAX'], func='max', outcard='XSIGMAX', after='XSIGMA'
           sxcombinepar_v2, iused_hdrarr, 'XSIGMA', h, Comment=key_match_dict['XSIGMIN'], func='min', outcard='XSIGMIN', after='XSIGMA'
           sxcombinepar_v2, iused_hdrarr, 'NGUIDE', h, Comment=key_match_dict['NGUIDE'], func='total'
           sxcombinepar_v2, iused_hdrarr, 'EXPTIME', h, Comment=key_match_dict['EXPTIME'], func='total'
           sxaddpar, h, 'NEXP', n_elements(added_exp), key_match_dict['NEXP']
           sxdelpar, h, ['date-obs','SHOPETIM', 'SHCLOTIM', 'ionpump']

           
        endelse
        modfits,fulloutname_coadd,0,h ;Update header only
      endif

    endfor
   endif
   return
end
;------------------------------------------------------------------------------
