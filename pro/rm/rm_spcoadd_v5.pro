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
function plmap_arr_to_strlist, arr, format=format, unique=unique
    if keyword_set(format) then begin
        sarr = strtrim(strcompress(string(arr, format=format)),2)
    endif else sarr = strtrim(strcompress(strtrim(arr, 2)), 2)

    if keyword_set(unique) then sarr = sarr[UNIQ(sarr, SORT(sarr))]

    sarr = strjoin(sarr, ' ')  ; Join array into a space-separated string

    return, sarr
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
function build_fmaprow, itarget, plugmap, expidx, combinedfibermask, ti ;$

         plugmap.fibermask = combinedfibermask
        
         plugmap.moon_dist = plmap_arr_to_strlist(ti.moon_target_rm[expidx],format='(f0.1)')
         plugmap.MOON_PHASE = plmap_arr_to_strlist(ti.moon_phasef_rm[expidx],format='(f0.2)')
         plugmap.FIBERID_LIST = plmap_arr_to_strlist(ti.fiberid_rm[expidx])
         plugmap.RA_LIST = plmap_arr_to_strlist(ti.finalra_rm[expidx],format='(f0.6)')
         plugmap.DEC_LIST = plmap_arr_to_strlist(ti.finaldec_rm[expidx],format='(f0.6)')
         plugmap.DELTA_RA_LIST = plmap_arr_to_strlist(ti.finaldra_rm[expidx])
         plugmap.DELTA_DEC_LIST = plmap_arr_to_strlist(ti.finalddec_rm[expidx])
         plugmap.FIRSTCARTON_LIST = plmap_arr_to_strlist(ti.firstcarton_rm[expidx], /unique)
         plugmap.CARTON_TO_TARGET_PK_LIST = plmap_arr_to_strlist(ti.carton2TarPK_rm[expidx], /unique)
         plugmap.ASSIGNED_LIST = plmap_arr_to_strlist(ti.Assigned_rm[expidx])
         plugmap.ON_TARGET_LIST = plmap_arr_to_strlist(ti.on_target_rm[expidx])
         plugmap.VALID_LIST = plmap_arr_to_strlist(ti.valid_rm[expidx])
         plugmap.DECOLLIDED_LIST = plmap_arr_to_strlist(ti.DECOLLIDED_rm[expidx])
         plugmap.TOO_LIST = plmap_arr_to_strlist(ti.TOO_rm[expidx])
         plugmap.XFOCAL_LIST = plmap_arr_to_strlist(ti.xfocal_rm[expidx],format='(f0.3)')
         plugmap.YFOCAL_LIST = plmap_arr_to_strlist(ti.yfocal_rm[expidx],format='(f0.3)')
         plugmap.TAI_LIST = plmap_arr_to_strlist(ti.tai_rm[expidx],format='(i15)')
         plugmap.MJDLIST = plmap_arr_to_strlist(ti.mjds_rm[expidx],format='(i5)')
         plugmap.DESIGNS = plmap_arr_to_strlist(ti.designs[expidx])
         plugmap.CONFIGS = plmap_arr_to_strlist(ti.config_rm[expidx])
         plugmap.AIRMASS_LIST = plmap_arr_to_strlist(ti.airmass_rm[expidx])
         plugmap.SEEING20_LIST = plmap_arr_to_strlist(ti.seeing20_rm[expidx])
         plugmap.SEEING50_LIST = plmap_arr_to_strlist(ti.seeing50_rm[expidx])
         plugmap.SEEING80_LIST = plmap_arr_to_strlist(ti.seeing80_rm[expidx])
         plugmap.RMSOFF20_LIST = plmap_arr_to_strlist(ti.rmsoff20_rm[expidx])
         plugmap.RMSOFF50_LIST = plmap_arr_to_strlist(ti.rmsoff50_rm[expidx])
         plugmap.RMSOFF80_LIST = plmap_arr_to_strlist(ti.rmsoff80_rm[expidx])
         
         plugmap.AIRMASS_pt  = ptr_new(ti.airmass_rm[expidx])
         plugmap.SEEING20_pt = ptr_new(ti.seeing20_rm[expidx])
         plugmap.SEEING50_pt = ptr_new(ti.seeing50_rm[expidx])
         plugmap.SEEING80_pt = ptr_new(ti.seeing80_rm[expidx])
         plugmap.RMSOFF20_pt = ptr_new(ti.rmsoff20_rm[expidx])
         plugmap.RMSOFF50_pt = ptr_new(ti.rmsoff50_rm[expidx])
         plugmap.RMSOFF80_pt = ptr_new(ti.rmsoff80_rm[expidx])
         
         plugmap.EXPTIME  = TOTAL(ti.exptime_rm[expidx])
         plugmap.tai_pt  = ptr_new(ti.tai_rm[expidx])
         plugmap.fiberid_pt  = ptr_new(ti.fiberid_rm[expidx])
         plugmap.exp_pt  = ptr_new(ti.expid_rm[expidx])
         plugmap.nexp = n_elements(expidx)
         plugmap.target_index = itarget+1
         
         return, plugmap
end
;------------------------------------------------------------------------------


function masked_mjd, mjd, sn2s
    sn2 = float(sn2s)
    idx = where(sn2 gt 0,ct)
    mjd_final = 0
    if ct gt 0 then begin
        m_mjd = mjd[idx]
        m_sn2 = sn2[idx]
        mjd_final = Total(mjd*sn2)/Total(sn2)
    endif else begin
        mjd_final = mean(mjd)
    endelse
    return, mjd_final
end
;------------------------------------------------------------------------------

pro rm_spcoadd_v5, spframes, outputname, obs=obs, $
 mjd=mjd, binsz=binsz, zeropoint=zeropoint, nord=nord, $
 wavemin=wavemin, wavemax=wavemax, $
 bkptbin=bkptbin, window=window, maxsep=maxsep, adderr=adderr, $
 docams=camnames, plotsnfile=plotsnfile, combinedir=combinedir, $
 bestexpnum=bestexpnum,nofcorr=nofcorr,nodist=nodist, save_novalid=save_novalid, $
 plates=plates, legacy=legacy, single_spectra=single_spectra, epoch=epoch,$
 radec_coadd=radec_coadd, no_reject=no_reject,  onestep_coadd=onestep_coadd

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
      if strmatch(obs,'APO',/fold_case) eq 1 then begin
        if NOT keyword_set(camnames) then camnames = ['b1', 'r1']
      endif else begin
        if NOT keyword_set(camnames) then camnames = ['b2', 'r2']
      endelse
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
                            'SEEING20',0.D,'SEEING50',0.D,'SEEING80',0.D, $
                            'RMSOFF20',0.D,'RMSOFF50',0.D,'RMSOFF80',0.D,'itai','')
   if keyword_set(legacy) then begin
     nexp_tmp2 = nfiles/4 ;Get data for each exposure
   endif else begin
     nexp_tmp2 = nfiles/2
   endelse
   rm_plugmap = replicate(plugmap_rm, nexp_tmp2)
   if keyword_set(epoch) then subdir = '..'
   for ifile=0, nfiles-1 do begin
     print, djs_filepath(filenames[ifile], subdirectory=subdir)
      spframe_read, djs_filepath(filenames[ifile], subdirectory=subdir), hdr=objhdr
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
        rm_plugmap[ifile].itai=strtrim(string(sxpar(objhdr,'TAI-BEG')+(double(sxpar(objhdr,'EXPTIME'))/2.d),format='(i15)'),2)
        rm_plugmap[ifile].mjd=sxpar(objhdr,'MJD')
        rm_plugmap[ifile].airmass=sxpar(objhdr,'AIRMASS')
        rm_plugmap[ifile].SEEING20=sxpar(objhdr,'SEEING20')
        rm_plugmap[ifile].SEEING50=sxpar(objhdr,'SEEING50')
        rm_plugmap[ifile].SEEING80=sxpar(objhdr,'SEEING80')
        rm_plugmap[ifile].RMSOFF20=sxpar(objhdr,'RMSOFF20')
        rm_plugmap[ifile].RMSOFF50=sxpar(objhdr,'RMSOFF50')
        rm_plugmap[ifile].RMSOFF80=sxpar(objhdr,'RMSOFF80')
        mjdt=sxpar(objhdr,'TAI-BEG')/(24.D*3600.D)
        mjd2datelist,mjdt,datelist=date
        rm_plugmap[ifile].date=date
        rm_plugmap[ifile].exptime=sxpar(objhdr,'EXPTIME')
      endif
   endfor
   npixmax = max(npixarr)
   nobj = sxpar(objhdr,'NAXIS2') ; Number of fibers per spectrograph

 ;struct_print, rm_plugmap, filename='rm_plugmap.html', /html
   for ifile=0, nfiles-1 do begin
      ;----------
      ; Read in all data from this input file.
      ; Reading the plug-map structure will fail if its structure is
      ; different between different files.

      splog, 'Reading file #', ifile, ': ', filenames[ifile], $
         prename=filenames[ifile]
      spframe_read, djs_filepath(filenames[ifile], subdirectory=subdir), objflux=tempflux, objivar=tempivar, $
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
      if not keyword_set(thisdesign) then thisdesign = sxpar(hdr,'PLATEID')
      thisconfig = sxpar(hdr,'CONFID')
      if not keyword_set(thisconfig) then thisconfig = sxpar(hdr, 'NAME')
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
         resolution = make_array(npixmax,nobj*nfiles,type=size(tempresolution,/type))
         pixelmask = make_array(npixmax,nobj*nfiles,type=size(temppixmask,/type))
         skyflux = make_array(npixmax,nobj*nfiles,type=size(tempsky,/type))
         ximg = make_array(npixmax,nobj*nfiles,type=size(tempximg,/type))
         superflat = make_array(npixmax,nobj*nfiles,type=size(tempsuperflat,/type))

         ; Append as vectors...
         camerasvec = cameras
         label = makelabel(hdr)
         filenum = lonarr(nfib) + ifile
         plugmap = tempplug
         if not tag_exist(plugmap,'fieldCadence') then Field_cadence = '' $
         else begin
            fc = plugmap.fieldCadence
            Field_cadence = strjoin(fc[UNIQ(fc, sort(fc))], ' ')
         endelse
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

         if keyword_set(plugmap) then begin
             nflags = max([n_elements(plugmap[0].SDSS5_TARGET_FLAGS),n_elements(tempplug[0].SDSS5_TARGET_FLAGS)])
             plugmap = rename_tags(plugmap, 'SDSS5_TARGET_FLAGS','SDSS5_TARGET_FLAGS_raw')
             tempplug = rename_tags(tempplug, 'SDSS5_TARGET_FLAGS','SDSS5_TARGET_FLAGS_raw')
            plugmap = struct_addtags(plugmap, $
                                     replicate(create_struct('SDSS5_TARGET_FLAGS', BYTARR(nflags)),$
                                               n_elements(plugmap)))
            plugmap.SDSS5_TARGET_FLAGS[0:n_elements(plugmap[0].SDSS5_TARGET_FLAGS_raw)-1] = plugmap.SDSS5_TARGET_FLAGS_raw
            plugmap = struct_trimtags(plugmap, except_tags='SDSS5_TARGET_FLAGS_RAW')
            tempplug = struct_addtags(tempplug, $
                                      replicate(create_struct('SDSS5_TARGET_FLAGS', BYTARR(nflags)),$
                                               n_elements(tempplug)))
            tempplug.SDSS5_TARGET_FLAGS[0:n_elements(tempplug[0].SDSS5_TARGET_FLAGS_raw)-1] = tempplug.SDSS5_TARGET_FLAGS_raw
            tempplug = struct_trimtags(tempplug, except_tags='SDSS5_TARGET_FLAGS_RAW')
        endif

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

    comb_hdrs       = []

   ;----------
   ; Scale the blue and red flux to the same flux level
   ; note that the filenames are sorted as b1[nexp],r1[nexp]
   if keyword_set(legacy) then begin
     nexp_tmp = nfiles/4
     for iexp=0, nexp_tmp - 1 do begin
       for iobj=0L,499L do begin
         ; for b1 and r1
         ifile = iexp
         comb_hdrs = [comb_hdrs, hdrarr[ifile]]
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
            ; for b and r
            ifile = iexp
            comb_hdrs = [comb_hdrs, hdrarr[ifile]]

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
            ifile = iexp + nexp_tmp
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


   if keyword_set(onestep_coadd) then shape = [nfiber, nexp_tmp] else shape = nfiber*nexp_tmp
        
       ti = {finalra_rm: fltarr(shape), finaldec_rm: fltarr(shape),$
             finaldra_rm: fltarr(shape), finalddec_rm: fltarr(shape),$
             fiberid_rm: lonarr(shape),$
             firstcarton_rm: strarr(shape), carton2TarPK_rm: LON64ARR(shape),$
             Assigned_rm: lonarr(shape), on_target_rm: lonarr(shape),$
             valid_rm: lonarr(shape), DECOLLIDED_rm: lonarr(shape),$
             TOO_rm: lonarr(shape),$
             xfocal_rm: fltarr(shape), yfocal_rm: fltarr(shape),$
             mjds_rm: lonarr(shape), mjds_rm_summ: lonarr(nfiber, nexp_tmp),$
             config_rm: lonarr(shape), designs:strarr(shape),$
             tai_rm: dblarr(shape), exptime_rm: fltarr(shape),$
             moon_target_rm:strarr(shape), moon_phasef_rm:strarr(shape),$
             snr2listG:strarr(nfiber, nexp_tmp), snr2listR:strarr(nfiber, nexp_tmp),$
             snr2listI:strarr(nfiber, nexp_tmp), airmass_rm:strarr(shape),$
             seeing20_rm:strarr(shape), seeing50_rm:strarr(shape),$
             seeing80_rm:strarr(shape), rmsoff20_rm:strarr(shape),$
             rmsoff50_rm:strarr(shape), rmsoff80_rm:strarr(shape),$
             expid_rm: lonarr(shape)}

   ;----------
   
   
   combinedwave    = dblarr(nfinalpix,nobj*nexp_tmp)
   combinedflux    = fltarr(nfinalpix,nobj*nexp_tmp)
   combinedivar    = fltarr(nfinalpix,nobj*nexp_tmp)
   combineddisp    = fltarr(nfinalpix,nobj*nexp_tmp)
   combinedskyflux = fltarr(nfinalpix,nobj*nexp_tmp)
   combinedresl    = fltarr(nfinalpix,nobj*nexp_tmp)
   combinedandmask = lonarr(nfinalpix,nobj*nexp_tmp)
   combinedormask  = lonarr(nfinalpix,nobj*nexp_tmp)
   combinedfibermask = LON64ARR(nobj,nexp_tmp)
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
         ffm = 0
         
         foreach indx1, indx do ffm = ffm or plugmap[indx1].fibermask
         plugmap[indx].fibermask = ffm
         
         ; The Following are used as input for the target centric coadd
         combinedwave[*,i]    = finalwave
         combinedflux[*,i]    = bestflux
         combinedivar[*,i]    = bestivar
         combineddisp[*,i]    = bestdispersion
         combinedskyflux[*,i] = bestsky
         combinedresl[*,i]    = bestresolution
         combinedandmask[*,i] = bestandmask
         combinedormask[*,i]  = bestormask
         combinedfibermask[i] = ffm
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
        if plugmap[indx[0]].delta_ra eq 0.0 then $ ;
             dratemp = string(plugmap[indx[0]].delta_ra,format='(f0.1)') $
        else dratemp = string(plugmap[indx[0]].delta_ra,format='(f0.6)')
        if plugmap[indx[0]].delta_dec eq 0.0 then $ ;
             ddectemp = string(plugmap[indx[0]].delta_dec,format='(f0.1)') $
        else ddectemp = string(plugmap[indx[0]].delta_dec,format='(f0.6)')
        moon_dist = djs_diff_angle(ra_moon, dec_moon, ratemp, dectemp)


         if keyword_set(onestep_coadd) then fiexp = [ifiber,iexp] else fiexp=[i-1,0]
         
         ti.finalra_rm[fiexp]=ratemp
         ti.finaldec_rm[fiexp]=dectemp
         ti.finaldra_rm[fiexp]=dratemp
         ti.finalddec_rm[fiexp]=ddectemp
         ti.mjds_rm[fiexp]=rm_plugmap[iexp].mjd
         ti.mjds_rm_summ[ifiber, iexp]=rm_plugmap[iexp].mjd
         ti.fiberid_rm[fiexp] = plugmap[indx[0]].fiberid
         ti.firstcarton_rm[fiexp] = plugmap[indx[0]].firstcarton
         ti.carton2TarPK_rm[fiexp] = plugmap[indx[0]].carton_to_target_pk
         ti.Assigned_rm[fiexp] = plugmap[indx[0]].assigned
         ti.on_target_rm[fiexp] = plugmap[indx[0]].on_target
         ti.valid_rm[fiexp] = plugmap[indx[0]].valid
         ti.DECOLLIDED_rm[fiexp] = plugmap[indx[0]].decollided
         ti.TOO_rm[fiexp] = plugmap[indx[0]].too
         ti.xfocal_rm[fiexp] = plugmap[indx[0]].xfocal
         ti.yfocal_rm[fiexp] = plugmap[indx[0]].yfocal
         ti.tai_rm[fiexp]=rm_plugmap[iexp].tai+double(rm_plugmap[iexp].exptime/2.0)
         ti.exptime_rm[fiexp]=rm_plugmap[iexp].exptime
         ; use expuse number instad of configuration number for legacy
         ti.config_rm[fiexp]=rm_plugmap[iexp].configuration
         ti.designs[fiexp] =rm_plugmap[iexp].design
         ti.airmass_rm[fiexp] = rm_plugmap[iexp].airmass
         ti.seeing20_rm[fiexp] =rm_plugmap[iexp].SEEING20
         ti.seeing50_rm[fiexp] =rm_plugmap[iexp].SEEING50
         ti.seeing80_rm[fiexp] =rm_plugmap[iexp].SEEING80
         ti.rmsoff20_rm[fiexp] =rm_plugmap[iexp].RMSOFF20
         ti.rmsoff50_rm[fiexp] =rm_plugmap[iexp].RMSOFF50
         ti.rmsoff50_rm[fiexp] =rm_plugmap[iexp].RMSOFF50
         ti.expid_rm[fiexp] = iexp
         ti.moon_target_rm[fiexp] =moon_dist
         ti.moon_phasef_rm[fiexp] =mfrac

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
      fname_ps_list = []
      for iexp=0, nexp_tmp - 1 do begin
        splog, 'EXPOSURE number ', expnumf[iexp]
        
        ; Plot S/N and throughput **before** this distortion-correction.
        splog, prelog='Initial'
        platesn, finalflux_rm[*,*,iexp], finalivar_rm[*,*,iexp], $
          finalandmask_rm[*,*,iexp], finalplugmap_rm[*,iexp], finalwave, obs=obs, $
          hdr=bighdr, legacy=legacy, plotfile=djs_filepath(repstr(plotsnfile+'.orig','X',string(iexp,format='(i2.2)')), root_dir=combinedir)
        splog, prelog=''
        
        splog, prelog='FluxDistort '+strtrim(iexp,2)
        fname_tmp=repstr(distortfitsfile,'.fits','-'+strtrim(string(expnumf[iexp],f='(i010.8)'),2)+'.fits')
        fname_ps_tmp=repstr(distortpsfile,'.ps','-'+strtrim(string(expnumf[iexp],f='(i010.8)'),2)+'.ps')
        corrimg = flux_distortion(finalflux_rm[*,*,iexp], finalivar_rm[*,*,iexp], $
          finalandmask_rm[*,*,iexp], finalormask_rm[*,*,iexp], $
          plugmap=finalplugmap_rm[*,iexp], loglam=finalwave, $
          plotfile=fname_ps_tmp, hdr=*(hdrarr[iexp]),legacy=legacy, /oneexp)

        if keyword_set(fname_ps_tmp) then fname_ps_list = [fname_ps_list,fname_ps_tmp]

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
        
        ; First write the file with the flux distortion vectors
        dist_hdr = *comb_hdrs[iexp]
        
        del_card = ['FILENAME', 'CAMERAS', 'CCD', 'CCDID', 'CCDTYPE']
        foreach card, del_card do sxdelpar, dist_hdr, card
        
        mwrfits_named, corrimg,fname_tmp , hr=dist_hdr, name='CORRIMG', /create

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
        finalflux_rm[*,*,iexp]=final_flux
        finalivar_rm[*,*,iexp]=final_ivar
        finalsky_rm[*,*,iexp]=final_sky
        finalandmask_rm[*,*,iexp]=finalandmask_t
        finalormask_rm[*,*,iexp]=finalormask_t


        i = where(ti.expid_rm eq iexp)
        fid = ti.fiberid_rm[i]-1
        combinedflux[*,i]    = final_flux[*,fid]
        combinedivar[*,i]    = final_ivar[*,fid]
        combinedskyflux[*,i] = final_sky[*,fid]
        combinedandmask[*,i] = finaland_mask[*,fid]
        combinedormask[*,i]  = finalor_mask[*,fid]


      endfor
      if n_elements(fname_ps_list) gt 0 then begin
          cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=ps2write -sOutputFile='+distortpsfile+' '+strjoin(fname_ps_list,' ')
          splog, 'SPAWN '+cmd
          spawn, cmd, sh_out, sh_err
          splog, 'SPAWN out=', sh_out
          splog, 'SPAWN err=', sh_err
          ps2pdf, distortpsfile
      endif

   endif
   splog, prelog='Final'

   if not keyword_set(onestep_coadd) then begin
        in_plugmap=plugmap
        plugmap=comb_plugmap
   endif
   indx_tar = []
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
             indx_tar=[ptr_new(nt1)]
           endif else begin
             indx_tar=[indx_tar,ptr_new(nt1)]
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
       
       ; Preallocate a string array of the same size
       sid_tmp = STRARR(N_ELEMENTS(plugmap.sdss_id))

       ; Convert each long64 element to string
       FOR i=0, N_ELEMENTS(plugmap.sdss_id)-1 DO BEGIN
          sid_tmp[i] = STRTRIM(STRING(plugmap[i].sdss_id), 2)
       ENDFOR
       
       catid_rm = plugmap.catalogid
       catid_ix = intarr(n_elements(catid_rm))
       cix = where(plugmap.sdss_id lt 0, ct)
       if ct gt 0 then sid_tmp[cix] = plugmap[cix].catalogid
       
       srt = SORT(sid_tmp)
       sid_tmp_id = UNIQ(sid_tmp, srt)
       
       plugmap = struct_addtags(plugmap, $
                                replicate(create_struct('catalogid_umod', string('')),$
                                          n_elements(plugmap)))
       plugmap.catalogid_umod = plugmap.catalogid
       
       foreach usid_ix, sid_tmp_id do begin
;            print, usid_ix
            sid_ix = where(sid_tmp eq sid_tmp[usid_ix])
;            print, sid_ix
;            print,sid_tmp[sid_ix]
            cid_vals = catid_rm[sid_ix]
            catid_rm[sid_ix] = max(catid_rm[sid_ix], xidx)
            if strmatch(strtrim(sid_tmp[usid_ix],2), '70099948') then begin
;                print, sid_ix
;                print, plugmap[sid_ix].catalogid
;                print, xidx
;                message, 'test'
            endif
            catid_ix[sid_ix] = xidx[0]
            plugmap[sid_ix].catalogid = max(catid_rm[sid_ix])
;#TODO rest of plugmap????
       endforeach

       ;catid_rm = plugmap.catalogid
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
            if not keyword_set(indx_tar) then indx_tar=[ptr_new(indx1)] $
            else indx_tar=[indx_tar,ptr_new(indx1)]
         endif
       endforeach
   endelse
   
   ntarget=n_elements(indx_tar)
   
   finalflux = []
   finalivar = []
   finalandmask = []
   finalormask = []
   finaldispersion = []
   finalresolution = []
   finalsky = []
   finalplugmap = 0
   weight = dblarr(nexp_tmp)
   mjd_targ = []
   ; FINALPLUGMAP structure.

   new_cols = {target_index:0, FIBERID_LIST:' ', FIBER_RA:0d, FIBER_DEC:0d,$
               RA_LIST:' ', DEC_LIST:' ', DELTA_RA_LIST:' ', DELTA_DEC_LIST:' ',$
               FIRSTCARTON_LIST:' ', CARTON_TO_TARGET_PK_LIST:' ', ASSIGNED_LIST:' ',$
               ON_TARGET_LIST:' ', VALID_LIST:' ', DECOLLIDED_LIST:' ',TOO_LIST:' ',$
               EXP_DISP_MED: 0.D, XFOCAL_LIST:' ', YFOCAL_LIST:' ', exptime:0, nexp:0,$
               MJD_FINAL: 0.D, TAI_LIST:' ', MJDLIST:' ', DESIGNS:' ', CONFIGS:' ',$
               MOON_DIST:' ',MOON_PHASE:' ',AIRMASS: 0.D, AIRMASS_LIST:' ', $
               SEEING20: 0.D, SEEING20_LIST: ' ', SEEING50: 0.D, SEEING50_LIST: ' ', $
               SEEING80: 0.D, SEEING80_LIST: ' ', RMSOFF20: 0.D, RMSOFF20_LIST: ' ', $
               RMSOFF50: 0.D, RMSOFF50_LIST: ' ', RMSOFF80: 0.D, RMSOFF80_LIST: ' ',$
               FIELDSNR2G_LIST:' ', FIELDSNR2R_LIST:' ', FIELDSNR2I_LIST:' ',$
               tai_pt:ptr_new(), fiberid_pt:ptr_new(), exp_pt:ptr_new(), $
               weight_pt:ptr_new(), AIRMASS_pt:ptr_new(),$
               SEEING20_pt:ptr_new(),SEEING50_pt:ptr_new(),SEEING80_pt:ptr_new(),$
               RMSOFF20_pt:ptr_new(),RMSOFF50_pt:ptr_new(),RMSOFF80_pt:ptr_new() }
    
   exclude_mask = (fibermask_bits('NOPLUG') OR $
                   fibermask_bits('BADTRACE') OR $
                   fibermask_bits('BADFLAT') OR $
                   fibermask_bits('BADARC') OR $
                   fibermask_bits('NODATA'))
   if keyword_set(no_reject) then begin
        ctype = "Two step coadd with no rejection of "
   endif else begin
        if not keyword_set(onestep_coadd) then begin
            ctype = 'Two step coadd with rejection of '
        endif else begin
            ctype = 'One step (legacy) coadd with rejection of '
        endelse
   endelse
 
   if keyword_set(radec_coadd) then begin
        c1type = 'exposures with the same coordinates'
   endif else c1type = 'exposures with the same CatalogID'
   
   itarget = -1
   for itarg=0, ntarget-1 do begin
      indx=*(indx_tar[itarg])
      
      mask = WHERE((plugmap[indx].fibermask AND exclude_mask) eq 0, nkeep, NCOMPLEMENT=ndrop)
      if nkeep eq 0 then begin
        targid_tar=plugmap[indx[0]].catalogid
        IF STRPOS(targid_tar, 'u', 0) NE 0 THEN begin
            indx[0] = -1
        endif
      endif
      
      if (indx[0] NE -1) then begin
      
         expidx = indx
         if keyword_set(onestep_coadd) then expidx = indx[0:(n_elements(indx)/2)-1]

         if nkeep gt 0 then begin
            indx = indx[mask]
            mask1 = WHERE((plugmap[expidx].fibermask AND exclude_mask) eq 0, nkeep1, NCOMPLEMENT=ndrop1)
            expidx = expidx[mask1]
         endif
         
      
      
         itarget = itarget+1
         splog, ctype+'all the ('+strtrim(n_elements(expidx),2)+') '+c1type+' ('+strtrim(itarg+1,2)+'/'+strtrim(ntarget,2)+')'
         if keyword_set(radec_coadd) then begin
           splog, 'Target', itarg+1, ' ', plugmap[indx[0]].objtype, $
                    plugmap[indx[0]].mag, format = '(a, i5.4, a, a, f6.2, 5f6.2)'
         endif else begin
           if not strmatch('u*', plugmap[indx[0]].catalogid) then begin
               catid =  max(plugmap[expidx].catalogid, xidx)
               catid_ix[expidx] = xidx[0]
           endif else begin
                catid = plugmap[indx[0]].catalogid
           endelse

           splog, 'Catalogid: ', strtrim(catid,2), ' ', plugmap[indx[0]].objtype, $
                    plugmap[indx[0]].mag, format = '(a, a, a, a, f6.2, 5f6.2)'
         endelse

         
         if nkeep gt 0 then begin
            if ndrop gt 0 then splog,'Dropping '+strtrim(ndrop,2)+' exposures with fibermask != 0'
         endif
    
         if keyword_set(ra_coadd) then begin
            thisplug = plugmap[expidx[0]]
         endif else begin
            ix = catid_ix[expidx[0]]
            print, expidx, ix
            thisplug = plugmap[expidx[ix]]
         endelse
;         help, plugmap
;         help, expidx
;         help, catid_ix
;         print, expidx
;         message, 'test'
         
         thisplug = struct_addtags(thisplug,new_cols)

         ffm = 0
         foreach eindx,expidx do begin
            if keyword_set(onestep_coadd) then begin
                ffm = ffm OR plugmap[eindx].fibermask
            endif else ffm = ffm OR combinedfibermask[eindx]
         endforeach

         thisplug = build_fmaprow(itarget, thisplug, expidx, ffm, ti)

         if not keyword_set(finalplugmap) then finalplugmap = thisplug $
         else finalplugmap = [finalplugmap, thisplug]

         bestresolution = 0
         if keyword_set(no_reject) or (not keyword_set(onestep_coadd)) then begin
           if keyword_set(no_reject) then begin
                maxiter=0
                upper=3d6
                lower=3d6
                maxrej=1
           endif else begin
                maxiter=50
                upper=3.0
                lower=3.0
                maxrej=1
           endelse
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
                     maxiter=maxiter, upper=upper, lower=lower, maxrej=maxrej
                     
           ; only calculate and report median exposure dispersion if valid target
           if (not keyword_set(legacy)) and (not keyword_set(plates)) then begin
            test=where((plugmap[indx].Assigned eq 1) $
                         AND (plugmap[indx].valid eq 1) $
                         AND (plugmap[indx].On_target eq 1), ctv)
           endif else ctv=1

           if ctv ne 0 then begin
             if n_elements(indx) gt 1 then begin
                finalplugmap[itarget].exp_disp_med = abs(median(STDDEV(combinedflux,DIMENSION=2,/double)/bestflux,$
                                                         /EVEN,/double))
             endif else finalplugmap[itarget].exp_disp_med = 0.0
           endif else finalplugmap[itarget].exp_disp_med = -1.d
         endif else begin
           ; This is for legacy to match older coadd algorithm
           ; 1 step coadd
           temppixmask = pixelmask[*,indx]
           combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
                     finalmask=temppixmask, indisp=dispersion[*,indx], $
                     skyflux=skyflux[*,indx], $
                     newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
                     andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
                     newsky=bestsky, inresl=resolution[*,indx], newresl=bestresolution, $
                     nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
                     maxiter=50, upper=3.0, lower=3.0, maxrej=1
            finalplugmap[itarget].exp_disp_med = -1.

           ; The following adds the COMBINEREJ bit to the input pixel masks
           pixelmask[*,indx] = temppixmask
         endelse

         finalflux       = [[finalflux],[bestflux]]
         finalivar       = [[finalivar],[bestivar]]
         finalandmask    = [[finalandmask],[bestandmask]]
         finalormask     = [[finalormask],[bestormask]]
         finaldispersion = [[finaldispersion], [bestdispersion]]
         finalresolution = [[finalresolution], [bestresolution]]
         finalsky        = [[finalsky], [bestsky]]
         mjd_targ        = [mjd_targ,LONG(max(ti.mjds_rm[expidx]))]
         
      endif else begin
         indx = *(indx_tar[itarg])
         sflag = 'Skipping '
         if indx[0] eq -1 then begin
            rflag = ' due to No Matching Fibers'
            rtflag = string('Target ', itarg+1, ' NO DATA')
            nexpf = ''
         endif else begin
            rflag = ' due to No Good Fibermasks'
            nexpf = 'all the ('+strtrim(n_elements(expidx),2)+') '
            if keyword_set(radec_coadd) then begin
                rtflag = string('Target', itarg+1, ' ', plugmap[indx[0]].objtype, $
                                plugmap[indx[0]].mag, format = '(a, i5.4, a, a, f6.2, 5f6.2)')
            endif else begin
                rtflag = string('Catalogid: ', strtrim(plugmap[indx[0]].catalogid,2), ' ', plugmap[indx[0]].objtype, $
                    plugmap[indx[0]].mag, format = '(a, a, a, a, f6.2, 5f6.2)')
            endelse
         endelse
          
         if keyword_set(save_novalid) then begin
            rflag = ''
            sflag = ''
         endif
          
         splog, sflag+ctype+nexpf+c1type+rflag+' ('+strtrim(itarg+1,2)+'/'+strtrim(ntarget,2)+')'
         splog, rtflag
         if indx[0] ne -1 then begin
            if keyword_set(onestep_coadd) then begin
                expidx = indx[0:(n_elements(indx)/2)-1]
                foreach idx, expidx, ie do $
                    splog, 'Exposure '+strtrim(ie+1,2)+' Fibermask: '+sdss_flagname('SPPIXMASK',plugmap[idx].fibermask, /concat)
            endif else begin
                foreach idx, indx, ie do $
                    splog, 'Exposure '+strtrim(ie+1,2)+' Fibermask: '+sdss_flagname('SPPIXMASK',combinedfibermask[idx], /concat)
            endelse
         endif

    
         if not keyword_set(save_novalid) then continue
         ; The code below creates spSpec files, and spField rows for targets with no valid exposures

         itarget = itarget+1
         if indx[0] eq -1 then begin
             thisplug = plugmap[0]
             thisplug = struct_addtags(thisplug,new_cols)
             foreach tag, tag_names(thisplug), itag do begin
                  type = size(thisplug.(itag), /TNAME)  ; Get field type as a string
                  ; Define custom defaults based on type
                  case type of
                     'STRING': thisplug.(itag) = ''
                     'FLOAT': thisplug.(itag) = !VALUES.F_NAN
                     'DOUBLE': thisplug.(itag) = !VALUES.D_NAN
                     'POINTER': thisplug.(itag) = ptr_new(0)  ; Reset to a new null pointer
                     'INT': thisplug.(itag) = -999  ; Example custom int value
                     'UINT': thisplug.(itag) = 0  ; Example custom int value
                     'LONG': thisplug.(itag) = -999L  ; Example custom long value
                     'ULONG': thisplug.(itag) = 0
                     'LONG64': thisplug.(itag) = -999LL
                     'ULONG64': thisplug.(itag) = 0
                     'OBJREF': thisplug.(itag) = OBJ_NEW(0)
                     'BYTE': thisplug.(itag) = 0B
                  endcase
             endforeach
             thisplug[0].fibermask = fibermask_bits('NODATA')
             thisplug[0].exp_disp_med = -1
             thisplug[0].target_index = itarget+1
             thisplug[0].CATALOGID = 'u'+strtrim(itarget+1,2)
         endif else begin
            expidx = indx
            if keyword_set(onestep_coadd) then expidx = indx[0:(n_elements(indx)/2)-1]
            if keyword_set(ra_coadd) then begin
                thisplug = plugmap[expidx[0]]
            endif else thisplug = plugmap[catid_ix[expidx[0]]]

            ;thisplug = plugmap[eindx[0]]
            thisplug = struct_addtags(thisplug,new_cols)
            ffm = 0
            foreach eindx,expidx do $
                ffm = ffm OR combinedfibermask[eindx]

            thisplug = build_fmaprow(itarget, thisplug, expidx, ffm, ti)

            thisplug[0].exp_disp_med = -1
         endelse
         
         finalflux       = [[finalflux],[fltarr(nfinalpix)]]
         finalivar       = [[finalivar],[fltarr(nfinalpix) - 1]]
         finalandmask    = [[finalandmask],[lonarr(nfinalpix)-1]]
         finalormask     = [[finalormask],[lonarr(nfinalpix)-1]]
         finaldispersion = [[finaldispersion], [fltarr(nfinalpix)]]
         finalresolution = [[finalresolution], [fltarr(nfinalpix)]]
         finalsky        = [[finalsky], [fltarr(nfinalpix)]]
         mjd_targ        = [mjd_targ,LONG(max(ti.mjds_rm[itarget]))]

        if not keyword_set(finalplugmap) then finalplugmap = thisplug $
        else finalplugmap = [finalplugmap, thisplug]
      endelse
   endfor
   
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
      platesn, finalflux, finalivar, finalandmask, finalplugmap, finalwave, obs=obs,$
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
   
   master_snr2=dblarr(2,3)
   master_snr2_dered=dblarr(2,3)



   for iexp=0, nexp_tmp - 1 do begin
        platesn, finalflux_rm[*,*,iexp], finalivar_rm[*,*,iexp], $
          finalandmask_rm[*,*,iexp], finalplugmap_rm[*,iexp], finalwave, hdr=bighdr, obs=obs, $
          legacy=legacy, plotfile=djs_filepath(repstr(plotsnfile,'X',string(iexp,format='(i2.2)')), root_dir=combinedir), $
          coeffs=coeffs, snplate=snplate, specsnlimit=specsnlimit, dered_snplate=dered_snplate
        splog, prelog=''
        bands = ['G','R','I']
        for ispec=1, 2 do begin
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

        if keyword_set(legacy) then begin
            for ifib=0, nfiber-1 do begin
                if ifib lt 500 then sp=0 else sp=1
                ti.snr2listG[ifib,iexp]=strtrim(strcompress(string(string(snplate[sp,0],format='(f0.2)'),format='(999a)')),2)
                ti.snr2listR[ifib,iexp]=strtrim(strcompress(string(string(snplate[sp,1],format='(f0.2)'),format='(999a)')),2)
                ti.snr2listI[ifib,iexp]=strtrim(strcompress(string(string(snplate[sp,2],format='(f0.2)'),format='(999a)')),2)
            endfor
            weight[iexp]=snplate[0,2]
        endif else begin
            if strmatch(obs,'APO',/fold_case) eq 1 then sp=0 else sp=1
            for ifib=0, nfiber-1 do begin
                ti.snr2listG[ifib,iexp]=strtrim(strcompress(string(string(snplate[sp,0],format='(f0.2)'),format='(999a)')),2)
                ti.snr2listR[ifib,iexp]=strtrim(strcompress(string(string(snplate[sp,1],format='(f0.2)'),format='(999a)')),2)
                ti.snr2listI[ifib,iexp]=strtrim(strcompress(string(string(snplate[sp,2],format='(f0.2)'),format='(999a)')),2)
            endfor
            weight[iexp]=snplate[sp,2]
        endelse
        if (NOT keyword_set(tailist)) then tailist = tai_t $
        else tailist = [tailist, tai_t]
        
        tai_flag=1
   endfor

   ; Plot S/N and throughput **after** this distortion-correction.
   ; (This over-writes header cards written in the first call.)
   splog, prelog='Final'
   platesn, finalflux, finalivar, finalandmask, finalplugmap, finalwave, obs=obs, $
   hdr=bighdr, legacy=legacy, plotfile=djs_filepath(repstr(plotsnfile,'-X',''), root_dir=combinedir), coeffs=coeffs
   splog, prelog=''
   bands = ['G','R','I']
   for ispec=1, 2 do begin
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
   for ispec=1, 2 do begin
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

   ntarget = n_elements(finalplugmap)

   for itarget=0, ntarget-1 do begin
         tai_t = *(finalplugmap[itarget].tai_pt)
         if not keyword_set(tai_t) then continue
         fiberid_t = *(finalplugmap[itarget].fiberid_pt)
         exp_t = *(finalplugmap[itarget].exp_pt)
         if n_elements(fiberid_t) eq 1 then begin
            match = [fiberid_t-1]
            finalplugmap[itarget].mjd_final = masked_mjd(tai_t/(24.D*3600.D),ti.snr2listI[match])
            finalplugmap[itarget].FIELDSNR2G_LIST = plmap_arr_to_strlist(ti.snr2listG[match],format='(f0.2)')
            finalplugmap[itarget].FIELDSNR2R_LIST = plmap_arr_to_strlist(ti.snr2listR[match],format='(f0.2)')
            finalplugmap[itarget].FIELDSNR2I_LIST = plmap_arr_to_strlist(ti.snr2listI[match],format='(f0.2)')
            finalplugmap[itarget].weight_pt = ptr_new(ti.snr2listI[match])
        endif else begin
            match=[[fiberid_t-1],[exp_t]]
            finalplugmap[itarget].mjd_final = masked_mjd(tai_t/(24.D*3600.D),ti.snr2listI[match[*, 0], match[*, 1]])
            finalplugmap[itarget].FIELDSNR2G_LIST = plmap_arr_to_strlist(ti.snr2listG[match[*, 0], match[*, 1]],format='(f0.2)')
            finalplugmap[itarget].FIELDSNR2R_LIST = plmap_arr_to_strlist(ti.snr2listR[match[*, 0], match[*, 1]],format='(f0.2)')
            finalplugmap[itarget].FIELDSNR2I_LIST = plmap_arr_to_strlist(ti.snr2listI[match[*, 0], match[*, 1]],format='(f0.2)')
            finalplugmap[itarget].weight_pt = ptr_new(ti.snr2listI[match[*, 0], match[*, 1]])
        endelse

        weights_target_f_tmp = *(finalplugmap[itarget].weight_pt)
        weights_target_f_tmp = Double(weights_target_f_tmp)
        airmass_target_f_tmp   = *(finalplugmap[itarget].AIRMASS_pt)
        seeing20_target_f_tmp  = *(finalplugmap[itarget].SEEING20_pt)
        seeing50_target_f_tmp  = *(finalplugmap[itarget].SEEING50_pt)
        seeing80_target_f_tmp  = *(finalplugmap[itarget].SEEING80_pt)
        rmsoff20_target_f_tmp  = *(finalplugmap[itarget].RMSOFF20_pt)
        rmsoff50_target_f_tmp  = *(finalplugmap[itarget].RMSOFF50_pt)
        rmsoff80_target_f_tmp  = *(finalplugmap[itarget].RMSOFF80_pt)
        if total(weights_target_f_tmp,/DOUBLE) le 0 then weights_target_f_tmp = DBLARR(n_elements(weights_target_f_tmp)) +1
        finalplugmap[itarget].airmass = total(airmass_target_f_tmp *weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
        finalplugmap[itarget].seeing20 = total(seeing20_target_f_tmp*weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
        finalplugmap[itarget].seeing50 = total(seeing50_target_f_tmp*weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
        finalplugmap[itarget].seeing80 = total(seeing80_target_f_tmp*weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
        finalplugmap[itarget].RMSOFF20 = total(rmsoff20_target_f_tmp*weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
        finalplugmap[itarget].RMSOFF50 = total(rmsoff50_target_f_tmp*weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
        finalplugmap[itarget].RMSOFF80 = total(rmsoff80_target_f_tmp*weights_target_f_tmp,/DOUBLE)/total(weights_target_f_tmp,/DOUBLE)
   endfor
   
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

      mwrfits_named, flux[*,indx], thisfile, hdr=hdr, name='FLUX', desc=' Flux', /create
      mwrfits_named, fluxivar[*,indx], thisfile, hdr=exthdr, name='IVAR', desc=' Inverse variance'
      mwrfits_named, pixelmask[*,indx], thisfile, hdr=exthdr, name='MASK', desc=' Pixel mask'
      mwrfits_named, wave[*,indx], thisfile, hdr=exthdr, name='WAVELENGTH', desc=' Wavelength solution'
      mwrfits_named, dispersion[*,indx], thisfile, hdr=exthdr, name='WAVEDISP', desc=' Wavelength dispersion'
      
      ;; need a different header for plugmap structure
      ;; undefine it first so that mwrfits doesn't duplicate comments
      ;; on successive writes
      undefine, plughdr
      mwrfits_named, plugmap[indx], thisfile, hdr=plughdr, name='PLUGMAP', desc=' Plugmap structure'
      mwrfits_named, skyflux[*,indx], thisfile, hdr=exthdr, name='SKY', desc=' Subtracted sky flux'
      mwrfits_named, ximg[*,indx], thisfile, hdr=exthdr, name='X', desc=' Trace X locations on CCD'
      mwrfits_named, superflat[*,indx], thisfile, hdr=exthdr, name= 'SUPERFLAT', desc=' Superflat'
      mwrfits_named, resolution[*,indx], thisfile, hdr = exthdr, name='SPRESL', desc=' Spectral resolution'
      
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

   fiber_rollcall, finalandmask, finalwave, legacy=legacy
;   fiber_rollcall, finalormask, finalwave, fibermask = finalplugmap.fibermask, legacy=legacy
;   fiber_rollcall, finalandmask, finalwave, fibermask = finalplugmap.fibermask, legacy=legacy

   ;----------
   ; Remove header cards that were specific to this first exposure
   ; (where we got the header).

   bighdr = delhdrsec(bighdr, 'EXPOSURE INFO')
   bighdr = delhdrsec(bighdr, 'EXPOSURE SETTINGS')
   bighdr = delhdrsec(bighdr, 'SPECTROGRAPH STATUS')
   bighdr = delhdrsec(bighdr, 'APO WEATHER')
   bighdr = delhdrsec(bighdr, 'LCO WEATHER')
   bighdr = delhdrsec(bighdr, 'CALIBRATION SETTINGS')

   Telescope_cards = ['M1ZROT','M1YTRAN','M1XTRAN','M1YTILT','M1XTILT','M1PISTON',$
                      'M2ZROT','M2YTRAN','M2XTRAN','M2YTILT','M2XTILT','M2PISTON',$
                      'GUIDOFFR','GUIDOFFY','GUIDOFFX','CALOFFR','CALOFFY','CALOFFX',$
                      'ARCOFFY','ARCOFFX','BOREOFFY','BOREOFFX','OBJSYS','SCALE',$
                      'FOCUS','ROTPOS', 'HA','GSEEING']
   sxdelpar, bighdr, Telescope_cards

   sxaddpar, bighdr, 'FIELDID', (strsplit(outputname,'-',/EXTRACT))[1]
   sxaddpar, bighdr, 'FIELDCAD', Field_cadence, after='FIELDID'
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
   for i=0, 4 do sxdelpar, bighdr, 'GAIN'+strtrim(string(i),2)
   for i=0, 4 do sxdelpar, bighdr, 'RDNOISE'+strtrim(string(i),2)
   sxdelpar, bighdr, ['CAMCOL', 'CAMROW']
   sxdelpar, bighdr, ['AMPLL', 'AMPLR', 'AMPUL', 'AMPUR']
   sxdelpar, bighdr, ['FFS', 'FF', 'NE', 'HGCD','HEAR']
   sxdelpar, bighdr, ['SPEC1', 'SPEC2']
   sxdelpar, bighdr, 'NBLEAD'
;   sxdelpar, bighdr, 'PIXFLAT'
;   sxdelpar, bighdr, 'PIXBIAS'
   sxdelpar, bighdr, 'FLATFILE'
   sxdelpar, bighdr, 'ARCFILE'
   sxdelpar, bighdr, 'OBJFILE'
   sxdelpar, bighdr, 'FRAMESN2'
   sxdelpar, bighdr, 'DEREDSN2'

   ;----------
   ; Average together some of the fields from the individual headers. fieldid

   cardnames = [ 'AZ', 'ALT', 'AIRMASS', 'TAI', 'WTIME', 'AIRTEMP', 'DEWPOINT', $
    'DEWDEP', 'DUSTA', 'DUSTB', 'DUSTC', 'DUSTD', 'GUSTS', 'HUMIDITY', $
    'HUMIDOUT', 'PRESSURE', 'WINDD', 'WINDS', 'TEMP01', 'TEMP02', $
    'TEMP03', 'TEMP04', 'HELIO_RV', 'V_RAD', 'SEEING20', 'SEEING50', 'SEEING80', $
    'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
    'WSIGMA', 'XSIGMA', 'CONFSFIL' ]
   sxdelpar, bighdr,cardnames
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
       sxdelpar, bighdr,cardnames
   endfor

   ;----------
   ; Use the MJD passed as a keyword, which will typically be for the most
   ; observation, and be consistent with the output file names

   if (keyword_set(mjd)) then $
    sxaddpar, bighdr, 'MJD', mjd, after = 'FIELDID'

   ; Get the list of MJD's used for these reductions, then convert to a string
   mjdlist = mjdlist[uniq(mjdlist, sort(mjdlist))]
   mjdlist = strtrim(strcompress(string(mjdlist,format='(999a)')),2)

   ; Get the list of Designs used for these reductions, then convert to a string
   designlist = designlist[uniq(designlist, sort(designlist))]
   designlist = strtrim(strcompress(string(designlist,format='(999a)')),2)
 
   ; Get the list of configurations used for these reductions, then convert to a string
   configlist = configlist[uniq(configlist, sort(configlist))]
   configlist = strtrim(strcompress(strjoin(configlist, ' ')),2)

   stailist = tailist[uniq(tailist, sort(tailist))]
   stailist = strtrim(strcompress(string(string(tailist,format='(i15)'),format='(999a)')),2)
   if keyword_set(tai_flag) then begin
    indtai=uniq(tailist, sort(tailist))
    tailist = tailist[indtai]
    snr2listG = ti.snr2listG[indtai]
    snr2listR = ti.snr2listR[indtai]
    snr2listI = ti.snr2listI[indtai]
    tailist = strtrim(strcompress(string(string(tailist,format='(i15)'),format='(999a)')),2)
    snr2listG = strtrim(strcompress(string(snr2listG,format='(999a)')),2)
    snr2listR = strtrim(strcompress(string(snr2listR,format='(999a)')),2)
    snr2listI = strtrim(strcompress(string(snr2listI,format='(999a)')),2)
   endif
   ;----------
   ; Add new header cards

   sxaddpar, bighdr, 'VERSCOMB', idlspec2d_version(), $
    ' Version of idlspec2d for combining multiple spectra', after='VERS2D'
    
   bighdr = addhdrsec(bighdr,"Exposure Info", 'NEXP',nfiles, $
                     ' Number of exposures in this file', before="FIELD/PLATE INFO")
   
    sxaddpar, bighdr, 'EXPTIME', min(exptimevec), $
    ' Minimum of exposure times for all cameras', after='NEXP'
    
;   sxaddpar, bighdr, 'NEXP', nfiles, $
;    ' Number of exposures in this file', before='EXPTIME'
   for ifile=0,nfiles-1 do $
    sxaddpar, bighdr, string('EXPID',ifile+1, format='(a5,i2.2)'), label[ifile], $
     ' ID string for exposure '+strtrim(ifile+1,2), before='EXPTIME'
   if (keyword_set(bestexpnum)) then sxaddpar, bighdr, 'BESTEXP', bestexpnum, before='EXPID01'


    
   foreach cam, ['B1','R1','B2','R2'], idx do begin
    icam = where(strmatch(camnames, cam, /fold_case),ct)
    if ct eq 0 then nexpcam = 0 else nexpcam = nexpvec[icam]
    sxaddpar, bighdr, 'NEXP_'+cam, nexpcam, $
     ' '+cam+' camera number of exposures', before='EXPTIME'
   endforeach
   foreach cam, ['B1','R1','B2','R2'], idx do begin
    icam = where(strmatch(camnames, cam, /fold_case),ct)
    if ct eq 0 then exptimecam = 0.0 else exptimecam = exptimevec[icam]
    sxaddpar, bighdr, 'EXPT_'+cam, exptimecam, $
     ' '+cam+' camera exposure time (seconds)', before='EXPTIME'
   endforeach

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
    fbadpix = fbadpix2 $
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
                    'ALPHA', 'BETA', 'FIBERID',$
                    'tai_pt','fiberid_pt','exp_pt', 'weight_pt','AIRMASS_pt',$
                    'SEEING20_pt','SEEING50_pt','SEEING80_pt',$
                    'RMSOFF20_pt','RMSOFF50_pt','RMSOFF80_pt' ]
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
   ; Add Condition Headers
    fieldhdr = bighdr

    cardnames_avg = ['SEEING20', 'SEEING50', 'SEEING80', 'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
                         'WSIGMA', 'XSIGMA','AIRTEMP','AIRMASS', 'TAI']

    sxaddpar, fieldhdr, 'MJDLIST', mjdlist, key_match_dict['MJDLIST'], after='EXPTIME'
    ;sxaddpar, fieldhdr, 'TAILIST', stailist, key_match_dict['TAILIST']
    sxaddpar, fieldhdr, 'DESIGNS', designlist, key_match_dict['DESIGNS'], after='MJDLIST'
    sxaddpar, fieldhdr, 'CONFIGS', configlist, key_match_dict['CONFIGS'], after='DESIGNS'

    foreach cardname, cardnames_avg, ci do begin
        if ci eq 0 then after = 'DESIGNS' else after = cardnames_avg[ci-1]
        sxcombinepar_v2, hdrarr, cardname, fieldhdr, Comment=key_match_dict[cardname],$
                         func='average',camnames=camnames, AFTER=after, /SaveComment
    endforeach
    sxcombinepar_v2, hdrarr, 'TAI-BEG', fieldhdr, Comment=key_match_dict['TAI'], func='average', outcard='TAI', after='CONFIGS'
    sxcombinepar_v2, hdrarr, 'TAI-BEG', fieldhdr, Comment=key_match_dict['TAIBEG'], func='min',  after='TAI'
    sxcombinepar_v2, hdrarr, 'TAI-END', fieldhdr, Comment=key_match_dict['TAIEND'], func='max',  after='TAI-BEG'

    sxaddpar, fieldhdr, 'DATE-OBS', sxpar(*hdrarr[0],'DATE-OBS'),$
            key_match_dict['DATEOBS'], before='TAI'
    ftai = sxpar(fieldhdr, 'TAI')
    jdtemp=ftai/(24.D*3600.D)
    jdtemp=jdtemp+2400000.5
    mphase,jdtemp,mfrac
    sxaddpar, fieldhdr, 'MOONFRAC',mfrac, 'Moon Phase', after='EXPTIME'


    sxcombinepar_v2, hdrarr, 'XCHI2', fieldhdr, Comment=key_match_dict['XCHI2MAX'], func='max', outcard='XCHI2MAX', after='XCHI2'
    sxcombinepar_v2, hdrarr, 'XCHI2', fieldhdr, Comment=key_match_dict['XCHI2MIN'], func='min', outcard='XCHI2MIN', after='XCHI2'
    sxcombinepar_v2, hdrarr, 'SKYCHI2', fieldhdr, Comment=key_match_dict['SCHI2MAX'], func='max', outcard='SCHI2MAX', after='SKYCHI2'
    sxcombinepar_v2, hdrarr, 'SKYCHI2', fieldhdr, Comment=key_match_dict['SCHI2MIN'], func='min', outcard='SCHI2MIN', after='SKYCHI2'
    sxcombinepar_v2, hdrarr, 'WSIGMA', fieldhdr, Comment=key_match_dict['WSIGMAX'], func='max', outcard='WSIGMAX', after='WSIGMA'
    sxcombinepar_v2, hdrarr, 'WSIGMA', fieldhdr, Comment=key_match_dict['WSIGMIN'], func='min', outcard='WSIGMIN', after='WSIGMA'
    sxcombinepar_v2, hdrarr, 'XSIGMA', fieldhdr, Comment=key_match_dict['XSIGMAX'], func='max', outcard='XSIGMAX', after='XSIGMA'
    sxcombinepar_v2, hdrarr, 'XSIGMA', fieldhdr, Comment=key_match_dict['XSIGMIN'], func='min', outcard='XSIGMIN', after='XSIGMA'
    sxcombinepar_v2, hdrarr, 'NGUIDE', fieldhdr, Comment=key_match_dict['NGUIDE'], func='total', camnames=camnames, after='MJD'
    sxcombinepar_v2, hdrarr, 'EXPTIME', fieldhdr, Comment=key_match_dict['EXPTIME'], func='total', camnames=camnames
    sxaddpar, fieldhdr, 'NEXP', n_elements(hdrarr)/n_elements(camnames), key_match_dict['NEXP']
   ;---------------------------------------------------------------------------
   ; Write combined output file
   ;---------------------------------------------------------------------------

   if keyword_set(onestep_coadd) then begin
        ; First write the file with the flux distortion vectors
        if not keyword_set(nodist) then $
            mwrfits_named, corrimg, distortfitsfile, hdr=fieldhdr, name= 'CORRIMG', /create $
        else mwrfits_named, corrimg, distortfitsfile, hdr=fieldhdr,  /create
   endif

   fulloutname = djs_filepath(outputname, root_dir=combinedir)

   ; HDU #0 is flux
   sxaddpar, fieldhdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
   mwrfits_named, finalflux, fulloutname, hdr=fieldhdr, name='FLUX', /create

   ; HDU #1 is inverse variance
   sxaddpar, hdrfloat, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
   mwrfits_named, finalivar, fulloutname, hdr=hdrfloat, name='IVAR', desc=' Inverse variance'

   ; HDU #2 is AND-pixelmask
   mwrfits_named, finalandmask, fulloutname, hdr=hdrlong, name='ANDMASK', desc=' AND Mask'

   ; HDU #3 is OR-pixelmask
   mwrfits_named, finalormask, fulloutname, hdr=hdrlong, name='ORMASK', desc=' OR Mask'

   ; HDU #4 is dispersion map
   sxaddpar, hdrfloat, 'BUNIT', 'pixels'
   mwrfits_named, finaldispersion, fulloutname, hdr=hdrfloat,name='WAVEDISP', desc=' Wavelength dispersion'

   ; HDU #5 is plugmap
   mwrfits_named, finalplugmap, fulloutname, hdr=hdrplug, name='PLUGMAP', desc=' Plugmap structure'

   ; HDU #6 is the sky
   mwrfits_named, finalsky, fulloutname, hdr=hdrsky, name='SKY', desc=' Subtracted sky flux'
   
   ; HDU #7 is the resolution map
   sxaddpar, hdrfloat, 'BUNIT', 'angstroms'
   mwrfits_named, finalresolution, fulloutname, hdr = hdrfloat, name ='SPECRESL', desc=' Spectral resolution'
   
   if keyword_set(single_spectra) then begin
    ;writing each individual coadd spectrum on the field
    sxdelpar, bighdr, 'NAXIS2'
    spawn,'mkdir -p '+combinedir+'coadd'
    for itarget=0, ntarget-1 do begin
       added_exp=[]

       finalvalues_t = {flux:0.0, loglam: 0.0, ivar: 0.0, $
                        and_mask:long(0), or_mask:long(0),$
                        wdisp:0.0, sky:0.0, wresl:0.0}

       finalvalues = replicate(finalvalues_t, n_elements(finalwave))
       finalvalues.flux=finalflux[*,itarget]
       finalvalues.loglam=finalwave
       finalvalues.ivar=finalivar[*,itarget]
       finalvalues.and_mask=finalandmask[*,itarget]
       finalvalues.or_mask=finalormask[*,itarget]
       finalvalues.wdisp=finaldispersion[*,itarget]
       finalvalues.sky=finalsky[*,itarget]
       finalvalues.wresl=finalresolution[*,itarget]
       
       if keyword_set(legacy) then begin
          targid_tar=string(finalplugmap[itarget].TARGET_INDEX,format='(i4.4)')
       endif else begin
          targid_tar=finalplugmap[itarget].catalogid
       endelse
       sxaddpar, bighdr, 'PLUG_RA', finalplugmap[itarget].ra, ' RA of Target'
       sxaddpar, bighdr, 'PLUG_DEC', finalplugmap[itarget].dec, ' DEC of Target'
       if keyword_set(legacy) or keyword_set(plates) then begin
         if keyword_set(mjd) then begin
           thismjd=mjd
         endif else begin
           thismjd=mjd_targ[itarget]
         endelse
       endif else begin
         if keyword_set(epoch) then begin
           thismjd = mjd
         endif else begin
           thismjd=mjd_targ[itarget]
         endelse
       endelse
       thismjd=strtrim(strcompress(string(thismjd,format='(99a)')),2)
       coadddir=combinedir+'coadd/'+thismjd
       spawn,'mkdir -p '+coadddir
       coaddname = repstr(repstr(outputname,'spField','spSpec'),'.fits', $
        '-'+strtrim(targid_tar,2)+'.fits') 
       fulloutname_coadd = djs_filepath(coaddname, root_dir=coadddir)
       ; HDU # 0 header
       mwrfits_named, junk_d, fulloutname_coadd, hdr=bighdr, /create
       ; HDU # 1 header
       mwrfits_named, finalvalues, fulloutname_coadd, hdr=coadd_val, name='COADD', desc=' Coadded spectrum'

       ; HDU #2 is plugmap
       mwrfits_named, finalplugmap[itarget], fulloutname_coadd, hdr = hdrplug, name='PLUGMAP', desc= ' Plugmap structure', /silent

       for ifiber=0, nfiber-1 do begin
        for iexp=0, nexp_tmp - 1 do begin
         if keyword_set(legacy) then begin
           targid_rm=string(finalplugmap_rm[ifiber,iexp].fiberid,format='(i4.4)');targid_tar
         endif else begin
           targid_rm=finalplugmap_rm[ifiber,iexp].catalogid
         endelse
         if targid_rm eq targid_tar then begin
           added_exp=[added_exp,iexp]
           
           finalvalues_rm = replicate(finalvalues_t, n_elements(finalwave))
           finalvalues_rm.flux=finalflux_rm[*,ifiber,iexp]
           finalvalues_rm.loglam=finalwave
           finalvalues_rm.ivar=finalivar_rm[*,ifiber,iexp]
           finalvalues_rm.and_mask=finalandmask_rm[*,ifiber,iexp]
           finalvalues_rm.or_mask=finalormask_rm[*,ifiber,iexp]
           finalvalues_rm.wdisp=finaldispersion_rm[*,ifiber,iexp]
           finalvalues_rm.sky=finalsky_rm[*,ifiber,iexp]
           finalvalues_rm.wresl=finalresolution_rm[*,ifiber,iexp]
           
           ; HDU # N header
           thisconf=ti.mjds_rm_summ[ifiber,iexp];  config_rm[ifiber,iexp]
           thisconf=string(thisconf,format='(i5.5)')+'-'+string(iexp,format='(i2.2)')
           
           indv_val=*hdrarr[iexp]
           sxdelpar,indv_val, ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND','CAMERAS']
           mwrfits_named, finalvalues_rm, fulloutname_coadd, hdr=indv_val, name='MJD_EXP_'+thisconf, desc=' Single exposure spectrum', /silent
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
                         'TEMP03', 'TEMP04', 'HELIO_RV', 'V_RAD', 'SEEING20', 'SEEING50', 'SEEING80', $
                         'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
                         'WSIGMA', 'XSIGMA' , 'CCDTEMP', 'LN2TEMP']
  
        h = headfits(fulloutname_coadd) ;Read primary header
    
        if n_elements(added_exp) eq 1 then begin
          sxaddpar, h, 'TAI',     SXPAR(*iused_hdrarr[0], 'TAI-BEG'),  key_match_dict['TAI'],    after='EXPTIME'

          foreach card, cardnames_avg, ic do begin
            val= SXPAR(*iused_hdrarr[0], card, COMMENT = cmt)
            sxaddpar, h, card, val, cmt, /SaveComment
          endforeach
          sxaddpar, h, 'EXPTIME',  SXPAR(*iused_hdrarr[0], 'EXPTIME'), key_match_dict['EXPTIME']
          sxaddpar, h, 'NEXP',    n_elements(added_exp),               key_match_dict['NEXP']
          sxaddpar, h, 'TAI-BEG', SXPAR(*iused_hdrarr[0], 'TAI-BEG'),  key_match_dict['TAIBEG'], after='TAI'
          sxaddpar, h, 'TAI-END', SXPAR(*iused_hdrarr[0], 'TAI-END'),  key_match_dict['TAIEND'], after='TAI-BEG'
          sxaddpar, h, 'XCHI2MAX', SXPAR(*iused_hdrarr[0], 'XCHI2'),   key_match_dict['XCHI2MAX'], after='XCHI2'
          sxaddpar, h, 'XCHI2MIN', SXPAR(*iused_hdrarr[0], 'XCHI2'),   key_match_dict['XCHI2MIN'], after='XCHI2'
          sxaddpar, h, 'SCHI2MAX', SXPAR(*iused_hdrarr[0], 'SKYCHI2'), key_match_dict['SCHI2MAX'], after='SKYCHI2'
          sxaddpar, h, 'SCHI2MIN', SXPAR(*iused_hdrarr[0], 'SKYCHI2'), key_match_dict['SCHI2MIN'], after='SKYCHI2'
          sxaddpar, h, 'WSIGMAX',  SXPAR(*iused_hdrarr[0], 'WSIGMA'),  key_match_dict['WSIGMAX'], after='WSIGMA'
          sxaddpar, h, 'WSIGMIN',  SXPAR(*iused_hdrarr[0], 'WSIGMA'),  key_match_dict['WSIGMIN'], after='WSIGMA'
          sxaddpar, h, 'XSIGMAX',  SXPAR(*iused_hdrarr[0], 'XSIGMA'),  key_match_dict['XSIGMAX'], after='XSIGMA'
          sxaddpar, h, 'XSIGMIN',  SXPAR(*iused_hdrarr[0], 'XSIGMA'),  key_match_dict['XSIGMIN'], after='XSIGMA'
          sxaddpar, h, 'NGUIDE',   SXPAR(*iused_hdrarr[0], 'NGUIDE'),  key_match_dict['NGUIDE'],  before='FIELDCAD'
          sxdelpar, h, ['SHOPETIM', 'SHCLOTIM', 'ionpump']
          
          
          
        endif else begin
           sxcombinepar_v2, iused_hdrarr, 'TAI-BEG', h, Comment=key_match_dict['TAI'], func='average', outcard='TAI', weights=used_weight,after='EXPTIME'
           foreach cardname, cardnames_avg do begin
               sxcombinepar_v2, iused_hdrarr, cardname, h, Comment=key_match_dict[cardname], func='average', weights=used_weight, /SaveComment
           endforeach

           sxcombinepar_v2, iused_hdrarr, 'EXPTIME', h, Comment=key_match_dict['EXPTIME'], func='total'
           sxaddpar, h, 'NEXP', n_elements(added_exp), key_match_dict['NEXP']
           sxcombinepar_v2, iused_hdrarr, 'TAI-BEG', h, Comment=key_match_dict['TAIBEG'], func='min', after='TAI'
           sxcombinepar_v2, iused_hdrarr, 'TAI-END', h, Comment=key_match_dict['TAIEND'], func='max', after='TAI-BEG'
           
           sxcombinepar_v2, iused_hdrarr, 'XCHI2', h, Comment=key_match_dict['XCHI2MAX'], func='max', outcard='XCHI2MAX', after='XCHI2'
           sxcombinepar_v2, iused_hdrarr, 'XCHI2', h, Comment=key_match_dict['XCHI2MIN'], func='min', outcard='XCHI2MIN', after='XCHI2'
           sxcombinepar_v2, iused_hdrarr, 'SKYCHI2', h, Comment=key_match_dict['SCHI2MAX'], func='max', outcard='SCHI2MAX', after='SKYCHI2'
           sxcombinepar_v2, iused_hdrarr, 'SKYCHI2', h, Comment=key_match_dict['SCHI2MIN'], func='min', outcard='SCHI2MIN', after='SKYCHI2'
           sxcombinepar_v2, iused_hdrarr, 'WSIGMA', h, Comment=key_match_dict['WSIGMAX'], func='max', outcard='WSIGMAX', after='WSIGMA'
           sxcombinepar_v2, iused_hdrarr, 'WSIGMA', h, Comment=key_match_dict['WSIGMIN'], func='min', outcard='WSIGMIN', after='WSIGMA'
           sxcombinepar_v2, iused_hdrarr, 'XSIGMA', h, Comment=key_match_dict['XSIGMAX'], func='max', outcard='XSIGMAX', after='XSIGMA'
           sxcombinepar_v2, iused_hdrarr, 'XSIGMA', h, Comment=key_match_dict['XSIGMIN'], func='min', outcard='XSIGMIN', after='XSIGMA'
           sxcombinepar_v2, iused_hdrarr, 'NGUIDE', h, Comment=key_match_dict['NGUIDE'], func='total',  before='FIELDCAD'
           sxdelpar, h, ['SHOPETIM', 'SHCLOTIM', 'ionpump']

           
        endelse
        modfits,fulloutname_coadd,0,h ;Update header only
      endif

    endfor
   endif
   return
end
;------------------------------------------------------------------------------
