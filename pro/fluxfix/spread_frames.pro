;+
; NAME:
;   spread_frames
;
; PURPOSE:
;   Read in multiple frame files produced by spectro2d and pass back arrays
;   of the various quantities.
;
; CALLING SEQUENCE:
;   spread_frames, spframes, [window=, binsz =,  adderr=, camnames=, $
;     tsobjname=, flux=, ivar=, loglam=, dispersion=,  pixelmask=, plugmap=, $
;     plugtag=, camerasvec=, filenum=,  expid=, sn2=, hdrarr=, merged_hdr=]
;
; INPUTS:
;   spframes -
;
; OPTIONAL INPUT KEYWORDS:
;   window   -  window size for apodizing the errors of the spectrum from each
;               individual frame -- default to 100 pixels on each end of the 
;               spectrum
;   binsz    -  bin size (in log10(ang) of the output spectra; default to 
;               1d-4 which is 69 km/s
;   adderr   -  additional error to add to the formal errors as a fraction 
;               of the flux
;   camnames  - camera names to combine -- default to ['b1', 'b2', 'r1', 'r2']
;  
;   tsobjname - full path name to tsobjfile if one is available. It is 
;               assumed that the target info is in HDU#1 and the latest 
;               photo re-run is in HDU#2
;
; OUTPUT:
;
; OPTIONAL OUTPUT:
;   flux       - flux array from all frames [npix, nfiber*nfames] 
;   ivar       - inverse variance array from all frames [npix, nfiber*nframes] 
;   loglam     - wavelength array (log10(Ang)) [npix, nfiber*nframes] 
;   dispersion - instrumental resolution array [npix, nfiber*nfrmes]
;   pixelmask  - mask array [npix, nfiber*nframes]
;   plugmap    - array of structures containing plugmap [nfiber*nframes]
;   plugtag    - like plugmap but with newer tsObj info if available and with
;                exposure ID info included (useful for bookkeeping)
;   camerasvec - vector of camera IDs (b1,b2,r1,r2) [nframes]
;   filenum    - file number of frame (used in bookkeeping) [nfiber*nframes]
;   expid      - exposure ID number [nframes]  
;   sn2        - signal-to-noise squared of each frame [nframes]
;   hdrarr     - array of pointers to header files from each frame [nframes]
;   merged_hdr - header to be used for combined frames (some quantities 
;                keywords contain averages of the keywords in the frames)
;
; COMMENTS:
;
; BUGS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   correct_dlam
;   djs_diff_angle
;   mrdfits
;   pixelmask_bits()
;   splog
;   sxaddpar
;   sxcombinepar
;   sxdelpar
;   traceset2xy
;   idlspec2d_version
;   
; INTERNAL SUPPORT ROUTINES:
;   makelabel()
;
; REVISION HISTORY:
;   12-Aug-2003  Made into a stand-alone routine by C. Tremonti, Steward Obs.
;   02-Jan-2000  Written as part of "spcoadd_frames" by D. Schlegel, Princeton
;
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

;-----------------------------------------------------------------------------

pro spread_frames, spframes, window=window, binsz = binsz, $
    adderr=adderr, camnames=camnames, tsobjname = tsobjname, $
    flux = flux, ivar = ivar, loglam = loglam, dispersion = dispersion, $
    pixelmask = pixelmask, plugmap = plugmap, plugtag = plugtag, $
    camerasvec = camerasvec, filenum = filenum, expid = expid, sn2 = sn2, $
    hdrarr = hdrarr, merged_hdr = merged_hdr
 
   ;---------------------------------------------------------------------------

   if (n_elements(window) EQ 0) then window = 100

   nfiles = n_elements(spframes)
   if (nfiles EQ 0) then return

   if NOT keyword_set(camnames) then camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)
   exptimevec = fltarr(ncam)

   plugtag_struct = {plateid: 0, mjd: 0, fiberid: 0, $
                     ra: 0.0, dec: 0.0, mag: fltarr(5), $
                     xfocal: 0.0, yfocal: 0.0, objtype: ' ', $
                     camcolor: ' ', spectrographid: 0, expid: ' ', $
                     tsobjid: lonarr(5)}

   ;-----------------
   ; Read in tsObjfile -- assume target info in HDU#1, new info in HDU#2
  
   if keyword_set(tsobjname) then tsobj = mrdfits(tsobjname, 2)

   ;---------------------------------------------------------------------------
   ; Loop through each 2D output and read in the data
   ;---------------------------------------------------------------------------

   for ifile=0, nfiles-1 do begin

      ;----------
      ; Read in all data from this input file.
      ; Reading the plug-map structure will fail if its structure is
      ; different between different files.

      splog, 'Reading file #', ifile, ': ', spframes[ifile]
      tempflux = mrdfits(spframes[ifile], 0, hdr)
      tempivar = mrdfits(spframes[ifile], 1)
      temppixmask = mrdfits(spframes[ifile], 2)

      ;  Zero out the following four mask bits
      bitval  = pixelmask_bits('COMBINEREJ') + pixelmask_bits('SMEARIMAGE') + $
                pixelmask_bits('SMEARHIGHSN') + pixelmask_bits('SMEARMEDSN')
      
      temppixmask = temppixmask AND (NOT bitval) 

      tempwset = mrdfits(spframes[ifile], 3)
      tempdispset = mrdfits(spframes[ifile], 4)
      tempplug = mrdfits(spframes[ifile], 5, structyp='PLUGMAPOBJ')
      if (NOT keyword_set(tempflux)) then $
       message, 'Error reading file ' + spframes[ifile]

      if (ifile EQ 0) then $
       hdrarr = ptr_new(hdr) $
      else $
       hdrarr = [hdrarr, ptr_new(hdr)]

      thismjd = sxpar(hdr, 'MJD')
      if (NOT keyword_set(mjdlist)) then mjdlist = thismjd $
       else mjdlist = [mjdlist, thismjd]

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
      expstr = string(sxpar(hdr, 'EXPOSURE'), format='(i8.8)')
      framesn2 = sxpar(hdr, 'FRAMESN2')

      ;----------
      ; Solve for wavelength and lambda-dispersion at each pixel in the image

      traceset2xy, tempwset, junk, temploglam
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

      ;----------
      ; Determine if this is a blue or red spectrum

      icam = (where(cameras EQ camnames))[0]
      if (icam EQ -1) then $
       message, 'Unknown camera ' + cameras
      exptimevec[icam] = exptimevec[icam] + sxpar(hdr, 'EXPTIME')

      ;----------
      ; Apodize the errors
      ; Do this only for the dichroic overlap region, which are the first
      ; rows in both the blue and red CCD's.

      if (keyword_set(window)) then begin
         swin = window < npix
         indx = lindgen(swin)
         tempivar[indx,*] = tempivar[indx,*] * (indx # replicate(1,nfib)) / swin
      endif
     
      ;------------------------
      ; Build structure like plugmap to carry all important tags

      tempplugtag = make_array(dim = nfib, val = plugtag_struct)
      struct_assign, tempplug, tempplugtag
      tempplugtag.plateid = 0  ; Fill in later!
      tempplugtag.mjd = 0  ; Fill in later!
      tempplugtag[*].camcolor = strmid(cameras, 0, 1)
      tempplugtag.expid = expstr
      ; spectrograph ID (1/2) is reset b/c it is -1 for unmapped fiberse
      tempplugtag[*].spectrographid = strmid(cameras, 1, 1)

      ;---------------
      ; Match tsObj to plugmap and update plugtag structure with fibermags
      ; from the tsObj

      if keyword_set(tsobjname) then begin
        tsobjid = [[tsobj.run], [tsobj.rerun], [tsobj.camcol], $
                   [tsobj.field], [tsobj.id]]

        for ifib = 0, nfib - 1 do begin
          adist = djs_diff_angle(tsobj.ra, tsobj.dec, $
                  tempplugtag[ifib].ra, tempplugtag[ifib].dec, $
                  units='degrees')
          match = where(adist LT 2./3600)

          ; No match will be found for unmapped fibers so set mags to zero
          if match[0] ne -1 then begin
             tempplugtag[ifib].mag = tsobj[match].fibercounts 
             tempplugtag[ifib].tsobjid = tsobjid[match,*]
          endif else begin
             splog, 'Fiber mags set to zero for unmapped fiber: ', $
                     tempplugtag[ifib].fiberid
             tempplugtag[ifib].mag = [0,0,0,0,0]
          endelse
        endfor
      endif

      ;----------
      ; Concatenate data from all images

      if (ifile EQ 0) then begin
         flux = tempflux
         ivar = tempivar
         loglam = temploglam
         dispersion = tempdispersion
         pixelmask = temppixmask

         camerasvec = cameras
         label = makelabel(hdr)
         filenum = lonarr(nfib) + ifile
         plugmap = tempplug
         expid = expstr 
         sn2 = framesn2
         plugtag = tempplugtag
      endif else begin
         ; Append as images...
         flux = [[flux], [tempflux]]
         ivar = [[ivar], [tempivar]]
         loglam = [[loglam], [temploglam]]
         dispersion = [[dispersion], [tempdispersion]]
         pixelmask = [[pixelmask], [temppixmask]]

         ; Append as vectors...
         camerasvec = [camerasvec, cameras]
         label = [label, makelabel(hdr)]
         filenum = [filenum, lonarr(nfib) + ifile]
         plugmap = [plugmap, tempplug]
         expid = [expid, expstr] 
         sn2 = [sn2, framesn2]
         plugtag = [plugtag, tempplugtag]
      endelse
   endfor

   ;---------------------------------------------------------------------------
   ; Create header for combined frames
   ;---------------------------------------------------------------------------

   ; Modify the 1st file's header to use for the combined plate header.
  
   merged_hdr = *hdrarr[0]

   ;----------
   ; Remove header cards that were specific to this first exposure
   ; (where we got the header).

   ncoeff = sxpar(merged_hdr, 'NWORDER')
   for i=2, ncoeff-1 do sxdelpar, merged_hdr, 'COEFF'+strtrim(string(i),2)

   sxdelpar, merged_hdr, ['SPA', 'IPA', 'IPARATE']
   sxdelpar, merged_hdr, 'EXPOSURE'
   sxdelpar, merged_hdr, 'SEQID'
   sxdelpar, merged_hdr, 'DARKTIME'
   sxdelpar, merged_hdr, 'CAMERAS'
   sxdelpar, merged_hdr, 'PLUGMAPO'
   for i=1, 4 do sxdelpar, merged_hdr, 'GAIN'+strtrim(string(i),2)
   for i=1, 4 do sxdelpar, merged_hdr, 'RDNOISE'+strtrim(string(i),2)
   sxdelpar, merged_hdr, ['CAMCOL', 'CAMROW']
   sxdelpar, merged_hdr, ['AMPLL', 'AMPLR', 'AMPUL', 'AMPUR']
   sxdelpar, merged_hdr, ['FFS', 'FF', 'NE', 'HGCD']
   sxdelpar, merged_hdr, ['SPEC1', 'SPEC2']
   sxdelpar, merged_hdr, 'NBLEAD'
   sxdelpar, merged_hdr, 'PIXFLAT'
   sxdelpar, merged_hdr, 'PIXBIAS'
   sxdelpar, merged_hdr, 'FLATFILE'
   sxdelpar, merged_hdr, 'ARCFILE'
   sxdelpar, merged_hdr, 'OBJFILE'
   sxdelpar, merged_hdr, 'FRAMESN2'

   ;----------
   ; Average together some of the fields from the individual headers.
   ; CT -- Weight by S/N^2 since this is the effective weighting of the
   ; exposures when they are combined

   cardname = [ 'AZ', 'ALT', 'TAI', 'WTIME', 'AIRTEMP', 'DEWPOINT', $
    'DEWDEP', 'DUSTA', 'DUSTB', 'DUSTC', 'DUSTD', 'GUSTS', 'HUMIDITY', $
    'HUMIDOUT', 'PRESSURE', 'WINDD', 'WINDS', 'TEMP01', 'TEMP02', $
    'TEMP03', 'TEMP04', 'HELIO_RV', 'SEEING20', 'SEEING50', 'SEEING80', $
    'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
    'WSIGMA', 'XSIGMA', 'AIRMASS']

   sxcombinepar, hdrarr, cardname, merged_hdr, func='average', weights=sn2

   sxcombinepar, hdrarr, 'TAI-BEG', merged_hdr, func='min'
   sxcombinepar, hdrarr, 'TAI-END', merged_hdr, func='max'

   sxcombinepar, hdrarr, 'XCHI2', merged_hdr, func='max', outcard='XCHI2MAX', $
    after='XCHI2'
   sxcombinepar, hdrarr, 'XCHI2', merged_hdr, func='min', outcard='XCHI2MIN', $
    after='XCHI2'

   sxcombinepar, hdrarr, 'SKYCHI2', merged_hdr, func='max', $
     outcard='SCHI2MAX', after='SKYCHI2'
   sxcombinepar, hdrarr, 'SKYCHI2', merged_hdr, func='min', $
     outcard='SCHI2MIN', after='SKYCHI2'

   sxcombinepar, hdrarr, 'WSIGMA', merged_hdr, func='max', outcard='WSIGMAX', $
    after='WSIGMA'
   sxcombinepar, hdrarr, 'WSIGMA', merged_hdr, func='min', outcard='WSIGMIN', $
    after='WSIGMA'

   sxcombinepar, hdrarr, 'XSIGMA', merged_hdr, func='max', outcard='XSIGMAX', $
    after='XSIGMA'
   sxcombinepar, hdrarr, 'XSIGMA', merged_hdr, func='min', outcard='XSIGMIN', $
    after='XSIGMA'

   ; Add the NGUIDE keywords for all headers of one flavor of CAMERAS
   ; (e.g., for all the 'b1' exposures if the first frame is 'b1'.)

   cardname = 'NGUIDE'
   sxcombinepar, hdrarr[0], cardname, merged_hdr, func='total'
   cameras0 = sxpar(*(hdrarr[0]), 'CAMERAS')
   for ihdr=1, n_elements(hdrarr)-1 do begin
      if (sxpar(*(hdrarr[ihdr]), 'CAMERAS') EQ cameras0) then $
       sxcombinepar, hdrarr[ihdr], cardname, merged_hdr, func='total'
   endfor

   ;----------
   ; Get the list of MJD's used for these reductions, then convert to a string
   ; The header keyword MJD will be that of the first exposure (this can
   ; be changed later if desired)

   mjdlist = mjdlist[uniq(mjdlist, sort(mjdlist))]
   mjdlist = strtrim(strcompress(string(mjdlist,format='(99a)')),2)
   sxaddpar, merged_hdr, 'MJDLIST', mjdlist, after='MJD'

   ;----------
   ; Add new header cards

   sxaddpar, merged_hdr, 'VERSCOMB', idlspec2d_version(), $
    ' Version of idlspec2d for combining multiple spectra', after='VERS2D'
   sxaddpar, merged_hdr, 'NEXP', nfiles, $
    ' Number of exposures in this file', before='EXPTIME'
   for ifile=0,nfiles-1 do $
    sxaddpar, merged_hdr, string('EXPID',ifile, format='(a5,i2.2)'), $
      label[ifile], ' ID string for exposure '+strtrim(string(ifile),2), $
      before='EXPTIME'

   sxaddpar, merged_hdr, 'EXPTIME', min(exptimevec), $
    ' Minimum of exposure times for all cameras'
   for icam=0, ncam-1 do $
    sxaddpar, merged_hdr, 'EXPT_'+camnames[icam], exptimevec[icam], $
     ' '+camnames[icam]+' camera exposure time (seconds)', before='EXPTIME'

   return
end
;------------------------------------------------------------------------------
