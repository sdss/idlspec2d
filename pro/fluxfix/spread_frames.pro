
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
    flux = flux, ivar = ivar, wave = wave, dispersion = dispersion, $
    pixelmask = pixelmask, plugmap = plugmap, plugtag = plugtag, $
    camerasvec = camerasvec, label = label, filenum = filenum,  $
    expid = expid, sn2 = sn2, exptimevec = exptimevec, mjdlist = mjdlist, $
    hdrarr = hdrarr
 
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
   ; Is there a way to keep track of where the photo mags are from???
   ; In tsObj header?
  
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

      traceset2xy, tempwset, junk, tempwave
      traceset2xy, tempwset, junk-0.5, lowerwave
      traceset2xy, tempwset, junk+0.5, upperwave
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
             tempplugtag[ifib].tsobjid = tsobjid[match]
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
         wave = tempwave
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
         wave = [[wave], [tempwave]]
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

   return
end
;------------------------------------------------------------------------------
