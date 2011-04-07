;+
; NAME:
;   write_uros
;
; PURPOSE:
;   Generate ASCII files of individual exposures on spectra for Uros Seljak
;
; CALLING SEQUENCE:
;   write_uros, plate, [ fiber, ] mjd=
;
; INPUTS:
;   plate      - Plate number(s)
;
; REQUIRED KEYWORDS:
;   mjd        - MJD number(s); if not set, then select the most recent
;                data for this plate (largest MJD).
;
; OPTIONAL INPUTS:
;   fiber      - Fiber number(s), 1-indexed; if not set, or zero, then
;                read all fibers for each plate.  We assume that there
;                are exactly 640 fibers.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The following file is written for each individual (CCD) spectrum
;   of each object:
;     spUros-$PLATE-$MJD-$FIBER-$EXPNUM-$CAMERA.fits
;
;   In general, there will be (3 exposures) x (2 cameras) = 6 spectra
;   for each plate+MJD+fiber combination.
;
;   Output the following for each individual spectrum of each object:
;     log-Wave   : Log10(wavelength [Ang])
;     Flux       : Flux in units of flat-fielded electrons
;     InverseVar : Inverse variance for the above
;     Sky        ; Sky
;     Pixelmask  : Pixel mask
;     Dispers    : Wavelength dispersion in the **native** pixel scale
;                  (e.g., relative to the local wavelength spacing)
;     Calibfac   : Flux-calibration vector for this plate+spectrograph
;     Invercorr  : Flux-correction vector for this exposure relative
;                  to the flux-calibrated (smear) exposure
;     Flatvec    : Flat-field vector, relative to a quartz lamp, combined
;                  with the atmospheric telluric correction (in the red)
;     Rnois      : Read noise in electrons (or photons)
;
;   The fluxes have already been divided by the fiber-to-fiber flats
;   (superfit) and the telluric-correction vector (telluricfactor).
;   Convert the fluxes back to photon number as follows:
;     Flux_photons = Flux * Flatvec
;     Sky_photons = Sky * Flatvec
;     InverseVar_photons = InverseVar / Flatvec^2
;   The read noise portion (Rnois) is already in units of electrons,
;   or equivalently photons.
;
;   Before combining multiple spectra, we re-normalize as follows:
;     Flux_ergs = Flux / (Calibfac * Invercorr)
;     InverseVar = InverseVar * (Calibfac * Invercorr)^2 * Window
;   where Window is a linear apodization of the first 100 pixels
;   of each spectrum (the blue/red overlap region, so the red side
;   of the blue camera, and the blue side of the red camera).
;   Once the flux has been normalized, it is in units of 10^(-17) erg/cm/s/Ang
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $BOSS_SPECTRO_REDUX/$PLATE/spFrame-$CAMERA-$EXPNUM.fits*
;   $BOSS_SPECTRO_REDUX/$PLATE/spFluxcalib-$PLATE-$MJD-$CAMERA.fits
;   $BOSS_SPECTRO_REDUX/$PLATE/spFluxcorr-$EXPNUM-$SPECID.fits
;   $BOSS_SPECTRO_REDUX/$PLATE/spPlate-$PLATE-$MJD.fits
;
; PROCEDURES CALLED:
;   bspline_valu()
;   mrdfits()
;   readspec
;   sxpar()
;   traceset2xy
;   yanny_free
;   yanny_read()
;
; INTERNAL SUPPORT ROUTINES:
;   write_uros1
;
; REVISION HISTORY:
;   14-Jun-2002  Written by David Schlegel, Princeton.
;                This is based upon the code in SPCOADD_FRAMES.
;-
;------------------------------------------------------------------------------
pro write_uros1, plate, fiber, mjd=mjd

   topdir = getenv('BOSS_SPECTRO_REDUX')
   if (NOT keyword_set(topdir)) then $
    message, 'Environment variable BOSS_SPECTRO_REDUX must be set!'

   if (fiber LE 320) then specid = 1 $
    else specid = 2
   platestr = string(plate, format='(i4.4)')

   ;----------
   ; Find the names of individual exposures

   readspec, plate, fiber, mjd=mjd, objhdr=objhdr
   if (NOT keyword_set(objhdr)) then begin
      print, 'WARNING: Object not found ', plate, mjd, fiber
      return
   endif

   icard = where(strmatch(objhdr, 'EXPID*'))
   expspec = long(strmid(objhdr[icard], 12, 1))
   expcams = strmid(objhdr[icard], 11, 2)
   expnums = long(strmid(objhdr[icard], 20, 8))
   iuse = where(expspec EQ specid, nuse)

   for ifile=0, nuse-1 do begin

      framefile = string(expcams[iuse[ifile]], expnums[iuse[ifile]], $
       format='("spFrame-",a2,"-",i8.8,".fits*")')
      corrfile = string(expnums[iuse[ifile]], expspec[iuse[ifile]], $
       format='("spFluxcorr-",i8.8,"-",i1,".fits*")')
      fcalibfile = string(plate,mjd,expcams[iuse[ifile]], $
       format='("spFluxcalib-",i4.4,"-",i5.5,"-",a2,".fits*")')
      outfile = string(plate,mjd,fiber,expnums[iuse[ifile]], $
       expcams[iuse[ifile]], $
       format='("spUros-",i4.4,"-",i5.5,"-",i4.4,"-",i8.8,"-",a2,".dat")')

      framefile = (findfile(filepath(framefile, root_dir=topdir, $
       subdir=platestr)))[0]
      corrfile = (findfile(filepath(corrfile, root_dir=topdir, $
       subdir=platestr)))[0]
      fcalibfile = (findfile(filepath(fcalibfile, root_dir=topdir, $
       subdir=platestr)))[0]

      splog, 'Reading file #', ifile, ': ', framefile
      rownum = (fiber-1) MOD 320
      range = [rownum,rownum]
      tempflux = mrdfits(framefile, 0, temphdr, range=range)
      tempivar = mrdfits(framefile, 1, range=range)
      temppixmask = mrdfits(framefile, 2, range=range)
      wsetall = mrdfits(framefile, 3)
      dispsetall = mrdfits(framefile, 4, range=range)
      tempsky = mrdfits(framefile, 6, range=range)
      tempflat = mrdfits(framefile, 7, range=range)
      temptelluric = mrdfits(framefile, 8, range=range)
      if (keyword_set(temptelluric)) then tempflat = tempflat * temptelluric

      ;----------
      ; Trim the wavelength + dispersion structures to only the requested object

      tempwset = create_struct('func', wsetall.func, $
       'xmin', wsetall.xmin, $
       'xmax', wsetall.xmax, $
       'coeff', wsetall.coeff[*,rownum])
      tempdispset = create_struct('func', dispsetall.func, $
       'xmin', dispsetall.xmin, $
       'xmax', dispsetall.xmax, $
       'coeff', dispsetall.coeff[*,rownum])

      ;----------
      ; Add an additional error term equal to ADDERR of the flux.

      if (keyword_set(adderr)) then begin
         gmask = tempivar NE 0 ; =1 for good points
         tempivar = 1.0 / ( 1.0/(tempivar + (1-gmask)) $
          + (adderr * (tempflux>0))^2 ) * gmask
      endif

      ;----------
      ; Solve for wavelength and lambda-dispersion at each pixel
      ; for this object.

      traceset2xy, tempwset, junk, tempwave

      traceset2xy, tempdispset, junk, tempdispersion

      ;----------
      ; Compute spectro-photometric calibration

      calibhdr = headfits(fcalibfile)
      cwavemin = sxpar(calibhdr, 'WAVEMIN')
      cwavemax = sxpar(calibhdr, 'WAVEMAX')
      calibset = mrdfits(fcalibfile, 1)
      calibfac = bspline_valu(tempwave, calibset)

      ; Set to bad any pixels whose wavelength is outside the known
      ; flux-calibration region.
      ibad = where(tempwave LT alog10(cwavemin) OR tempwave GT alog10(cwavemax))
      if (ibad[0] NE -1) then calibfac[ibad] = 0

      ; Set lower-threshhold on calib-factor
      ; This isn't exactly what is done in SPCOADD_FRAMES,
      ; since the median is taken over all fibers.
      minval = 0.05 * mean(calibfac)
      calibfac = calibfac * (calibfac GT minval)

      ;----------
      ; Compute flux-correction factors

      corrall = mrdfits(corrfile, 1)
      corrset = create_struct('func', corrall.func, $
       'xmin', corrall.xmin, $
       'xmax', corrall.xmax, $
       'coeff', corrall.coeff[*,rownum])
      traceset2xy, corrset, tempwave, corrimg
      invertcorr = 1.0 / corrimg

      ; Set lower-threshhold on calib-factor
      ; This isn't exactly what is done in SPCOADD_FRAMES,
      ; since the median is taken over all fibers.
      medcorr = median(corrimg)
      medcorr = medcorr + (medcorr LE 0) ; prevent divide-by-zeros
      minval = 0.05 / medcorr
      invertcorr = invertcorr * (invertcorr GT minval)

      ;----------
      ; Read the opECalib file to get the gain + read noise

      config_dir = filepath('', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
      ecalibfile = findopfile('opECalib*par', sxpar(temphdr,'MJD'), $
       config_dir, /abort_notfound)
      yanny_read, filepath(ecalibfile, root_dir=config_dir), pdata
      ecalib = *pdata[0]
      yanny_free, pdata
      jj = where(ecalib.camrow EQ sxpar(temphdr,'CAMROW') $
       AND ecalib.camcol EQ sxpar(temphdr,'CAMCOL'))
      ecalib = ecalib[jj[0]]
      if (fiber LE 320) then begin
         rdnoise = ecalib.readnoisedn2 * ecalib.gain2
      endif else begin
         rdnoise = ecalib.readnoisedn3 * ecalib.gain3
      endelse

      ;----------
      ; Write the output file

      splog, 'Writing file ', outfile
      openw, lun, outfile, /get_lun
      printf, lun, '#log-Wave Flux          InverseVar    Sky          ' $
       + ' Pixelmask   Dispers   Calibfac  Invercorr Flatvec   Rnois'
      printf, lun, '#-------- ------------- ------------- -------------' $
       + ' ----------- --------- --------- --------- --------- -----'
      for ipix=0, n_elements(tempflux)-1 do $
       printf, lun, tempwave[ipix], tempflux[ipix], tempivar[ipix], $
        tempsky[ipix], temppixmask[ipix], tempdispersion[ipix], $
        calibfac[ipix], invertcorr[ipix], tempflat[ipix], rdnoise, $
        format='(1x,f8.6,1x,e13.6,1x,e13.6,1x,e13.6,1x,i11,1x,f9.4,1x,f9.4,1x,f9.4,f9.4,f6.2)'
      close, lun
      free_lun, lun
   endfor

   return
end
;------------------------------------------------------------------------------
pro write_uros, plate, fiber, mjd=mjd

   if (n_params() LT 1 OR NOT keyword_set(mjd)) then begin
      doc_library, 'write_uros'
      return
   endif
   if (n_elements(plate) NE n_elements(mjd)) then begin
      print, 'Number of elements in PLATE and MJD must agree'
      return
   endif

   if (NOT keyword_set(fiber)) then begin
      thisplate = (transpose(plate # replicate(1,640)))[*]
      thismjd = (transpose(mjd # replicate(1,640)))[*]
      thisfiber = ((lindgen(640)+1) # replicate(1,n_elements(plate)))[*]
   endif else begin
      if (n_elements(plate) NE n_elements(fiber)) then begin
         print, 'Number of elements in PLATE and FIBER must agree'
         return
      endif
      thisplate = plate
      thismjd = mjd
      thisfiber = fiber
   endelse

   for iobj=0, n_elements(plate)-1 do begin
      write_uros1, thisplate[iobj], thisfiber[iobj], mjd=thismjd[iobj]
   endfor

   return
end
;------------------------------------------------------------------------------
