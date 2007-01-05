;+
; NAME:
;   readonespec
;
; PURPOSE:
;   Routine for reading single exposure spectra from Spectro-2D outputs
;
; CALLING SEQUENCE:
;   readonespec, plate, fiber, [mjd=, cameras=, flux=, flerr=, invvar=, $
;    mask=, disp=, sky=, loglam=, wave=, ximg=, synflux=, lineflux=, $
;    objhdr=, framehdr=, expnum=, topdir=, path=, /silent ]
;
; INPUTS:
;   plate      - Plate number (scalar)
;   fiber      - Fiber number (scalar)
;
; OPTIONAL INPUTS:
;   mjd        - MJD number; if not set, then select the most recent
;                data for this plate (largest MJD).
;   cameras    - If specified, then only match to either the blue ('b')
;                or red ('r')
;   topdir     - Top-level directory for data; default to the environment
;                variable $SPECTRO_DATA.
;   path       - Override all path information with this directory name.
;   silent     - If set, then call MRDFITS with /SILENT.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   mjd        - If not specified, then this is returned
;   flux       - Flux [NPIXEL,NFILE]
;   flerr      - Flux error [NPIXEL,NFILE]
;   invvar     - Inverse variance [NPIXEL,NFILE]
;   mask       - AND-mask [NPIXEL,NFILE]
;   disp       - Wavelength dispersion [NPIXEL,NFILE]
;   sky        - Sky flux [NPIXEL,NFILE]
;   loglam     - Log10-wavelength in log10-Angstroms [NPIXEL,NFILE]
;   wave       - Wavelength in Angstroms [NPIXEL,NFIBER]
;   ximg       - X position on the CCD [NPIXEL,NFIBER]
;   synflux    - Best-fit synthetic eigen-spectrum [NPIXEL,NFILE];
;                return vectors of zeros if the Spectro-1D files
;                cannot be found
;   lineflux   - Best-fit emission lines fits + background terms [NPIXEL,NFILE];
;                return vectors of zeros if the Spectro-1D files
;                cannot be found
;   objhdr     - The FITS header from the spPlate file
;   framehdr   - Pointer array to the FITS headers from all the spCFrame files
;   expnum     - Exposure number for each NFILE spectra
;
; COMMENTS:
;   The environment variable SPECTRO_DATA must be set to tell this routine
;   where to find the data.  The list of spCFrame files to read are determined
;   from the EXPID header keywords in the spPlate file.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $SPECTRO_DATA/$PLATE/spPlate-$PLATE-$MJD.fits
;   $SPECTRO_DATA/$PLATE/spCFrame-$CAMERA-$EXPOSURE.fits*
;
; PROCEDURES CALLED:
;   combine1fiber()
;   readspec
;   spframe_read
;   splog
;
; REVISION HISTORY:
;   10-Feb-2004  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro readonespec, plate, fiber, mjd=mjd, cameras=cameras, $
 flux=flux, flerr=flerr, invvar=invvar, $
 mask=mask, disp=disp, sky=sky, loglam=loglam, wave=wave, ximg=ximg, $
 synflux=synflux, lineflux=lineflux, objhdr=objhdr, framehdr=framehdr, $
 expnum=expnum, topdir=topdir, path=path, silent=silent

   ; Set default return values
   flux = 0
   flerr = 0
   invvar = 0
   mask = 0
   disp = 0
   sky = 0
   loglam = 0
   wave = 0
   ximg = 0
   synflux = 0
   lineflux = 0
   objhdr = 0
   framehdr = 0

   if (n_elements(plate) NE 1 OR n_elements(fiber) NE 1 $
    OR n_elements(mjd) GT 1) then begin
      print, 'PLATE, FIBER, MJD must be scalars'
      return
   endif
   if (keyword_set(cameras)) then begin
      if (cameras[0] NE 'b' AND cameras[0] NE 'r') then begin
         print, 'Invalid value for CAMERAS'
         return
      endif
   endif

   if (NOT keyword_set(topdir) AND NOT keyword_set(path)) then begin
      topdir = getenv('SPECTRO_DATA')
      if (NOT keyword_set(topdir)) then $
       message, 'Environment variable SPECTRO_DATA must be set!'
   endif

   readspec, plate, mjd=mjd, fiber, objhdr=objhdr, topdir=topdir, path=path, $
    silent=silent
   if (NOT keyword_set(objhdr)) then begin
      print, 'spPlate file not found'
      return
   endif
   if (keyword_set(sxpar(objhdr, 'EXPID00'))) then begin
      print, 'Reductions before Spectro-2D v5 cannot access individual spectra'
      return
   endif

   expid = sxpar(objhdr, 'EXPID*')
   filename = 'spCFrame-' + strmid(expid,0,11) + '.fits*'

   if (fiber LE 320) then begin
      spectroid = '1'
      indx = fiber - 1
   endif else begin
      spectroid = '2'
      indx = fiber - 321
   endelse

   ; Trim to the files that correspond to the spectrograph with this fiber
   itrim = where(strmid(expid,1,1) EQ spectroid, nfile)
   if (nfile GT 0 AND keyword_set(cameras)) then $
    itrim = itrim[where(strmid(expid[itrim],0,1) EQ cameras[0], nfile)>0]
   if (nfile EQ 0) then begin
      print, 'No files on this for spectrograph #', spectroid
      return
   endif
   expid = expid[itrim]
   filename = filename[itrim]

   ; Get the fully-qualified path names for the files, and make
   ; sure the files exist
   platestr = string(plate,format='(i4.4)')
   for ifile=0L, nfile-1 do begin
      if (keyword_set(path)) then $
       tmpname = lookforgzip(filepath(filename[ifile], $
        root_dir=path)) $
      else $
       tmpname = lookforgzip(filepath(filename[ifile], $
        root_dir=topdir, subdirectory=platestr))
      if (NOT keyword_set(tmpname)) then begin
         print, 'File not found: ' + filename[ifile]
         return
      endif
      filename[ifile] = tmpname
   endfor

   if (arg_present(framehdr)) then begin
      framehdr = ptrarr(nfile)
      for ifile=0L, nfile-1 do begin
         spframe_read, filename[ifile], hdr=framehdr1
         framehdr[ifile] = ptr_new(framehdr1)
      endfor
   endif

   if (arg_present(flux)) then begin
      flux = fltarr(2048,nfile)
      for ifile=0L, nfile-1 do begin
         spframe_read, filename[ifile], indx, objflux=flux1
         flux[*,ifile] = flux1
      endfor
   endif
   itrim = where(filename NE '', nfile)
   if (nfile EQ 0) then begin
      print, 'No spCFrame files exist on disk'
      return
   endif
   filename = filename[itrim]

   if (arg_present(invvar) OR arg_present(flerr)) then begin
      invvar = fltarr(2048,nfile)
      for ifile=0L, nfile-1 do begin
         spframe_read, filename[ifile], indx, objivar=invvar1
         invvar[*,ifile] = invvar1
      endfor
      if (arg_present(flerr)) then begin
         qgood = invvar GT 0
         flerr = qgood / sqrt(invvar * qgood + (1-qgood))
      endif
   endif

   if (arg_present(mask)) then begin
      mask = lonarr(2048,nfile)
      for ifile=0L, nfile-1 do begin
         spframe_read, filename[ifile], indx, mask=mask1
         mask[*,ifile] = mask1
      endfor
   endif

   if (arg_present(loglam)) then begin
      loglam = fltarr(2048,nfile)
      for ifile=0L, nfile-1 do begin
         spframe_read, filename[ifile], indx, loglam=loglam1
         loglam[*,ifile] = loglam1
      endfor
   endif

   if (arg_present(loglam) OR arg_present(wave) OR arg_present(synflux) $
    OR arg_present(lineflux)) then begin
      loglam = fltarr(2048,nfile)
      for ifile=0L, nfile-1 do begin
         spframe_read, filename[ifile], indx, loglam=loglam1
         loglam[*,ifile] = loglam1
      endfor
      if (arg_present(wave)) then wave = 10.d^loglam
   endif

   if (arg_present(ximg)) then begin
      ximg = fltarr(2048,nfile)
      for ifile=0L, nfile-1 do begin
         spframe_read, filename[ifile], indx, ximg=ximg1
         ximg[*,ifile] = ximg1
      endfor
   endif

   if (arg_present(disp)) then begin
      disp = fltarr(2048,nfile)
      for ifile=0L, nfile-1 do begin
         spframe_read, filename[ifile], indx, dispimg=disp1
         disp[*,ifile] = disp1
      endfor
   endif

   if (arg_present(sky)) then begin
      sky = fltarr(2048,nfile)
      for ifile=0L, nfile-1 do begin
         spframe_read, filename[ifile], indx, sky=sky1
         sky[*,ifile] = sky1
      endfor
   endif

   if (arg_present(synflux)) then begin
      readspec, plate, fiber, mjd=mjd, loglam=loglam1, synflux=synflux1, $
       topdir=topdir, path=path, silent=silent
      if (keyword_set(synflux1)) then begin
         combine1fiber, loglam1, synflux1, newloglam=loglam, newflux=synflux
         synflux = reform(synflux, size(loglam,/dimens))
      endif else begin
         splog, 'WARNING: Setting missing SYNFLUX vectors to zero'
         synflux = 0 * loglam
      endelse
   endif

   if (arg_present(lineflux)) then begin
      readspec, plate, fiber, mjd=mjd, loglam=loglam1, lineflux=lineflux1, $
       topdir=topdir, path=path, silent=silent
      if (keyword_set(lineflux1)) then begin
         combine1fiber, loglam1, lineflux1, newloglam=loglam, newflux=lineflux
         lineflux = reform(lineflux, size(loglam,/dimens))
      endif else begin
         splog, 'WARNING: Setting missing SYNFLUX vectors to zero'
         lineflux = 0 * loglam
      endelse
   endif

   if (arg_present(expnum)) then expnum = long(strmid(expid,3,8))

   return
end
;------------------------------------------------------------------------------
