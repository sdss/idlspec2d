;+
; NAME:
;   readspec
;
; PURPOSE:
;   Routine for reading 2D/1D spectro outputs at Princeton
;
; CALLING SEQUENCE:
;   readspec, plate, fiber, [mjd=, flux=, flerr=, invvar=, $
;    andmask=, ormask=, disp=, plugmap=, loglam=, wave=, tsobj=, $
;    zans=, /silent ]
;
; INPUTS:
;   plate      - Plate number(s)
;
; OPTIONAL INPUTS:
;   fiber      - Fiber number(s), 1-indexed; if not set, or zero, then
;                read all fibers for each plate.  We assume that there
;                are exactly 640 fibers.
;   mjd        - MJD number(s); if not set, then select the most recent
;                data for this plate (largest MJD).
;   silent     - If set, then call MRDFITS with /SILENT.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   mjd        - If not specified, then this is returned as an array of one
;                MJD per object.
;   flux       - Flux [NPIXEL,NFIBER]
;   flerr      - Flux error [NPIXEL,NFIBER]
;   invvar     - Inverse variance [NPIXEL,NFIBER]
;   andmask    - AND-mask [NPIXEL,NFIBER]
;   ormask     - OR-mask [NPIXEL,NFIBER]
;   disp       - Wavelength dispersion [NPIXEL,NFIBER]
;   plugmap    - Plug-map entries [NFIBER]
;   loglam     - Log10-wavelength in log10-Angstroms [NPIXEL,NFIBER]
;   wave       - Wavelength in Angstroms [NPIXEL,NFIBER]
;   tsobj      - tsObj-structure output [NFIBER]
;   zans       - Redshift output structure [NFIBER]
;
; COMMENTS:
;   One can input PLATE and FIBER as vectors, in which case there must
;   be a one-to-one correspondence between them.  Or, one can input FIBER
;   numbers as a vector, in which case the same PLATE is used for all.
;   Or, one can input PLATE as a vector, in which case the same FIBER is
;   read for all.
;
;   The environment variable SPECTRO_DATA must be set to tell this routine
;   where to find the data.  The reduced spectro data files are assumed to
;   be $SPECTRO_DATA/pppp/spPlate-pppp-mmmmm.fits, where pppp=plate number
;   and mmmm=MJD.
;
;   The tsObj files are assumed to be in the directory $SPECTRO_DATA/plates.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   copy_struct_inx
;   headfits()
;   mrdfits
;   plug2tsobj()
;   spec_append
;   struct_append()
;
; INTERNAL SUPPORT ROUTINES:
;   readspec1
;
; REVISION HISTORY:
;   25-Jun-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro readspec1, plate, range, mjd=mjd, flux=flux, flerr=flerr, invvar=invvar, $
 andmask=andmask, ormask=ormask, disp=disp, plugmap=plugmap, $
 loglam=loglam, wave=wave, tsobj=tsobj, zans=zans, $
 root_dir=root_dir, silent=silent

   common com_readspec, q_flux, q_flerr, q_invvar, q_andmask, q_ormask, $
    q_disp, q_plugmap, q_loglam, q_wave, q_tsobj, q_zans, q_mjd

   platestr = string(plate,format='(i4.4)')
   if (NOT keyword_set(mjd)) then mjdstr = '*' $
    else mjdstr = string(mjd,format='(i5.5)')

   filename = 'spPlate-' + platestr + '-' + mjdstr + '.fits'
   filename = findfile(filepath(filename, root_dir=root_dir, $
    subdirectory=platestr), count=ct)

   if (ct GT 1) then filename = filename[ (reverse(sort(filename)))[0] ] $
    else filename = filename[0]

   if (NOT keyword_set(filename)) then begin
      flux = 0
      flerr = 0
      invvar = 0
      andmask = 0
      ormask = 0
      disp = 0
      plugmap = 0
      loglam = 0
      wave = 0
      tsobj = 0
      zans = 0
      return
   end

   if (q_flux) then begin
      flux = mrdfits(filename, 0, hdr, range=range, silent=silent)
   endif

   if (q_invvar OR q_flerr) then begin
      invvar = mrdfits(filename, 1, range=range, silent=silent)
      if (q_flerr) then begin
         i = where(invvar GT 0)
         flerr = 0 * invvar
         if (i[0] NE -1) then flerr[i] = 1 / sqrt(invvar[i])
      endif
   endif

   if (q_andmask) then begin
      andmask = mrdfits(filename, 2, range=range, silent=silent)
   endif

   if (q_ormask) then begin
      ormask = mrdfits(filename, 3, range=range, silent=silent)
   endif

   if (q_disp) then begin
      disp = mrdfits(filename, 4, range=range, silent=silent)
   endif

   if (q_plugmap) then begin
      plugmap = mrdfits(filename, 5, range=range, silent=silent, $
       structyp='PLUGMAPOBJ')
   endif

   if (q_loglam OR q_wave) then begin
      if (NOT keyword_set(hdr)) then hdr = headfits(filename)
      naxis1 = sxpar(hdr, 'NAXIS1')
      naxis2 = sxpar(hdr, 'NAXIS2')
      coeff0 = sxpar(hdr, 'COEFF0')
      coeff1 = sxpar(hdr, 'COEFF1')
      loglam = coeff0 + coeff1 * findgen(naxis1)
      if (keyword_set(range)) then begin
         loglam = rebin(loglam, naxis1, range[1]-range[0]+1)
      endif else begin
         loglam = rebin(loglam, naxis1, naxis2)
      endelse
      if (q_wave) then wave = 10^loglam
   endif

   if (q_tsobj) then begin
      tsobj = plug2tsobj(plate, plugmap=plugmap)
   endif

   if (q_zans) then begin
      zfile = 'spZbest-' + platestr + '-' + mjdstr + '.fits'
      zfile = findfile(filepath(zfile, root_dir=root_dir, $
       subdirectory=platestr), count=ct)

      if (ct GT 1) then zfile = zfile[ (reverse(sort(zfile)))[0] ] $
       else zfile = zfile[0]

      zans = mrdfits(zfile, 1, range=range)
   endif

   if (q_mjd) then begin
      if (NOT keyword_set(hdr)) then hdr = headfits(filename)
      mjd = sxpar(hdr, 'MJD')
   endif

   return
end

;------------------------------------------------------------------------------
pro readspec, plate, fiber, mjd=mjd, flux=flux, flerr=flerr, invvar=invvar, $
 andmask=andmask, ormask=ormask, disp=disp, plugmap=plugmap, $
 loglam=loglam, wave=wave, tsobj=tsobj, zans=zans, silent=silent

   if (n_params() LT 1) then begin
      print, 'Syntax: readspec, plate, [ fiber, mjd=, flux=, flerr=, invvar=, $'
      print, ' andmask=, ormask=, disp=, plugmap=, loglam=, wave=, tsobj=, $'
      print, ' zans=, /silent ] '
      return
   endif

   ; This common block specifies which keywords will be returned.
   common com_readspec, q_flux, q_flerr, q_invvar, q_andmask, q_ormask, $
    q_disp, q_plugmap, q_loglam, q_wave, q_tsobj, q_zans, q_mjd

   root_dir = getenv('SPECTRO_DATA')
   if (NOT keyword_set(root_dir)) then $
    message, 'Environment variable SPECTRO_DATA must be set!'

   q_flux = arg_present(flux)
   q_flerr = arg_present(flerr)
   q_invvar = arg_present(invvar)
   q_andmask = arg_present(andmask)
   q_ormask = arg_present(ormask)
   q_disp = arg_present(disp)
   q_plugmap = arg_present(plugmap) OR arg_present(tsobj)
   q_loglam = arg_present(loglam)
   q_wave = arg_present(wave)
   q_tsobj = arg_present(tsobj)
   q_zans = arg_present(zans)
   q_mjd = arg_present(mjd) AND (keyword_set(mjd) EQ 0)

   if (NOT keyword_set(fiber)) then fiber = lindgen(640) + 1
;      nfiber = n_elements(fiber)
;      if (nplate GT 1 AND nfiber GT 1 AND nplate NE nfiber) then $
;       message, 'Number of elements in PLATE and FIBER must agree or be 1'
;   endif else begin
;      nfiber = 1
;      fiber = 0 ; This will force a read of all fiber numbers
;   endelse

   nvec = n_elements(plate) > n_elements(fiber)
   if (n_elements(plate) GT 1) then platevec = plate $
    else platevec = lonarr(nvec) + plate[0]
   if (n_elements(fiber) GT 1) then fibervec = fiber $
    else fibervec = lonarr(nvec) + fiber[0]
   if (keyword_set(mjd)) then mjdvec = lonarr(nvec) + mjd $
    else mjdvec = lonarr(nvec)

   ; Find unique plate+MJD combinations, since each has its own data file
   sortstring = strtrim(string(platevec),2) + '-' + strtrim(string(mjdvec),2)
   isort = sort(sortstring)
   iuniq = uniq(sortstring[isort])
   platenums = platevec[ isort[iuniq] ]
   mjdnums = mjdvec[ isort[iuniq] ]
   nfile = n_elements(platenums)

   for ifile=0, nfile-1 do begin
      flux1 = 0
      flerr1 = 0
      invvar1 = 0
      andmask1 = 0
      ormask1 = 0
      disp1 = 0
      plugmap1 = 0
      loglam1 = 0
      wave1 = 0
      tsobj1 = 0
      zans1 = 0

      indx = where(platevec EQ platenums[ifile] AND mjdvec EQ mjdnums[ifile])
      irow = fibervec[indx] - 1
      i1 = min(irow)
      i2 = max(irow)

      if (keyword_set(silent)) then print, '+', format='(A,$)'

      mjd1 = mjdnums[ifile]
      readspec1, platenums[ifile], [i1,i2], mjd=mjd1, $
       flux=flux1, flerr=flerr1, invvar=invvar1, andmask=andmask1, $
       ormask=ormask1, disp=disp1, plugmap=plugmap1, loglam=loglam1, $
       wave=wave1, tsobj=tsobj1, zans=zans1, root_dir=root_dir, silent=silent

      if (ifile EQ 0) then begin
         allindx = indx
         if (q_flux) then flux = flux1[*,irow-i1]
         if (q_flerr) then flerr = flerr1[*,irow-i1]
         if (q_invvar) then invvar = invvar1[*,irow-i1]
         if (q_andmask) then andmask = andmask1[*,irow-i1]
         if (q_ormask) then ormask = ormask1[*,irow-i1]
         if (q_disp) then disp = disp1[*,irow-i1]
         if (q_plugmap) then plugmap = plugmap1[irow-i1]
         if (q_loglam) then loglam = loglam1[*,irow-i1]
         if (q_wave) then wave = wave1[*,irow-i1]
         if (q_tsobj) then tsobj = tsobj1[irow-i1]
         if (q_zans) then zans = zans1[irow-i1]
         if (q_mjd) then mjd = mjd1[irow-i1]
      endif else begin
         allindx = [allindx, indx]
         if (q_flux) then spec_append, flux, flux1[*,irow-i1]
         if (q_flerr) then spec_append, flerr, flerr1[*,irow-i1]
         if (q_invvar) then spec_append, invvar, invvar1[*,irow-i1]
         if (q_andmask) then spec_append, andmask, andmask1[*,irow-i1]
         if (q_ormask) then spec_append, ormask, ormask1[*,irow-i1]
         if (q_disp) then spec_append, disp, disp1[*,irow-i1]
         if (q_plugmap) then plugmap = struct_append(plugmap, plugmap1[irow-i1])
         if (q_loglam) then spec_append, loglam, loglam1[*,irow-i1]
         if (q_wave) then spec_append, wave, wave1[*,irow-i1]
         if (q_tsobj) then tsobj = struct_append(tsobj, tsobj1[irow-i1])
         if (q_zans) then zans = struct_append(zans, zans1[irow-i1])
         if (q_mjd) then mjd = [mjd, mjd1[irow-i1]]
      endelse
   endfor

   ; Re-sort the data
   if (q_flux) then flux[*,[allindx]] = flux[*]
   if (q_flerr) then flerr[*,[allindx]] = flerr[*]
   if (q_invvar) then invvar[*,[allindx]] = invvar[*]
   if (q_andmask) then andmask[*,[allindx]] = andmask[*]
   if (q_ormask) then ormask[*,[allindx]] = ormask[*]
   if (q_disp) then disp[*,[allindx]] = disp[*]
   if (q_plugmap) then begin
      if (keyword_set(plugmap[0])) then $
       copy_struct_inx, plugmap, plugmap, index_to=allindx
   endif
   if (q_loglam) then loglam[*,[allindx]] = loglam[*]
   if (q_wave) then wave[*,[allindx]] = wave[*]
   if (q_tsobj) then begin
      if (keyword_set(tsobj[0])) then $
       copy_struct_inx, tsobj, tsobj, index_to=allindx
   endif
   if (q_zans) then begin
      if (keyword_set(zans[0])) then $
       copy_struct_inx, zans, zans, index_to=allindx
   endif
   if (q_mjd) then mjd[allindx] = mjd[*]

   if (keyword_set(silent)) then print

   return
end
;------------------------------------------------------------------------------
