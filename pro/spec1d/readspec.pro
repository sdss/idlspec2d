;+
; NAME:
;   readspec
;
; PURPOSE:
;   Routine for reading 2D/1D spectro outputs at Princeton
;
; CALLING SEQUENCE:
;   readspec, plate, fiber, [mjd=, flux=, flerr=, invvar=, $
;    andmask=, ormask=, plugmap=, loglam=, wave=, tsobj=, root_dir=, /silent ]
;
; INPUTS:
;   plate      - Plate number(s)
;
; OPTIONAL INPUTS:
;   fiber      - Fiber number(s), 1-indexed; if not set, or zero, then
;                read all fibers for each plate.
;   mjd        - MJD number(s); if not set, then select the most recent
;                data for this plate (largest MJD).
;   silent     - Set to read files in a (more) silent way.
;   root_dir   - Root directory for reduced spectro data files, which
;                are assumed to be ROOT_DIR/pppp/spPlate-pppp-mmmmm.fits,
;                where pppp=plate number and mmmm=MJD.
;                Default to '/data/spectro/2d_3c'.
;   silent     - If set, then call MRDFITS with /SILENT.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   flux       - Flux [NPIXEL,NFIBER]
;   flerr      - Flux error [NPIXEL,NFIBER]
;   invvar     - Inverse variance [NPIXEL,NFIBER]
;   andmask    - AND-mask [NPIXEL,NFIBER]
;   ormask     - OR-mask [NPIXEL,NFIBER]
;   plugmap    - Plug-map entries [NFIBER]
;   loglam     - Log10-wavelength in log10-Angstroms [NPIXEL,NFIBER]
;   wave       - Wavelength in Angstroms [NPIXEL,NFIBER]
;   tsobj      - tsObj-structure output [NFIBER]
;
; COMMENTS:
;   One can input PLATE and FIBER as vectors, in which case there must
;   be a one-to-one correspondence between them.  Or, one can input FIBER
;   numbers as a vector, in which case the same PLATE is used for all.
;   Or, one can input PLATE as a vector, in which case the same FIBER is
;   read for all.
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
;   27-Jun-2000  silent keyword added - Doug Finkbeiner
;-
;------------------------------------------------------------------------------
pro readspec1, plate, range, mjd=mjd, silent=silent, flux=flux, flerr=flerr, $
 invvar=invvar, andmask=andmask, ormask=ormask, plugmap=plugmap, $
 loglam=loglam, wave=wave, tsobj=tsobj, root_dir=root_dir

   common com_readspec, q_flux, q_flerr, q_invvar, q_andmask, q_ormask, $
    q_plugmap, q_loglam, q_wave, q_tsobj

   platestr = string(plate,format='(i4.4)')
   if (NOT keyword_set(mjd)) then mjdstr = '*' $
    else mjdstr = string(mjd,format='(i5.5)')

   filename = 'spPlate-' + platestr + '-' + mjdstr + '.fits'
   filename = findfile(filepath(filename, root_dir=root_dir, $
    subdirectory=platestr), count=ct)

   if (ct GT 1) then begin
      i = reverse(sort(filename))
      filename = filename[i[0]]
   endif else begin
      filename = filename[0]
   endelse

   if (NOT keyword_set(filename)) then begin
      flux = 0
      flerr = 0
      invvar = 0
      andmask = 0
      ormask = 0
      plugmap = 0
      loglam = 0
      wave = 0
      tsobj = 0
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

   if (q_plugmap) then begin
      plugmap = mrdfits(filename, 4, range=range, silent=silent, $
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

   return
end

;------------------------------------------------------------------------------
pro readspec, plate, fiber, mjd=mjd, silent=silent, flux=flux, flerr=flerr, $
 invvar=invvar, andmask=andmask, ormask=ormask, plugmap=plugmap, $
 loglam=loglam, wave=wave, tsobj=tsobj, root_dir=root_dir

   if (n_params() LT 1) then begin
      print, 'Syntax: readspec, plate, [ fiber, mjd=, /silent, flux=, flerr=, invvar=, $'
      print, ' andmask=, ormask=, plugmap=, loglam=, wave=, tsobj= ] '
      return
   endif

   ; This common block specifies which keywords will be returned.
   common com_readspec, q_flux, q_flerr, q_invvar, q_andmask, q_ormask, $
    q_plugmap, q_loglam, q_wave, q_tsobj

   if (NOT keyword_set(root_dir)) then root_dir = '/data/spectro/2d_3c'

   q_flux = arg_present(flux)
   q_flerr = arg_present(flerr)
   q_invvar = arg_present(invvar)
   q_andmask = arg_present(andmask)
   q_ormask = arg_present(ormask)
   q_plugmap = arg_present(plugmap) OR arg_present(tsobj)
   q_loglam = arg_present(loglam)
   q_wave = arg_present(wave)
   q_tsobj = arg_present(tsobj)

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
      plugmap1 = 0
      loglam1 = 0
      wave1 = 0
      tsobj1 = 0

      indx = where(platevec EQ platenums[ifile] AND mjdvec EQ mjdnums[ifile])
      irow = fibervec[indx] - 1
      i1 = min(irow)
      i2 = max(irow)

      if (keyword_set(silent)) then print, '+', format='(A,$)'

      readspec1, platenums[ifile], [i1,i2], mjd=mjdnums[ifile], $
       silent=silent, flux=flux1, flerr=flerr1, invvar=invvar1, $
        andmask=andmask1, ormask=ormask1, plugmap=plugmap1, $
        loglam=loglam1, wave=wave1, tsobj=tsobj1, root_dir=root_dir

      if (ifile EQ 0) then begin
         allindx = indx
         if (q_flux) then flux = flux1[*,irow-i1]
         if (q_flerr) then flerr = flerr1[*,irow-i1]
         if (q_invvar) then invvar = invvar1[*,irow-i1]
         if (q_andmask) then andmask = andmask1[*,irow-i1]
         if (q_ormask) then ormask = ormask1[*,irow-i1]
         if (q_plugmap) then plugmap = plugmap1[irow-i1]
         if (q_loglam) then loglam = loglam1[*,irow-i1]
         if (q_wave) then wave = wave1[*,irow-i1]
         if (q_tsobj) then tsobj = tsobj1[irow-i1]
      endif else begin
         allindx = [allindx, indx]
         if (q_flux) then spec_append, flux, flux1[*,irow-i1]
         if (q_flerr) then spec_append, flerr, flerr1[*,irow-i1]
         if (q_invvar) then spec_append, invvar, invvar1[*,irow-i1]
         if (q_andmask) then spec_append, andmask, andmask1[*,irow-i1]
         if (q_ormask) then spec_append, ormask, ormask1[*,irow-i1]
         if (q_plugmap) then plugmap = struct_append(plugmap, plugmap1[irow-i1])
         if (q_loglam) then spec_append, loglam, loglam1[*,irow-i1]
         if (q_wave) then spec_append, wave, wave1[*,irow-i1]
         if (q_tsobj) then tsobj = struct_append(tsobj, tsobj1[irow-i1])
      endelse
   endfor

   ; Re-sort the data
   if (q_flux) then flux[*,[allindx]] = flux[*]
   if (q_flerr) then flerr[*,[allindx]] = flerr[*]
   if (q_invvar) then invvar[*,[allindx]] = invvar[*]
   if (q_andmask) then andmask[*,[allindx]] = andmask[*]
   if (q_ormask) then ormask[*,[allindx]] = ormask[*]
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

   if (keyword_set(silent)) then print

   return
end
;------------------------------------------------------------------------------
