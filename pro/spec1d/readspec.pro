;+
; NAME:
;   readspec
;
; PURPOSE:
;   Routine for reading 2D/1D spectro outputs at Princeton
;
; CALLING SEQUENCE:
;   readspec, plate, fiber, [mjd=, flux=, flerr=, invvar=, $
;    andmask=, ormask=, plugmap=, loglam=, wave=, tsobj= ]
;
; INPUTS:
;   plate      - Plate number(s)
;
; OPTIONAL INPUTS:
;   fiber      - Fiber number(s), 1-indexed; if not set, or zero, then
;                read all fibers for each plate.
;   mjd        - MJD number(s); if not set, then select the most recent
;                data for this plate (largest MJD).
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
;   headfits()
;   mrdfits
;   plug2tsobj()
;   struct_append()
;
; INTERNAL SUPPORT ROUTINES:
;   spec_append
;   readspec1
;
; REVISION HISTORY:
;   25-Jun-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
; Append the array ARG2 to the array ARG1.
; If the first dimension of these arrays is the same, then append
; as [[ARG1],[ARG2]].  If the first dimension of one array is larger
; than the other, then increase the first dimension of the smaller array.

pro spec_append, arg1, arg2

   if (n_elements(arg1) EQ 0) then begin
      arg1 = arg2
      return
   endif
   if (n_elements(arg2) EQ 0) then return

   dims1 = size(arg1, /dimens)
   dims2 = size(arg2, /dimens)
   itype = size(arg1, /type)

   if (dims2[0] EQ dims1[0]) then begin
      ; Case where the new data has the same number of points
      arg1 = [[arg1], [arg2]]
   endif else if (dims2[0] LT dims1[0]) then begin
      ; Case where the new data has fewer points
      dims3 = dims2
      dims3[0] = dims1[0]
      arg3 = make_array(dimension=dims3, type=itype)
      arg3[0:dims2[0]-1,*] = arg2
      arg1 = [[arg1], [arg3]]
   endif else begin
      ; Case where the new data has more points
      dims3 = dims1
      dims3[0] = dims2[0]
      arg3 = make_array(dimension=dims3, type=itype)
      arg3[0:dims1[0]-1,*] = arg1
      arg1 = [[arg3], [arg2]]
   endelse

   return
end

;------------------------------------------------------------------------------
pro readspec1, plate, range, mjd=mjd, flux=flux, flerr=flerr, invvar=invvar, $
 andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=loglam, wave=wave, $
 tsobj=tsobj

   common com_readspec, q_flux, q_flerr, q_invvar, q_andmask, q_ormask, $
    q_plugmap, q_loglam, q_wave, q_tsobj

   platestr = string(plate,format='(i4.4)')
   if (NOT keyword_set(mjd)) then mjdstr = '*' $
    else mjdstr = string(mjd,format='(i5.5)')

   dirname = '/data/spectro/2d_3c/' + platestr + '/2dnew'
   filename = 'spPlate-' + platestr + '-' + mjdstr + '.fits'
   filename = findfile(filepath(filename, root_dir=dirname), count=ct)

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
      flux = mrdfits(filename, 0, hdr, range=range)
   endif

   if (q_invvar OR q_flerr) then begin
      invvar = mrdfits(filename, 1, range=range)
      if (q_flerr) then begin
         i = where(invvar GT 0)
         flerr = 0 * invvar
         if (i[0] NE -1) then flerr[i] = 1 / sqrt(invvar[i])
      endif
   endif

   if (q_andmask) then begin
      andmask = mrdfits(filename, 2, range=range)
   endif

   if (q_ormask) then begin
      ormask = mrdfits(filename, 3, range=range)
   endif

   if (q_plugmap) then begin
      plugmap = mrdfits(filename, 4, range=range, $
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
pro readspec, plate, fiber, mjd=mjd, flux=flux, flerr=flerr, invvar=invvar, $
 andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=loglam, wave=wave, $
 tsobj=tsobj

   if (n_params() LT 1) then begin
      print, 'Syntax: readspec, plate, [ fiber, mjd=, flux=, flerr=, invvar=, $'
      print, ' andmask=, ormask=, plugmap=, loglam=, wave=, tsobj= ] '
      return
   endif

   ; This common block specifies which keywords will be returned.
   common com_readspec, q_flux, q_flerr, q_invvar, q_andmask, q_ormask, $
    q_plugmap, q_loglam, q_wave, q_tsobj

   q_flux = arg_present(flux)
   q_flerr = arg_present(flerr)
   q_invvar = arg_present(invvar)
   q_andmask = arg_present(andmask)
   q_ormask = arg_present(ormask)
   q_plugmap = arg_present(plugmap) OR arg_present(tsobj)
   q_loglam = arg_present(loglam)
   q_wave = arg_present(wave)
   q_tsobj = arg_present(tsobj)

   nplate = n_elements(plate)
   if (keyword_set(fiber)) then begin
      nfiber = n_elements(fiber)
      if (nplate GT 1 AND nfiber GT 1 AND nplate NE nfiber) then $
       message, 'Number of elements in PLATE and FIBER must agree or be 1'
   endif else begin
      nfiber = 1
      fiber = 0 ; This will force a read of all fiber numbers
   endelse

   if (nplate GT 1 OR nfiber GT 1) then begin
      ; Call this routine recursively...
      nvec = nplate > nfiber
      platevec = lonarr(nvec) + plate
      fibervec = lonarr(nvec) + fiber

      if (keyword_set(mjd)) then mjdvec = lonarr(nvec) + mjd $
       else mjdvec = lonarr(nvec)

      for ifiber=0, nvec-1 do begin
         flux1 = 0
         flerr1 = 0
         invvar1 = 0
         andmask1 = 0
         ormask1 = 0
         plugmap1 = 0
         loglam1 = 0
         wave1 = 0
         tsobj1 = 0

         if (fibervec[ifiber] EQ 0) then begin
            range = 0
         endif else begin
            irow = fibervec[ifiber] -1
            range = [irow,irow]
         endelse
         readspec1, platevec[ifiber], range, mjd=mjdvec[ifiber], $
          flux=flux1, flerr=flerr1, invvar=invvar1, $
           andmask=andmask1, ormask=ormask1, plugmap=plugmap1, $
           loglam=loglam1, wave=wave1, tsobj=tsobj1
         if (ifiber EQ 0) then begin
            if (q_flux) then flux = flux1
            if (q_flerr) then flerr = flerr1
            if (q_invvar) then invvar = invvar1
            if (q_andmask) then andmask = andmask1
            if (q_ormask) then ormask = ormask1
            if (q_plugmap) then plugmap = plugmap1
            if (q_loglam) then loglam = loglam1
            if (q_wave) then wave = wave1
            if (q_tsobj) then tsobj = tsobj1
         endif else begin
            if (q_flux) then spec_append, flux, flux1
            if (q_flerr) then spec_append, flerr, flerr1
            if (q_invvar) then spec_append, invvar, invvar1
            if (q_andmask) then spec_append, andmask, andmask1
            if (q_ormask) then spec_append, ormask, ormask1
            if (q_plugmap) then outstruct = struct_append(plugmap, plugmap1)
            if (q_loglam) then spec_append, loglam, loglam1
            if (q_wave) then spec_append, wave, wave1
            if (q_tsobj) then tsobj = struct_append(tsobj, tsobj1)
         endelse
      endfor
   endif else begin
      if (fiber EQ 0) then begin
         range = 0
      endif else begin
         irow = fiber -1
         range = [irow,irow]
      endelse
      readspec1, plate, range, mjd=mjd, $
       flux=flux, flerr=flerr, invvar=invvar, $
        andmask=andmask, ormask=ormask, plugmap=plugmap, $
        loglam=loglam, wave=wave, tsobj=tsobj
   endelse

   return
end
;------------------------------------------------------------------------------
