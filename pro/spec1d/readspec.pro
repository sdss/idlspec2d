;+
; NAME:
;   readspec
;
; PURPOSE:
;   Routine for reading 2D/1D spectro outputs at Princeton
;
; CALLING SEQUENCE:
;   readspec, plate, fiber, [mjd=, flux=, flerr=, invvar=, $
;    andmask=, ormask=, plugmap=, loglam=, wave= ]
;
; INPUTS:
;   plate      - Plate number(s)
;   fiber      - Fiber number(s)
;
; OPTIONAL INPUTS:
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
pro readspec1, plate, fiber, mjd=mjd, flux=flux, flerr=flerr, invvar=invvar, $
 andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=loglam, wave=wave

   common com_readspec, q_flux, q_flerr, q_invvar, q_andmask, q_ormask, $
    q_plugmap, q_loglam, q_wave

   irow = fiber - 1

   platestr = string(plate,format='(i4.4)')
   if (NOT keyword_set(mjd)) then mjdstr = '*' $
    else mjdstr = string(mjd,format='(i5.5)')

   dirname = '/data/spectro/2d_3c/' + platestr + '/2dnew'
   filename = 'spPlate-' + platestr + mjdstr + '.fits'
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
      return
   end

   if (q_flux) then begin
      flux = mrdfits(filename, 0, hdr, range=[irow,irow])
   endif

   if (q_invvar OR q_flerr) then begin
      invvar = mrdfits(filename, 1, range=[irow,irow])
      if (q_flerr) then begin
         i = where(invvar GT 0)
         flerr = 0 * invvar
         if (i[0] NE -1) then flerr[i] = 1 / sqrt(invvar[i])
      endif
   endif

   if (q_andmask) then begin
      andmask = mrdfits(filename, 2, range=[irow,irow])
   endif

   if (q_ormask) then begin
      ormask = mrdfits(filename, 3, range=[irow,irow])
   endif

   if (q_plugmap) then begin
      plugmap = mrdfits(filename, 4, range=[irow,irow], $
       structyp='PLUGMAPOBJ')
   endif

   if (q_loglam OR q_wave) then begin
      if (NOT keyword_set(hdr)) then hdr = headfits(filename)
      naxis1 = sxpar(hdr, 'NAXIS1')
      coeff0 = sxpar(hdr, 'COEFF0')
      coeff1 = sxpar(hdr, 'COEFF1')
      loglam = coeff0 + coeff1 * findgen(naxis1)
      if (q_wave) then wave = 10^loglam
   endif

   return
end

;------------------------------------------------------------------------------
pro readspec, plate, fiber, mjd=mjd, flux=flux, flerr=flerr, invvar=invvar, $
 andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=loglam, wave=wave

   if (n_params() LT 2) then begin
      print, 'Syntax: readspec, plate, fiber, [mjd=, flux=, flerr=, invvar=, $'
      print, ' andmask=, ormask=, plugmap=, loglam=, wave= ] '
      return
   endif

   ; This common block specifies which keywords will be returned.
   common com_readspec, q_flux, q_flerr, q_invvar, q_andmask, q_ormask, $
    q_plugmap, q_loglam, q_wave

   q_flux = arg_present(flux)
   q_flerr = arg_present(flerr)
   q_invvar = arg_present(invvar)
   q_andmask = arg_present(andmask)
   q_ormask = arg_present(ormask)
   q_plugmap = arg_present(plugmap)
   q_loglam = arg_present(loglam)
   q_wave = arg_present(wave)

   nplate = n_elements(plate)
   nfiber = n_elements(fiber)
   if (nplate GT 1 AND nfiber GT 1 AND nplate NE nfiber) then $
    message, 'Number of elements in PLATE and FIBER must agree or be 1'

   if (nplate GT 1 OR nfiber GT 1) then begin
      ; Call this routine recursively...
      nvec = nplate > nfiber
      platevec = lonarr(nvec) + plate
      fibervec = lonarr(nvec) + fiber
      if (keyword_set(mjd)) then mjdvec = lonarr(nvec) + mjd $
       else mjdvec = lonarr(nvec)

      for ifiber=0, nfiber-1 do begin
         flux1 = 0
         flerr1 = 0
         invvar1 = 0
         andmask1 = 0
         ormask1 = 0
         plugmap1 = 0
         loglam1 = 0
         wave1 = 0
         readspec1, platevec[ifiber], fibervec[ifiber], mjd=mjdvec[ifiber], $
          flux=flux1, flerr=flerr1, invvar=invvar1, $
           andmask=andmask1, ormask=ormask1, plugmap=plugmap1, $
           loglam=loglam1, wave=wave1
         if (ifiber EQ 0) then begin
            if (q_flux) then flux = flux1
            if (q_flerr) then flerr = flerr1
            if (q_invvar) then invvar = invvar1
            if (q_andmask) then andmask = andmask1
            if (q_ormask) then ormask = ormask1
            if (q_plugmap) then plugmap = plugmap1
            if (q_loglam) then loglam = loglam1
            if (q_wave) then wave = wave1
         endif else begin
            if (q_flux) then spec_append, flux, flux1
            if (q_flerr) then spec_append, flerr, flerr1
            if (q_invvar) then spec_append, invvar, invvar1
            if (q_andmask) then spec_append, andmask, andmask1
            if (q_ormask) then spec_append, ormask, ormask1
            if (q_plugmap) then begin
               if (keyword_set(plugmap)) then begin
                  tmp_plug= plugmap[0]
                  if (size(plugmap1,/tname) EQ 'STRUCT') then $
                   struct_assign, plugmap1, tmp_plug $
                  else $
                   struct_assign, {junk:0}, tmp_plug
                  plugmap = [plugmap, tmp_plug]
               endif
            endif
            if (q_loglam) then spec_append, loglam, loglam1
            if (q_wave) then spec_append, wave, wave1
         endelse
      endfor
   endif else begin
      readspec1, plate, fiber, mjd=mjd, $
       flux=flux, flerr=flerr, invvar=invvar, $
        andmask=andmask, ormask=ormask, plugmap=plugmap, $
        loglam=loglam, wave=wave
   endelse

   return
end
;------------------------------------------------------------------------------
