;+
; NAME:
;   readspec
;
; PURPOSE:
;   Routine for reading 2D/1D spectro outputs at Princeton
;
; CALLING SEQUENCE:
;   readspec, plate, fiber, [mjd=, znum=, flux=, flerr=, invvar=, $
;    andmask=, ormask=, disp=, plugmap=, loglam=, wave=, tsobj=, $
;    zans=, synflux=synflux, /silent ]
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
;   znum       - If set, then return not the best-fit redshift, but the
;                ZNUM-th best-fit; e.g., set ZNUM=2 for second-best fit.
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
;   synflux    - Best-fit synthetic eigen-spectrum [NPIXEL,NFIBER]
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
; DATA FILES:
;   $SPECTRO_DATA/$PLATE/spPlate-$PLATE-$MJD.fits
;   $SPECTRO_DATA/$PLATE/spZbest-$PLATE-$MJD.fits
;   $IDLSPEC2D_DIR/templates/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   copy_struct_inx
;   headfits()
;   lookforgzip()
;   mrdfits
;   plug2tsobj()
;   spec_append
;   struct_append()
;   synthspec()
;
; INTERNAL SUPPORT ROUTINES:
;   rspec_mrdfits()
;   readspec1
;
; REVISION HISTORY:
;   25-Jun-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function rspec_mrdfits, fcb, exten_no, rownums=rownums, _EXTRA=EXTRA

   if (exten_no GT fcb.nextend) then return, 0

   nrows = n_elements(rownums)
   if (nrows EQ 1) then begin
      nchunks = 1
      row_start = rownums[0]
      row_end = rownums[0]
   endif else begin
      rowdiff = rownums[1:nrows-1] - rownums[0:nrows-2]
      i0 = where(rowdiff NE 1, nchunks)
      nchunks = nchunks + 1
      if (nchunks EQ 1) then begin
         row_start = rownums[0]
         row_end = rownums[nrows-1]
      endif else begin
         row_start = rownums[ [0, i0 + 1] ]
         row_end = rownums[ [i0, nrows-1] ]
      endelse
   endelse

   naxis1 = fcb.axis[0,exten_no]
;   naxis2 = fcb.axis[1,exten_no]
   iadd = 0
   for ichunk=0L, nchunks-1 do begin
      nadd = row_end[ichunk] - row_start[ichunk] + 1
      if (exten_no EQ 0 OR fcb.xtension[exten_no] EQ 'IMAGE') then begin
         fits_read, fcb, data1, exten_no=exten_no, $
          first=naxis1*row_start[ichunk], $
          last=naxis1*(row_end[ichunk]+1)-1, _EXTRA=EXTRA
         if (ichunk EQ 0) then $
          alldata = make_array(naxis1, nrows, size=size(data1))
         alldata[*,iadd:iadd+nadd-1] = data1
      endif else if (fcb.xtension[exten_no] EQ 'BINTABLE') then begin
         data1 = mrdfits(fcb.filename, exten_no, $
          range=[row_start[ichunk], row_end[ichunk]], _EXTRA=EXTRA)
         if (ichunk EQ 0) then $
          alldata = replicate(data1[0], nrows)
         if (keyword_set(tag_names(data1,/structure_name))) then $
          alldata[iadd:iadd+nadd-1] = data1 $ ; named structure
         else $
          copy_struct_inx, data1, alldata, index_to=iadd+lindgen(nadd)
      endif
      iadd = iadd + nadd
   endfor

   return, alldata
end

;------------------------------------------------------------------------------
pro readspec1, plate, rownums, mjd=mjd, flux=flux, flerr=flerr, invvar=invvar, $
 andmask=andmask, ormask=ormask, disp=disp, plugmap=plugmap, $
 loglam=loglam, wave=wave, tsobj=tsobj, zans=zans, synflux=synflux, $
 znum=znum, root_dir=root_dir, silent=silent

   common com_readspec, q_flux, q_flerr, q_invvar, q_andmask, q_ormask, $
    q_disp, q_plugmap, q_loglam, q_wave, q_tsobj, q_zans, q_synflux, q_mjd

   platestr = string(plate,format='(i4.4)')
   if (NOT keyword_set(mjd)) then mjdstr = '*' $
    else mjdstr = string(mjd,format='(i5.5)')

   filename = 'spPlate-' + platestr + '-' + mjdstr + '.fits'
   filename = lookforgzip(filepath(filename, root_dir=root_dir, $
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
      synflux = 0
      return
   end

   fits_open, filename, fcb

   if (q_flux) then begin
      flux = rspec_mrdfits(fcb, 0, rownums=rownums, silent=silent)
   endif

   if (q_invvar OR q_flerr) then begin
      invvar = rspec_mrdfits(fcb, 1, rownums=rownums, silent=silent)
      if (q_flerr) then begin
         i = where(invvar GT 0)
         flerr = 0 * invvar
         if (i[0] NE -1) then flerr[i] = 1 / sqrt(invvar[i])
      endif
   endif

   if (q_andmask) then begin
      andmask = rspec_mrdfits(fcb, 2, rownums=rownums, silent=silent)
   endif

   if (q_ormask) then begin
      ormask = rspec_mrdfits(fcb, 3, rownums=rownums, silent=silent)
   endif

   if (q_disp) then begin
      disp = rspec_mrdfits(fcb, 4, rownums=rownums, silent=silent)
   endif

   if (q_plugmap) then begin
      plugmap = rspec_mrdfits(fcb, 5, rownums=rownums, silent=silent, $
       structyp='PLUGMAPOBJ')
   endif

   if (q_loglam OR q_wave) then begin
      if (NOT keyword_set(hdr)) then hdr = headfits(filename)
      naxis1 = sxpar(hdr, 'NAXIS1')
;      naxis2 = sxpar(hdr, 'NAXIS2')
      coeff0 = sxpar(hdr, 'COEFF0')
      coeff1 = sxpar(hdr, 'COEFF1')
      loglam = coeff0 + coeff1 * findgen(naxis1)
      loglam = rebin(loglam, naxis1, n_elements(rownums))
      if (q_wave) then wave = 10^loglam
   endif

   if (q_tsobj) then begin
      tsobj = plug2tsobj(plate, plugmap=plugmap)
   endif

   if (q_zans OR q_synflux) then begin
      if (NOT keyword_set(znum)) then $
       zfile = 'spZbest-' + platestr + '-' + mjdstr + '.fits' $
      else $
       zfile = 'spZall-' + platestr + '-' + mjdstr + '.fits'

      zfile = lookforgzip(filepath(zfile, root_dir=root_dir, $
       subdirectory=platestr), count=ct)
      if (ct GT 1) then zfile = zfile[ (reverse(sort(zfile)))[0] ] $
       else zfile = zfile[0]

      if (keyword_set(zfile)) then begin
         if (NOT keyword_set(znum)) then begin
            fits_open, zfile, zfcb
            zans = rspec_mrdfits(zfcb, 1, rownums=rownums, silent=silent)
            fits_close, zfcb
         endif else begin
            zhdr = headfits(zfile, exten=1)
            nper = sxpar(zhdr,'NAXIS2') / 640L
            fits_open, zfile, zfcb
            zans = rspec_mrdfits(zfcb, 1, rownums=rownums*nper+znum-1, $
             silent=silent)
            fits_close, zfcb
         endelse
      endif
   endif

   if (q_synflux) then begin
      ; Read the synthetic spectrum from the Zbest file if ZNUM is not set.
      if (NOT keyword_set(znum)) then begin
         fits_open, zfile, zfcb
         synflux = rspec_mrdfits(zfcb, 2, rownums=rownums, silent=silent)
         fits_close, zfcb
      endif else begin
         if (NOT keyword_set(hdr)) then hdr = headfits(filename)
         if (keyword_set(zans)) then $
          synflux = synthspec(zans, hdr=hdr)
      endelse
   endif

   if (q_mjd) then begin
      if (NOT keyword_set(hdr)) then hdr = headfits(filename)
      mjd = sxpar(hdr, 'MJD')
   endif

   fits_close, fcb

   return
end

;------------------------------------------------------------------------------
pro readspec, plate, fiber, mjd=mjd, flux=flux, flerr=flerr, invvar=invvar, $
 andmask=andmask, ormask=ormask, disp=disp, plugmap=plugmap, $
 loglam=loglam, wave=wave, tsobj=tsobj, zans=zans, synflux=synflux, $
 znum=znum, silent=silent

   if (n_params() LT 1) then begin
      print, 'Syntax: readspec, plate, [ fiber, mjd=, flux=, flerr=, invvar=, $'
      print, ' andmask=, ormask=, disp=, plugmap=, loglam=, wave=, tsobj=, $'
      print, ' zans=, synflux=, /silent ] '
      return
   endif

   ; This common block specifies which keywords will be returned.
   common com_readspec, q_flux, q_flerr, q_invvar, q_andmask, q_ormask, $
    q_disp, q_plugmap, q_loglam, q_wave, q_tsobj, q_zans, q_synflux, q_mjd

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
   q_synflux = arg_present(synflux)
   q_mjd = arg_present(mjd) AND (keyword_set(mjd) EQ 0)

   nplate = n_elements(plate)
   if (nplate EQ 0) then $
    message, 'PLATE must be defined'
   if (n_elements(mjd) GT 0 AND n_elements(mjd) NE nplate) then $
    message, 'Number of elements in PLATE and MJD must agree'

   if (NOT keyword_set(fiber)) then begin
      ; Special case to read all 640 fibers of each plate
      platevec = reform( [plate] ## replicate(1,640L), 640L*nplate)
      fibervec = reform( (lindgen(640L) + 1) # replicate(1,nplate), 640L*nplate)
      if (keyword_set(mjd)) then $
       mjdvec = reform( [mjd] ## replicate(1,640L), 640L*nplate) $
      else $
       mjdvec = lonarr(640L*nplate)
   endif else begin
      nfiber = n_elements(fiber)
      if (nplate GT 1 AND nfiber GT 1 AND nplate NE nfiber) then $
       message, 'Number of elements in PLATE and FIBER must agree or be 1'

      nvec = nplate > nfiber
      if (nplate GT 1) then platevec = plate $
       else platevec = lonarr(nvec) + plate[0]
      if (nfiber GT 1) then fibervec = fiber $
       else fibervec = lonarr(nvec) + fiber[0]
      if (keyword_set(mjd)) then mjdvec = lonarr(nvec) + mjd $
       else mjdvec = lonarr(nvec)
   endelse

   ; Find unique plate+MJD combinations, since each has its own data file
   sortstring = strtrim(string(platevec),2) + '-' + strtrim(string(mjdvec),2)
   isort = sort(sortstring)
   iuniq = uniq(sortstring[isort])
   platenums = platevec[ isort[iuniq] ]
   mjdnums = mjdvec[ isort[iuniq] ]
   nfile = n_elements(platenums)

   for ifile=0L, nfile-1 do begin
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
      synflux1 = 0

      indx = where(platevec EQ platenums[ifile] AND mjdvec EQ mjdnums[ifile])
      irow = fibervec[indx] - 1

      if (keyword_set(silent)) then print, '+', format='(A,$)'

      mjd1 = mjdnums[ifile]
      readspec1, platenums[ifile], irow, mjd=mjd1, $
       flux=flux1, flerr=flerr1, invvar=invvar1, andmask=andmask1, $
       ormask=ormask1, disp=disp1, plugmap=plugmap1, loglam=loglam1, $
       wave=wave1, tsobj=tsobj1, zans=zans1, synflux=synflux1, $
       znum=znum, root_dir=root_dir, silent=silent

      if (ifile EQ 0) then begin
         allindx = indx
         if (q_flux) then flux = flux1
         if (q_flerr) then flerr = flerr1
         if (q_invvar) then invvar = invvar1
         if (q_andmask) then andmask = andmask1
         if (q_ormask) then ormask = ormask1
         if (q_disp) then disp = disp1
         if (q_plugmap) then plugmap = plugmap1
         if (q_loglam) then loglam = loglam1
         if (q_wave) then wave = wave1
         if (q_tsobj) then tsobj = tsobj1
         if (q_zans) then zans = zans1
         if (q_synflux) then synflux = synflux1
         if (q_mjd) then mjd = mjd1
      endif else begin
         allindx = [allindx, indx]
         if (q_flux) then spec_append, flux, flux1
         if (q_flerr) then spec_append, flerr, flerr1
         if (q_invvar) then spec_append, invvar, invvar1
         if (q_andmask) then spec_append, andmask, andmask1
         if (q_ormask) then spec_append, ormask, ormask1
         if (q_disp) then spec_append, disp, disp1
         if (q_plugmap) then plugmap = struct_append(plugmap, [plugmap1])
         if (q_loglam) then spec_append, loglam, loglam1
         if (q_wave) then spec_append, wave, wave1
         if (q_tsobj) then tsobj = struct_append(tsobj, [tsobj1])
         if (q_zans) then zans = struct_append(zans, [zans1])
         if (q_synflux) then spec_append, synflux, synflux1
         if (q_mjd) then mjd = [mjd, mjd1]
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
   if (q_synflux) then synflux[*,[allindx]] = synflux[*]
   if (q_mjd) then mjd[allindx] = mjd[*]

   if (keyword_set(silent)) then print

   return
end
;------------------------------------------------------------------------------
