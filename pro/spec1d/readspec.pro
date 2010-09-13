;+
; NAME:
;   readspec
;
; PURPOSE:
;   Routine for reading SDSS/BOSS Spectro-2D and Spectro-1D files
;
; CALLING SEQUENCE:
;   readspec, plate, fiber, [mjd=, znum=, flux=, flerr=, invvar=, $
;    andmask=, ormask=, disp=, sky=, plugmap=, loglam=, wave=, tsobj=, $
;    zans=, zmanual=, zline=, synflux=, lineflux=, objhdr=, zhdr=, nfiber=, $
;    topdir=, path=, /align, /silent ]
;
; INPUTS:
;   plate      - Plate number(s)
;
; OPTIONAL INPUTS:
;   fiber      - Fiber number(s), 1-indexed; if not set, or zero, then
;                read all fibers for each plate.
;   mjd        - MJD number(s); if not set, then select the most recent
;                data for this plate (largest MJD).
;   znum       - If set, then return not the best-fit redshift, but the
;                ZNUM-th best-fit; e.g., set ZNUM=2 for second-best fit.
;   topdir     - Optional override value for the environment
;                variable $BOSS_SPECTRO_REDUX.
;   path       - Override all path information with this directory name.
;   align      - If set, then align all the spectra in wavelength.
;                Also, LOGLAM and WAVE will be output as single vectors
;                (since they are all the same) rather than as one per object.
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
;   sky        - Sky flux [NPIXEL,NFIBER]
;   plugmap    - Plug-map entries [NFIBER]
;   loglam     - Log10-wavelength in log10-Angstroms [NPIXEL,NFIBER],
;                or the vector [NPIXEL] if /ALIGN is set
;   wave       - Wavelength in Angstroms [NPIXEL,NFIBER],
;                or the vector [NPIXEL] if /ALIGN is set
;   tsobj      - tsObj-structure output [NFIBER]
;   zans       - Redshift output structure [NFIBER]
;   zmanual    - Manual inspection structure [NFIBER]
;   zline      - Line-fit output structure [NFIBER,NLINE]
;   synflux    - Best-fit synthetic eigen-spectrum [NPIXEL,NFIBER]
;   lineflux   - Best-fit emission line fits + background terms  [NPIXEL,NFIBER]
;   objhdr     - The FITS header from the first object spPlate file read.
;                If spectra from multiple plates are read, then it is
;                indeterminant which header this will be.
;   zhdr       - The FITS header from the first object spZ file read.
;                If spectra from multiple plates are read, then it is
;                indeterminant which header this will be.
;   nfiber     - Number of fibers per plate [NFIBER]
;
; COMMENTS:
;   One can input PLATE and FIBER as vectors, in which case there must
;   be a one-to-one correspondence between them.  Or, one can input FIBER
;   numbers as a vector, in which case the same PLATE is used for all.
;   Or, one can input PLATE as a vector, in which case the same FIBER is
;   read for all.
;
;   The environment variable BOSS_SPECTRO_REDUX must be set to tell this routine
;   where to find the data.  The reduced spectro data files are assumed to
;   be $BOSS_SPECTRO_REDUX/$RUN2D/pppp/spPlate-pppp-mmmmm.fits,
;   where pppp=plate number and mmmm=MJD.
;
;   The tsObj files are assumed to be in the directory $BOSS_SPECTRO_REDUX/plates.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $BOSS_SPECTRO_REDUX/$RUN2D/$PLATE/spPlate-$PLATE-$MJD.fits
;   $BOSS_SPECTRO_REDUX/$RUN2D/$PLATE/$RUN1D/spZbest-$PLATE-$MJD.fits
;   $BOSS_SPECTRO_REDUX/$RUN2D/$PLATE/$RUN1D/spZall-$PLATE-$MJD.fits
;   $BOSS_SPECTRO_REDUX/$RUN2D/$PLATE/$RUN1D/spZline-$PLATE-$MJD.fits
;   $BOSS_SPECTRO_REDUX/plates/tsObj*-$PLATE.fit
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
;   rspec_zline_append()
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
 andmask=andmask, ormask=ormask, disp=disp, sky=sky, plugmap=plugmap, $
 loglam=loglam, tsobj=tsobj, zans=zans, zmanual=zmanual, zline=zline, $
 synflux=synflux, lineflux=lineflux, znum=znum, objhdr=objhdr, zhdr=zhdr, $
 nfiber=nfiber, coeffzero=coeff0, coeffone=coeff1, npix=npix, $
 topdir=topdir1, run2d=run2d, run1d=run1d, path=path, $
 align=align, silent=silent, qread=qread, sdss=sdss

   platestr = string(plate,format='(i4.4)')
   if (NOT keyword_set(mjd)) then mjdstr = '*' $
    else mjdstr = string(mjd,format='(i5.5)')
   if (keyword_set(path)) then begin
      topdir = path
   endif else begin
      if keyword_set(topdir1) then topdir = topdir1[0] else $
          topdir = keyword_set(sdss) ? getenv('SPECTRO_REDUX') : getenv('BOSS_SPECTRO_REDUX')
      ; topdir = keyword_set(topdir1) ? topdir1[0] : getenv('BOSS_SPECTRO_REDUX')
      oneddir = n_elements(run1d) GT 0 ? strtrim(run1d[0],2) : getenv('RUN1D')
      twoddir = n_elements(run2d) GT 0 ? strtrim(run2d[0],2) : getenv('RUN2D')
   endelse


   filename = 'spPlate-' + platestr + '-' + mjdstr + '.fits'
   if (keyword_set(path)) then $
    filename = lookforgzip(filepath(filename, root_dir=path), count=ct) $
   else $
    filename = lookforgzip(filepath(filename, root_dir=topdir, $
     subdirectory=[twoddir,platestr]), count=ct)

   if (ct GT 1) then filename = filename[ (reverse(sort(filename)))[0] ] $
    else filename = filename[0]

   nrows = n_elements(rownums)

   ; Set default return values
   flux = fltarr(1,nrows)
   flerr = fltarr(1,nrows)
   invvar = fltarr(1,nrows)
   andmask = lonarr(1,nrows)
   ormask = lonarr(1,nrows)
   disp = fltarr(1,nrows)
   sky = fltarr(1,nrows)
   plugmap = lonarr(nrows)
   loglam = fltarr(1,nrows)
   tsobj = lonarr(nrows)
   zans = lonarr(nrows)
   zmanual = lonarr(nrows)
   zline = lonarr(nrows)
   synflux = fltarr(1,nrows)
   lineflux = fltarr(1,nrows)
   coeff0 = 0
   coeff1 = 0
   npix = 0
   nfiber = (mjd LT 55025) ? replicate(640L,nrows) : replicate(1000L,nrows)

   if (NOT keyword_set(filename)) then begin
      return
   end

   fits_open, filename, fcb

   if (qread.flux) then begin
      flux = rspec_mrdfits(fcb, 0, rownums=rownums, silent=silent)
   endif

   if (qread.invvar OR qread.flerr) then begin
      invvar = rspec_mrdfits(fcb, 1, rownums=rownums, silent=silent)
      if (qread.flerr) then begin
         i = where(invvar GT 0)
         flerr = 0 * invvar
         if (i[0] NE -1) then flerr[i] = 1 / sqrt(invvar[i])
      endif
   endif

   if (qread.andmask) then begin
      andmask = rspec_mrdfits(fcb, 2, rownums=rownums, silent=silent)
   endif

   if (qread.ormask) then begin
      ormask = rspec_mrdfits(fcb, 3, rownums=rownums, silent=silent)
   endif

   if (qread.disp) then begin
      disp = rspec_mrdfits(fcb, 4, rownums=rownums, silent=silent)
   endif

   if (qread.sky) then begin
      sky = rspec_mrdfits(fcb, 6, rownums=rownums, silent=silent)
   endif

   if (qread.plugmap) then begin
      plugmap = rspec_mrdfits(fcb, 5, rownums=rownums, silent=silent)
   endif

   if (qread.needwave) then begin
      if (NOT keyword_set(objhdr)) then objhdr = headfits(filename)
      coeff0 = replicate(sxpar(objhdr, 'COEFF0'), nrows)
      coeff1 = replicate(sxpar(objhdr, 'COEFF1'), nrows)
      npix = sxpar(objhdr, 'NAXIS1')
   endif

   if (qread.tsobj) then begin
      if keyword_set(sdss) then $
         matchdir=getenv('SPECTRO_MATCH')+'/'+twodir+'/'+ $
            file_basename(getenv('PHOTO_RESOLVE'))+'/'+platestr $
      else $
         matchdir=topdir+'/'+twoddir+'/'+platestr
      tsobj1 = plug2tsobj(plate, mjd=mjd, $
       indir=matchdir, silent=silent)
      if (keyword_set(tsobj1)) then tsobj = tsobj1[rownums]
   endif

   if (qread.zmanual) then begin
      zmanual1 = spmanual(plate, mjd=mjd)
      if (keyword_set(zmanual1)) then zmanual = zmanual1[rownums]
   endif

   if (qread.zans OR qread.synflux OR qread.zhdr) then begin
      if (NOT keyword_set(znum)) then $
       zfile = 'spZbest-' + platestr + '-' + mjdstr + '.fits' $
      else $
       zfile = 'spZall-' + platestr + '-' + mjdstr + '.fits'

      if (keyword_set(path)) then $
       zfile = lookforgzip(filepath(zfile, root_dir=path), count=ct) $
      else $
       zfile = lookforgzip(filepath(zfile, root_dir=topdir, $
        subdirectory=[twoddir,platestr,oneddir]), count=ct)
      if (ct GT 1) then zfile = zfile[ (reverse(sort(zfile)))[0] ] $
       else zfile = zfile[0]

      if (keyword_set(zfile)) then begin
         if (NOT keyword_set(znum)) then begin
            fits_open, zfile, zfcb
            zans = rspec_mrdfits(zfcb, 1, rownums=rownums, silent=silent)
            fits_close, zfcb
         endif else begin
            zhdr0 = headfits(zfile, exten=0)
            zhdr = headfits(zfile, exten=1)
            nper = sxpar(zhdr0, 'DIMS0') ; number of fits per object
            fits_open, zfile, zfcb
            zans = rspec_mrdfits(zfcb, 1, rownums=rownums*nper+znum-1, $
             silent=silent)
            fits_close, zfcb
         endelse
      endif
   endif

   if (qread.synflux AND keyword_set(zfile)) then begin
      ; Read the synthetic spectrum from the Zbest file if ZNUM is not set.
      if (NOT keyword_set(znum)) then begin
         fits_open, zfile, zfcb
         synflux = rspec_mrdfits(zfcb, 2, rownums=rownums, silent=silent)
         fits_close, zfcb
      endif else begin
         if (NOT keyword_set(objhdr)) then objhdr = headfits(filename)
         if (keyword_set(zans)) then $
          synflux = synthspec(zans, hdr=objhdr)
      endelse
   endif

   if (qread.zline OR qread.lineflux) then begin
      linefile = 'spZline-' + platestr + '-' + mjdstr + '.fits'

      if (keyword_set(path)) then $
       linefile = lookforgzip(filepath(linefile, root_dir=path), count=ct) $
      else $
       linefile = lookforgzip(filepath(linefile, root_dir=topdir, $
        subdirectory=[twoddir,platestr,oneddir]), count=ct)
      if (ct GT 1) then linefile = linefile[ (reverse(sort(linefile)))[0] ] $
       else linefile = linefile[0]
   endif

   if (qread.zline AND keyword_set(linefile)) then begin
      linehdr = headfits(linefile)
      nlines = sxpar(linehdr, 'DIMS0') ; number of emission lines per object

      fits_open, linefile, linefcb
      allrows = reform( rebin(reform(rownums*nlines,1,nrows), nlines, nrows), $
       nlines*nrows ) $
       + reform( rebin(lindgen(nlines), nlines, nrows), nlines*nrows)
      allrows = (rebin(reform(rownums*nlines,1,nrows), nlines, nrows))[*] $
       + (rebin(lindgen(nlines), nlines, nrows))[*]
      zline = rspec_mrdfits(linefcb, 1, $
       rownums=allrows, silent=silent)
      fits_close, linefcb

      zline = reform(zline, nlines, nrows)
   endif

   if (qread.lineflux AND keyword_set(linefile)) then begin
      ; Read the line fit flux from the Zline file
      fits_open, linefile, linefcb
      lineflux = rspec_mrdfits(linefcb, 2, rownums=rownums, silent=silent)
      fits_close, linefcb
   endif

   if (qread.mjd) then begin
      if (NOT keyword_set(objhdr)) then objhdr = headfits(filename)
      mjd = sxpar(objhdr, 'MJD')
   endif

   fits_close, fcb

   if (qread.objhdr AND (NOT keyword_set(objhdr))) then objhdr = headfits(filename)
   if (qread.zhdr AND (NOT keyword_set(zhdr))) then zhdr = headfits(zfile)
   if (qread.nfiber) then begin
      if (keyword_set(objhdr)) then nfiber = long(sxpar(objhdr, 'NAXIS2')) $
       else nfiber = 0L
      nfiber = replicate(nfiber, nrows)
   endif

   return
end

;------------------------------------------------------------------------------
; Append two ZLINE structures, even if they are different structures
; or have different number of lines measured.
; If the structure definitions are different, use the definition from ZLINE1.
; If the number of lines are different, use the larger line list, and set
; those lines to all zeros in the shorter line list.

function rspec_zline_append, zline1, zline2

   qstruct1 = size(zline1,/tname) EQ 'STRUCT'
   qstruct2 = size(zline2,/tname) EQ 'STRUCT'

   if (qstruct1 EQ 0 AND qstruct2 EQ 0) then $
    return, lonarr(n_elements(zline1) + n_elements(zline2))

   ndim1 = size(zline1, /n_dimen)
   ndim2 = size(zline2, /n_dimen)
   dims1 = size(zline1, /dimens)
   dims2 = size(zline2, /dimens)
   nline = 1L
   if (qstruct1) then begin
      if (ndim1 EQ 1) then nobj1 = 1L $
       else nobj1 = dims1[1]
      nline = nline > dims1[0]
   endif else begin
      nobj1 = n_elements(zline1)
   endelse
   if (qstruct2) then begin
      if (ndim2 EQ 1) then obj2 = 1L $
       else nobj2 = dims2[1]
      nline = nline > dims2[0]
   endif else begin
      nobj2 = n_elements(zline2)
   endelse
   nobjtot = nobj1 + nobj2

   if (qstruct1) then blankline = zline1[0] $
    else blankline = zline2[0]
   struct_assign, {junk:0}, blankline
   zlinetot = replicate(blankline, nline, nobjtot)
   if (qstruct1) then $
    copy_struct_inx, zline1, zlinetot, $
     index_to=(lindgen(dims1[0]) # (lonarr(nobj1)+1) $
     + (lonarr(dims1[0])+1) # (lindgen(nobj1)*nline))[*]
   if (qstruct2) then $
    copy_struct_inx, zline2, zlinetot, $
     index_to=(lindgen(dims2[0]) # (lonarr(nobj2)+1) $
     + (lonarr(dims2[0])+1) # (lindgen(nobj2)*nline))[*] + nline*nobj1

   return, zlinetot
end

;------------------------------------------------------------------------------
pro readspec, plate, fiber, mjd=mjd, flux=flux, flerr=flerr, invvar=invvar, $
 andmask=andmask, ormask=ormask, disp=disp, sky=sky, plugmap=plugmap, $
 loglam=loglam, wave=wave, tsobj=tsobj, zans=zans, zmanual=zmanual, $
 zline=zline, synflux=synflux, lineflux=lineflux, objhdr=objhdr, zhdr=zhdr, $
 nfiber=nfiber, znum=znum, align=align, silent=silent, _EXTRA=Extra

   if (n_params() LT 1) then begin
      doc_library, 'readspec'
      return
   endif

   ; This common block specifies which keywords will be returned.
   qread_blank = { flux: 0b, flerr: 0b, invvar: 0b, andmask: 0b, ormask: 0b, $
    disp: 0b, sky: 0b, plugmap: 0b, loglam: 0b, wave: 0b, tsobj: 0b, $
    zans: 0b, zmanual: 0b, zline: 0b, synflux: 0b, lineflux: 0b, mjd: 0b, $
    objhdr: 0b, zhdr: 0b, nfiber: 0b, needwave: 0b }

; ???
;   if (keyword_set(topdir) EQ 0 AND keyword_set(path) EQ 0) then begin
;      topdir = getenv('BOSS_SPECTRO_REDUX')
;      if (NOT keyword_set(topdir)) then $
;       message, 'Environment variable BOSS_SPECTRO_REDUX must be set!'
;   endif

   qread = qread_blank
   qread.flux = arg_present(flux)
   qread.flerr = arg_present(flerr)
   qread.invvar = arg_present(invvar)
   qread.andmask = arg_present(andmask)
   qread.ormask = arg_present(ormask)
   qread.disp = arg_present(disp)
   qread.sky = arg_present(sky)
   qread.plugmap = arg_present(plugmap) OR arg_present(tsobj)
   qread.loglam = arg_present(loglam)
   qread.wave = arg_present(wave)
   qread.tsobj = arg_present(tsobj)
   qread.zans = arg_present(zans)
   qread.zmanual = arg_present(zmanual)
   qread.zline = arg_present(zline)
   qread.synflux = arg_present(synflux)
   qread.lineflux = arg_present(lineflux)
   qread.mjd = arg_present(mjd) AND (keyword_set(mjd) EQ 0)
   qread.objhdr = arg_present(objhdr) OR (arg_present(nfiber))
   qread.zhdr = arg_present(zhdr)
   qread.nfiber = arg_present(nfiber)
   objhdr = ''
   zhdr = ''
   qread.needwave = qread.loglam OR qread.wave OR keyword_set(align)

   nplate = n_elements(plate)
   if (nplate EQ 0) then $
    message, 'PLATE must be defined'
   if (keyword_set(mjd) AND n_elements(mjd) NE nplate) then $
    message, 'Number of elements in PLATE and MJD must agree'

   if (NOT keyword_set(fiber)) then begin
      ; Special case to read all fibers of each plate
      ; We need to know how many there are on each plate!
      ; Find MJD if not specified...
      if (keyword_set(mjd)) then mjd_tmp=mjd $
       else readspec, plate, plate*0+1, mjd=mjd_tmp, $
        silent=silent, _EXTRA=Extra
      readspec, plate, plate*0+1, mjd=mjd_tmp, nfiber=nfiber_tmp, $
       silent=silent, _EXTRA=Extra
      nfiber_tot = long(total(nfiber_tmp))
      platevec = lonarr(nfiber_tot>1)
      fibervec = lonarr(nfiber_tot>1)
      mjdvec = lonarr(nfiber_tot>1)
      j = 0L
      for i=0L, nplate-1L do begin
         if (nfiber_tmp[i] GT 0) then begin
            platevec[j:j+nfiber_tmp[i]-1] = plate[i]
            fibervec[j:j+nfiber_tmp[i]-1] = lindgen(nfiber_tmp[i]) + 1
            mjdvec[j:j+nfiber_tmp[i]-1] = mjd_tmp[i]
            j += nfiber_tmp[i]
         endif
      endfor
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

   for ifile=0L, nfile-1L do begin
      objhdr1 = 0
      zhdr1 = 0
      flux1 = 0
      flerr1 = 0
      invvar1 = 0
      andmask1 = 0
      ormask1 = 0
      disp1 = 0
      sky1 = 0
      plugmap1 = 0
      tsobj1 = 0
      zans1 = 0
      zmanual1 = 0
      zline1 = 0
      synflux1 = 0
      lineflux1 = 0
      nfiber1 = 0

      indx = where(platevec EQ platenums[ifile] AND mjdvec EQ mjdnums[ifile])
      irow = fibervec[indx] - 1

;      if (keyword_set(silent)) then print, '+', format='(A,$)'

      mjd1 = mjdnums[ifile]
      readspec1, platenums[ifile], irow, mjd=mjd1, $
       flux=flux1, flerr=flerr1, invvar=invvar1, andmask=andmask1, $
       ormask=ormask1, disp=disp1, sky=sky1, plugmap=plugmap1, $
       tsobj=tsobj1, zans=zans1, zmanual=zmanual1, zline=zline1, $
       synflux=synflux1, lineflux=lineflux1, objhdr=objhdr1, zhdr=zhdr1, $
       nfiber=nfiber1, znum=znum, coeffzero=coeff0, coeffone=coeff1, npix=npix, $
       align=align, silent=silent, qread=qread, _EXTRA=Extra
      coeff0 = double(coeff0)
      coeff1 = double(coeff1)

      if (qread.objhdr AND NOT keyword_set(objhdr)) then objhdr = objhdr1
      if (qread.zhdr AND NOT keyword_set(zhdr)) then zhdr = zhdr1
      if (ifile EQ 0) then begin
         allindx = indx
         if (qread.needwave) then begin
            allcoeff0 = coeff0
            allcoeff1 = coeff1
            npixmax = npix
         endif
         if (qread.flux) then flux = flux1
         if (qread.flerr) then flerr = flerr1
         if (qread.invvar) then invvar = invvar1
         if (qread.andmask) then andmask = andmask1
         if (qread.ormask) then ormask = ormask1
         if (qread.disp) then disp = disp1
         if (qread.sky) then sky = sky1
         if (qread.plugmap) then plugmap = plugmap1
         if (qread.tsobj) then tsobj = tsobj1
         if (qread.zans) then zans = zans1
         if (qread.zmanual) then zmanual = zmanual1
         if (qread.zline) then zline = zline1
         if (qread.synflux) then synflux = synflux1
         if (qread.lineflux) then lineflux = lineflux1
         if (qread.mjd) then mjd = mjd1
         if (qread.nfiber) then nfiber = nfiber1
      endif else begin
         allindx = [allindx, indx]
         if (qread.needwave) then begin
            ; If pixshift > 0, then this newly-read spectrum starts at
            ; bigger wavelengths than the previously-read spectra.
            ; Adjust the starting wavelengths appropriately to correspond
            ; to how spec_append will shift the spectra.
            mincoeff0 = min(allcoeff0)
            if (keyword_set(align)) then begin
               ; The following few lines of code deal with setting a
               ; wavelength scale for objects with missing data.
               if (mincoeff0 EQ 0 AND coeff0[0] GT 0) then begin
                  allcoeff0[*] = coeff0[0]
                  allcoeff1[*] = coeff1[0]
               endif
               if (mincoeff0 GT 0 AND coeff0[0] EQ 0) then begin
                  coeff0[*] = mincoeff0
                  coeff1[*] = allcoeff1[0]
               endif

               ; Add 0.5 below to take care of round-off errors
               pixshift = floor( (coeff0[0] - mincoeff0) / coeff1[0] + 0.5 )
               if (pixshift GT 0) then begin
                  coeff0 = coeff0 - pixshift * coeff1
                  npixmax = max([npixmax, npix+pixshift])
               endif else begin
                  allcoeff0 = allcoeff0 + pixshift * allcoeff1
                  npixmax = max([npixmax-pixshift, npix])
               endelse
            endif else begin
               npixmax = max([npix, npixmax])
               pixshift = 0
            endelse
            allcoeff0 = [allcoeff0, coeff0]
            allcoeff1 = [allcoeff1, coeff1]
         endif
         if (qread.flux) then spec_append, flux, flux1, pixshift
         if (qread.flerr) then spec_append, flerr, flerr1, pixshift
         if (qread.invvar) then spec_append, invvar, invvar1, pixshift
         if (qread.andmask) then spec_append, andmask, andmask1, pixshift
         if (qread.ormask) then spec_append, ormask, ormask1, pixshift
         if (qread.disp) then spec_append, disp, disp1, pixshift
         if (qread.sky) then spec_append, sky, sky1, pixshift
         if (qread.plugmap) then plugmap = struct_append(plugmap, plugmap1, /force)
         if (qread.tsobj) then tsobj = struct_append(tsobj, tsobj1, /force)
         if (qread.zans) then zans = struct_append(zans, zans1, /force)
         if (qread.zmanual) then zmanual = struct_append(zmanual, zmanual1, /force)
; The first two attempts below can fail if the ZLINE structure changes.
;         if (qread.zline) then zline = struct_append(zline, [zline1])
;         if (qread.zline) then zline = [[zline], [zline1]]
         if (qread.zline) then zline = rspec_zline_append(zline, zline1)
         if (qread.synflux) then spec_append, synflux, synflux1, pixshift
         if (qread.lineflux) then spec_append, lineflux, lineflux1, pixshift
         if (qread.mjd) then mjd = [mjd, mjd1]
         if (qread.nfiber) then nfiber = [nfiber, nfiber1]
      endelse
   endfor

   ; Re-sort the data
   if (qread.flux) then flux[*,[allindx]] = flux[*]
   if (qread.flerr) then flerr[*,[allindx]] = flerr[*]
   if (qread.invvar) then invvar[*,[allindx]] = invvar[*]
   if (qread.andmask) then andmask[*,[allindx]] = andmask[*]
   if (qread.ormask) then ormask[*,[allindx]] = ormask[*]
   if (qread.disp) then disp[*,[allindx]] = disp[*]
   if (qread.sky) then sky[*,[allindx]] = sky[*]
   if (qread.plugmap) then begin
      if (keyword_set(plugmap[0])) then $
       copy_struct_inx, plugmap, plugmap, index_to=allindx
   endif
   if (qread.needwave) then begin
      allcoeff0[allindx] = allcoeff0[*] ; Must de-reference before assignment
      allcoeff1[allindx] = allcoeff1[*] ; Must de-reference before assignment
   endif
   if (qread.tsobj) then begin
      if (keyword_set(tsobj[0])) then $
       copy_struct_inx, tsobj, tsobj, index_to=allindx
   endif
   if (qread.zans) then begin
      if (keyword_set(zans[0])) then $
       copy_struct_inx, zans, zans, index_to=allindx
   endif
   if (qread.zmanual) then begin
      if (keyword_set(zmanual[0])) then $
       copy_struct_inx, zmanual, zmanual, index_to=allindx
   endif
   if (qread.zline) then begin
      if (keyword_set(zline[0])) then begin
         ; Logically, we want to make the assignment ZLINE[*,ALLINDX] = ZLINE
         nlines = (size(zline, /dimen))[0]
         index_to = make_array(size=size(zline), /long)
         for iline=0L, nlines-1 do index_to[iline,*] = allindx * nlines + iline
         copy_struct_inx, zline, zline, index_to=index_to[*]
      endif
   endif
   if (qread.synflux) then synflux[*,[allindx]] = synflux[*]
   if (qread.lineflux) then lineflux[*,[allindx]] = lineflux[*]
   if (qread.mjd) then mjd[allindx] = mjd[*]
   if (qread.nfiber) then nfiber[allindx] = nfiber[*]

   ;----------
   ; Construct the output wavelength solutions

   if (qread.loglam OR qread.wave) then begin
      npixmax = npixmax > 1
      if (keyword_set(align)) then begin
         loglam = allcoeff0[0] + allcoeff1[0] * lindgen(npixmax)
      endif else begin
         nobj = n_elements(allcoeff0)
         loglam = dblarr(npixmax,nobj)
         for iobj=0L, nobj-1 do $
          loglam[*,iobj] = allcoeff0[iobj] + allcoeff1[iobj] * lindgen(npixmax)
      endelse
      if (qread.wave AND qread.loglam) then wave = 10^loglam $
       else if (qread.wave) then wave = temporary(10^loglam)
   endif

   return
end
;------------------------------------------------------------------------------
