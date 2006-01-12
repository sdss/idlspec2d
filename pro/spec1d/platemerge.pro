;+
; NAME:
;   platemerge
;
; PURPOSE:
;   Merge all Spectro-1D outputs with tsObj files.
;
; CALLING SEQUENCE:
;   platemerge, [ plate, mjd=, outroot=, public= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   plate       - Plates to include; default to all files
;                 specified by the PLATELIST routine.
;   mjd         - Optional MJDs corresponding to the specified PLATEs;
;                 if specified, then PLATE and MJD should have the same
;                 number of elements.
;   outroot     - Root name for output files; default to '$SPECTRO_DATA/spAll';
;                 the files are then 'spAll.fits', 'spAll.dat', 'spAllLine.dat'.
;                 If /PUBLIC is set, then add '-public' to the root name.
;                 If PUBLIC is set to a string, then add that string to the
;                 root name.
;   public      - If set with /PUBLIC, then limit to plates that have
;                 any entry in the PUBLIC field of the plate list.
;                 If set to a string, then select those plates that contain
;                 the substring PUBLIC within their PUBLIC field.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The SPECPRIMARY output element is used to select a unique set of
;   objects in the case of duplicate observations.  Any objects observed
;   multiple times will have SPECPRIMARY=1 for one instance only, and =0
;   for all other instances.  The criteria (in order of importance) are
;   as follows:
;     1) Prefer observations with positive SN_MEDIAN
;     2) Prefer PLATEQUALITY='good' over any other plate quality
;     3) Prefer observations with ZWARNING=0
;     4) Prefer objects with larger SN_MEDIAN
;
;   Temporary files are created first, such as 'spAll.fits.tmp', which
;   are renamed at the end of the routine to 'spAll.fits', etc.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   copy_struct
;   copy_struct_inx
;   djs_filepath()
;   headfits()
;   hogg_mrdfits()
;   mrdfits()
;   mwrfits_chunks
;   plug2tsobj()
;   platelist
;   readspec
;   repstr
;   spheregroup
;   splog
;   struct_print
;   sxaddpar
;   sxpar()
;
; REVISION HISTORY:
;   30-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platemerge, plate, mjd=mjd, outroot=outroot1, public=public

   dtheta = 2.0 / 3600.

   if (keyword_set(outroot1)) then begin
      outroot = [outroot1, outroot1+'Line']
   endif else begin
      outroot = ['spAll','spAllLine']
      if (keyword_set(public)) then begin
         if (size(public,/tname) EQ 'STRING') then $
          outroot = outroot + '-' + public $
          else outroot = outroot + '-public'
      endif
      outroot = djs_filepath(outroot, root_dir=getenv('SPECTRO_DATA'))
   endelse

   t1 = systime(1)
   thismem = memory()

   ;----------
   ; Find the list of spZ files.

   platelist, plist=plist
   if (NOT keyword_set(plist)) then return

   if (keyword_set(plate)) then begin
      nplate = n_elements(plist)
      if (keyword_set(mjd) AND n_elements(mjd) NE nplate) then $
       message, 'Number of elements in PLATE and MJD must agree'

      qkeep = bytarr(nplate)
      if (keyword_set(mjd)) then begin
         for i=0L, n_elements(plate)-1 do $
          qkeep = qkeep OR plist.plate EQ plate[i]
      endif else begin
         for i=0L, n_elements(plate)-1 do $
          qkeep = qkeep OR (plist.plate EQ plate[i] AND plist.mjd EQ mjd[i])
      endelse
      ikeep = where(qkeep, nkeep)
      if (nkeep EQ 0) then return
      plist = plist[ikeep]
   endif

   indx = where(strtrim(plist.status1d,2) EQ 'Done' AND $
    (strtrim(plist.platequality,2) EQ 'good' $
    OR strtrim(plist.platequality,2) EQ 'marginal' $
    OR strtrim(plist.public,2) NE ''), ct)
   if (ct EQ 0) then return
   if (keyword_set(public)) then begin
      if (size(public,/tname) EQ 'STRING') then begin
         itrim = where(strmatch(plist[indx].public,'*'+public+'*'), ntrim)
      endif else begin
         itrim = where(strtrim(plist[indx].public) NE '', ntrim)
      endelse
      if (ntrim EQ 0) then return
      indx = indx[itrim]
   endif
   plist = plist[indx]

   nfile = n_elements(plist)
   fullzfile = strarr(nfile)
   fullzfile = 'spZbest-' + string(plist.plate, format='(i4.4)') $
    + '-' + string(plist.mjd, format='(i5.5)') + '.fits'
   zsubdir = string(plist.plate, format='(i4.4)')
   for i=0L, nfile-1 do $
    fullzfile[i] = djs_filepath(fullzfile[i], $
     root_dir=getenv('SPECTRO_DATA'), subdirectory=zsubdir[i])

   splog, 'Found ', nfile, ' files'
   if (nfile EQ 0) then return
   fullzfile = fullzfile[ sort(fullzfile) ]

   nout = nfile * 640L
   splog, 'Total number of objects = ', nout

   ;----------
   ; Find the corresponding spPlate files (needed only for PRIMTARGET+SECTARGET
   ; flags, which are incorrect in the tsObj files.

   fullplatefile = repstr(fullzfile, 'spZbest', 'spPlate')

   ;----------
   ; Find the first tsObj file that exists for use in constructing the
   ; output structure.

   ifile = 0
   while (NOT keyword_set(tsobj0)) do begin
      tsobj0 = plug2tsobj(plist[ifile].plate, 0, 0)
      ifile = ifile + 1
      if (ifile EQ nfile) then $
       message, 'No calibObj files found!'
   endwhile

   ;----------
   ; Create the additional tags to add to the output structure

   pstuff = create_struct( $
    'progname'    , ' ', $
    'chunkname'   , ' ', $
    'platequality', ' ', $
    'platesn2'    , 0.0, $
    'primtarget'  ,  0L, $
    'sectarget'   ,  0L, $
    'specprimary' ,  0B, $
    'specobj_id'  ,  0L, $
    'nspecobs'    ,  0, $
    'calibflux'   , fltarr(5), $
    'calibflux_ivar', fltarr(5) )

   ;----------
   ; Loop through each file

   splog, 'Reading ZANS files'
   for ifile=0L, nfile-1 do begin
      print, 'Reading ZANS file ',ifile+1, ' of ', nfile

      hdr = headfits(fullzfile[ifile])
      plate = sxpar(hdr, 'PLATEID')
      zans = mrdfits(fullzfile[ifile], 1, /silent)
      plugmap = mrdfits(fullplatefile[ifile], 5, /silent)

      if (ifile EQ 0) then begin
         outdat1 = create_struct(pstuff, zans[0])
         struct_assign, {junk:0}, outdat1 ; Zero-out all elements
         outdat = replicate(outdat1, nout)
      endif

      indx = (ifile * 640L) + lindgen(640)

      ; The following is very slow, so we do this differently...
;      copy_struct_inx, zans, outdat, index_to=indx
      tmpdat = outdat[indx]
      copy_struct, zans, tmpdat
      outdat[indx] = tmpdat

      ; Fill in the first columns of this output structure
      outdat[indx].progname = plist[ifile].progname
      outdat[indx].chunkname = plist[ifile].chunkname
      outdat[indx].platequality = plist[ifile].platequality
      outdat[indx].platesn2 = plist[ifile].platesn2

      ; Get PRIMTARGET+SECTARGET with those values from
      ; the plug-map structure in spPlate file.
      outdat[indx].primtarget = plugmap.primtarget
      outdat[indx].sectarget = plugmap.sectarget

      ; Read the following from the plug-map if those tags exist
      if (tag_exist(plugmap,'CALIBFLUX')) then $
       outdat[indx].calibflux = plugmap.calibflux
      if (tag_exist(plugmap,'CALIBFLUX_IVAR')) then $
       outdat[indx].calibflux_ivar = plugmap.calibflux_ivar
   endfor

   splog, 'Time to read data = ', systime(1)-t1, ' sec'

   ;----------
   ; Set the SPECPRIMARY flag to 0 or 1

   t2 = systime(1)

   ; Determine the score for each object
   ; 1) Prefer observations with positive SN_MEDIAN
   ; 2) Prefer PLATEQUALITY='good' over any other plate quality
   ; 3) Prefer observations with ZWARNING=0
   ; 4) Prefer objects with larger SN_MEDIAN
   score = 4 * (outdat.sn_median GT 0) $
    + 2 * (strmatch(outdat.platequality,'good*') EQ 1) $
    + 1 * (outdat.zwarning EQ 0) $
    + (outdat.sn_median>0) / max(outdat.sn_median+1.)

   ingroup = spheregroup(outdat.plug_ra, outdat.plug_dec, dtheta, $
    multgroup=multgroup, firstgroup=firstgroup, nextgroup=nextgroup)

   ; Set the unique object IDs
   outdat.specobj_id = ingroup + 1L

   for j=0L, n_elements(firstgroup)-1L do begin
      if (firstgroup[j] NE -1) then begin
         if (multgroup[j] EQ 1) then begin
            outdat[firstgroup[j]].specprimary = 1
            outdat[firstgroup[j]].nspecobs = 1
         endif else begin
            indx = lonarr(multgroup[j])
            indx[0] = firstgroup[j]
            for k=0L, multgroup[j]-2L do indx[k+1] = nextgroup[indx[k]]
            foo = max(score[indx], ibest)
            outdat[indx[ibest]].specprimary = 1
            outdat[indx].nspecobs = multgroup[j]
         endelse
      endif
   endfor

   splog, 'Time to assign primaries = ', systime(1)-t2, ' sec'

   ;----------
   ; Pre-condition to FITS structure to have same-length strings
   ; (for any given tag name) by concatenating spaces.

   ntag = n_tags(outdat)
   tags = tag_names(outdat)
   for itag=0L, ntag-1L do begin
      if (size(outdat[0].(itag), /tname) EQ 'STRING') then begin
         if (NOT keyword_set(silent)) then $
          print, 'Padding whitespace for string array ' + tags[itag]
         taglen = strlen(outdat.(itag))
         maxlen = max(taglen)
         padspace = string('', format='(a'+string(maxlen)+')')
         outdat.(itag) = strmid(outdat.(itag) + padspace, 0, maxlen)
      endif
   endfor

   ;----------
   ; Write the output FITS file, writing one plate at a time

   ; Don't allow duplicate tags between the tsObj structure and what
   ; is already in the output structure.  For ex, MJD is in both.
   tsobj0 = struct_selecttags(tsobj0, except_tags=tag_names(outdat))
   platedat1 = create_struct(outdat[0], tsobj0)
   struct_assign, {junk:0}, platedat1 ; Zero-out all elements

   splog, 'Writing FITS file ' + outroot[0]+'.fits'
   for ifile=0L, nfile-1 do begin
      print, 'Writing plate ', ifile+1, ' of ', nfile

      platedat = replicate(platedat1, 640)
      indx = (ifile * 640L) + lindgen(640)
      tsobj = plug2tsobj(plist[ifile].plate, outdat[indx].plug_ra, $
       outdat[indx].plug_dec)
      if (keyword_set(tsobj)) then $
       copy_struct_inx, tsobj, platedat $
      else $
       splog, 'WARNING: No tsObj file found for plate ', outdat[indx[0]].plate
;      copy_struct_inx, outdat, platedat, index_from=indx
      copy_struct, outdat[indx], platedat

      ; All strings must be the same length, or appending to the FITS file
      ; will result in corruption.  The only string in the tsobj structure
      ; is for the RERUN.
      platedat.rerun = string(platedat.rerun+'   ',format='(a3)')

      mwrfits_chunks, platedat, outroot[0]+'.fits.tmp', $
       create=(ifile EQ 0), append=(ifile GT 0)
   endfor

   outdat = 0 ; Clear memory

   ;----------
   ; Create the structure for ASCII output

   adat1 = create_struct( $
    'plate'      ,  0L, $
    'mjd'        ,  0L, $
    'fiberid'    ,  0L, $
    'class'      ,  '', $
    'subclass'   ,  '', $
    'z'          , 0.0, $
    'z_err'      , 0.0, $
    'zwarning'   ,  0L, $
    'rchi2'      , 0.0, $
    'plug_ra'    , 0.0d, $
    'plug_dec'   , 0.0d, $
    'platesn2'   ,  0.0, $
    'modelflux'  ,  fltarr(5), $
    'objc_type'  ,  '', $
    'primtarget' ,  0L, $
    'sectarget'  ,  0L, $
    'progname',     '', $
    'specprimary',  0L, $
    'objtype'    ,  '' )

   ascii_tags = [ $
    'plate', 'mjd', 'fiberid', 'class', 'subclass', 'z', 'z_err', $
    'zwarning', 'rchi2', 'plug_ra', 'plug_dec', 'platesn2', $
    'modelflux', 'objc_type', 'primtarget', 'sectarget', 'progname', $
    'specprimary']

   ; Read the tags that we need from the FITS file
   outdat = hogg_mrdfits(outroot[0]+'.fits.tmp', 1, nrowchunk=10000L, $
    columns=ascii_tags)
   adat = replicate(adat1, n_elements(outdat))
   copy_struct, outdat, adat

   ; Replace any blank strings for CLASS with "".
   ii = where(strtrim(adat.class,2) EQ '')
   if (ii[0] NE -1) then adat[ii].class = '""'

   ; Replace any blank strings for SUBCLASS with "".
   ; If SUBCLASS contains several words, then use a plus sign between
   ; the words rather than a space.
   adat.subclass = strtrim(adat.subclass,2)
   ii = where(adat.subclass EQ '')
   if (ii[0] NE -1) then adat[ii].subclass = '""'
   adat.subclass = repstr(adat.subclass, ' ', '+')

   objtypes = ['UNKNOWN', 'CR', 'DEFECT', 'GALAXY', 'GHOST', 'KNOWNOBJ', $
    'STAR', 'TRAIL', 'SKY']
   adat.objc_type = objtypes[outdat.objc_type]

   outdat = 0 ; Clear memory

   splog, 'Writing ASCII file ' + outroot[0]+'.dat'
   struct_print, adat, filename=outroot[0]+'.dat.tmp'

   adat = 0 ; Clear memory

   ;----------
   ; Create the merged line data

   splog, 'Writing FITS zline file ' + outroot[1]+'.fits'
   for ifile=0L, nfile-1 do begin
      splog, 'Writing zline ', ifile+1, ' of ', nfile
      readspec, plist[ifile].plate, mjd=plist[ifile].mjd, zline=linedat

      if (ifile EQ 0) then begin
         nobj = nfile * 640L
         nper = n_elements(linedat) / 640L
         sxaddpar, linehdr, 'DIMS0', nper, ' Number of emission lines'
         sxaddpar, linehdr, 'DIMS1', nobj, ' Number of objects'
         linedat1 = linedat[0]
         struct_assign, {junk:0}, linedat1
      endif

      ; Demand that the structure has the same format as the first
      ; one written.
      linedat_out = replicate(linedat1, n_elements(linedat))
      struct_assign, linedat, linedat_out

      mwrfits_chunks, linedat_out, outroot[1]+'.fits.tmp', $
       create=(ifile EQ 0), append=(ifile GT 0)
   endfor

   ;----------
   ; Rename temporary files

   spawn, ['mv', outroot[0]+'.fits.tmp', outroot[0]+'.fits'], /noshell
   spawn, ['mv', outroot[0]+'.dat.tmp', outroot[0]+'.dat'], /noshell
   spawn, ['mv', outroot[1]+'.fits.tmp', outroot[1]+'.fits'], /noshell

   thismem = memory()
   maxmem = thismem[3]
   splog, 'Maximum memory usage = ', maxmem/1.d6, ' MB'
   splog, 'Total time = ', systime(1)-t1, ' sec'

   return
end
;------------------------------------------------------------------------------
