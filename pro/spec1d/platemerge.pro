;+
; NAME:
;   platemerge
;
; PURPOSE:
;   Merge all Spectro-1D outputs with tsObj files.
;
; CALLING SEQUENCE:
;   platemerge, [ zfile, outroot=, public= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   zfile       - Redshift file(s) from spectro-1D; default to all files
;                 specified by the PLATELIST routine.
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
;     1) Prefer PROGNAME='main' over any other program names
;     2) Prefer PLATEQUALITY='good' over any other plate quality
;     3) Prefer observations with ZWARNING=0
;     4) Prefer the observation with the larger PLATESN2
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   djs_angle_match()
;   djs_diff_angle()
;   headfits()
;   mrdfits()
;   mwrfits
;   mwrfits_chunks
;   plug2tsobj()
;   platelist
;   readspec
;   repstr
;   struct_print
;   sxaddpar
;   sxpar()
;
; REVISION HISTORY:
;   30-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platemerge, zfile, outroot=outroot1, public=public

   dtheta = 2.0 / 3600.
   tags_exclude = ['FIRST*','ROSAT*','MATCHID']

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

   ;----------
   ; Find the list of spZ files.

   if (NOT keyword_set(zfile)) then begin
      platelist, plist=plist
      if (NOT keyword_set(plist)) then return

      indx = where(strtrim(plist.status1d,2) EQ 'Done' AND $
       (strtrim(plist.platequality,2) EQ 'good' $
       OR strtrim(plist.platequality,2) EQ 'marginal'), ct)
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
   endif else begin
      fullzfile = findfile(zfile, count=nfile)
   endelse

   print, 'Found ', nfile, ' files'
   if (nfile EQ 0) then return
   fullzfile = fullzfile[ sort(fullzfile) ]

   nout = nfile * 640L
   print, 'Total number of objects = ', nout

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
       message, 'No tsObj files found!'
   endwhile

   ;----------
   ; Loop through each file

   for ifile=0L, nfile-1 do begin
      print,'File ',ifile+1, ' of ', nfile,': '+fullzfile[ifile]

      hdr = headfits(fullzfile[ifile])
      plate = sxpar(hdr, 'PLATEID')
      zans = mrdfits(fullzfile[ifile], 1, /silent)
      tsobj = plug2tsobj(plate, zans.plug_ra, zans.plug_dec)
      if (NOT keyword_set(tsobj)) then $
       splog, 'WARNING: No tsObj file found for plate ', plate

      if (NOT keyword_set(outdat)) then begin
         pstuff = create_struct( $
          'progname'    , ' ', $
          'chunkname'   , ' ', $
          'platequality', ' ', $
          'platesn2'    , 0.0, $
          'smearuse'    , ' ', $
          'primtarget'  ,  0L, $
          'sectarget'   ,  0L, $
          'specprimary' ,  0B )
         ; Trim away tag names from the tsObj structure that
         ; is also in the ZANS structure.  At the moment (May 2003),
         ; the only such tag is "MJD".
         tsobj0 = struct_trimtags(tsobj0, except_tags=tag_names(zans[0]))
         tmpout = create_struct(pstuff, zans[0], tsobj0)

         ; Exclude tags we don't want
         tags = tag_names(tmpout)
         ntag = n_elements(tags)
         qkeep = bytarr(ntag) + 1B
         for itag=0, ntag-1 do begin
            for jtag=0, n_elements(tags_exclude)-1 do begin
               if (strmatch(tags[itag], tags_exclude[jtag])) then $
                qkeep[itag] = 0B
            endfor
         endfor
         ikeep = where(qkeep, nkeep)
         outdat1 = create_struct(tags[ikeep[0]], tmpout.(ikeep[0]))
         for ii=1, nkeep-1 do $
          outdat1 = create_struct(outdat1, tags[ikeep[ii]], tmpout.(ikeep[ii]))

         struct_assign, {junk:0}, outdat1 ; Zero-out all elements
         sz1 = n_tags(outdat1, /length)
         splog, 'Size of one FITS structure = ', sz1, ' bytes'
         splog, 'Number of objects = ', nout
         splog, 'Total size of FITS structure = ', float(sz1)*nout/1.e6, ' Mbyte'
         outdat = replicate(outdat1, nout)
      endif

      ; Assign the tsObj entries first, then over-write any
      ; duplicate tags with ZANS (such as MJD).
      thisdat = replicate(outdat1, 640)
      if (keyword_set(tsobj)) then $
       struct_assign, tsobj, thisdat, /nozero
      struct_assign, zans, thisdat, /nozero

      ; Fill in the first columns of this output structure
      thisdat.progname = plist[ifile].progname
      thisdat.chunkname = plist[ifile].chunkname
      thisdat.platequality = plist[ifile].platequality
      thisdat.platesn2 = plist[ifile].platesn2
      thisdat.smearuse = plist[ifile].smearuse

      ; Get PRIMTARGET+SECTARGET with those values from
      ; the plug-map structure in spPlate file.
      plugmap = mrdfits(fullplatefile[ifile], 5, /silent)
      thisdat.primtarget = plugmap.primtarget
      thisdat.sectarget = plugmap.sectarget

      ; Over-write the MJD with that from the plate file name ???
      ; Early versions of 2D (such as v4_3_1) could have an inconsistent value.
;      thismjd = long( strmid(fileandpath(fullplatefile[ifile]), 13, 5) )
;      thisdat.mjd = thismjd

      ; Copy the data for this plate into the big output structure
      indx = lindgen(640)+640L*ifile
      outdat[indx] = thisdat
   endfor

   splog, 'Time to read data = ', systime(1)-t1, ' sec'

   ;----------
   ; Set the SPECPRIMARY flag to 0 or 1

   t2 = systime(1)

   outdat.specprimary = 1 ; Start as all objects set to primary

   ; Loop through each possible pairing of plates, paying attention
   ; only to those within 4.5 deg of eachother on the sky.
   ; (This is a rather generous match distance; 3.0 deg should be enough
   ; unless there is a mistake somewhere.)

   for ifile1=0, nfile-1 do begin
      for ifile2=ifile1+1, nfile-1 do begin
         adist = djs_diff_angle(plist[ifile1].ra, plist[ifile1].dec, $
          plist[ifile2].ra, plist[ifile2].dec)
         if (adist LT 4.5) then begin
            print, 'Matching plate #', ifile1+1, ' and ', ifile2+1, $
             ' (of ', nfile, ')'
            indx1 = ifile1 * 640L + lindgen(640)
            indx2 = ifile2 * 640L + lindgen(640)
            nn = djs_angle_match(outdat[indx1].plug_ra, outdat[indx1].plug_dec, $
             outdat[indx2].plug_ra, outdat[indx2].plug_dec, dtheta=dtheta, $
             mcount=mcount, mindx=mindx, mmax=1)
            for i1=0, n_elements(indx1)-1 do begin
               if (mcount[i1] GT 1) then $
                message, 'More than 1 match found between two plates!'
               if (mcount[i1] EQ 1) then begin
                  ; Resolve a conflict between object indx1[i1]
                  ; and indx2[mindx[i1]]
                  ; 1) Prefer PROGNAME='main' over any other program names
                  ; 2) Prefer PLATEQUALITY='good' over any other plate quality
                  ; 3) Prefer observations with ZWARNING=0
                  ; 4) Prefer the observation with the larger PLATESN2
                  j1 = indx1[i1]
                  j2 = indx2[mindx[i1]]
                  if ((strmatch(outdat[j1].progname,'main*') EQ 1) $
                   AND (strmatch(outdat[j2].progname,'main*') EQ 0)) then begin
                     outdat[j2].specprimary = 0
                  endif else if ((strmatch(outdat[j1].progname,'main*') EQ 0) $
                   AND (strmatch(outdat[j2].progname,'main*') EQ 1)) then begin
                     outdat[j1].specprimary = 0
                  endif else if ((strmatch(outdat[j1].platequality,'good*') EQ 1) $
                   AND (strmatch(outdat[j2].platequality,'good*') EQ 0)) then begin
                     outdat[j2].specprimary = 0
                  endif else if ((strmatch(outdat[j1].platequality,'good*') EQ 0) $
                   AND (strmatch(outdat[j2].platequality,'good*') EQ 1)) then begin
                     outdat[j1].specprimary = 0
                  endif else if (outdat[j1].zwarning EQ 0 $
                   AND outdat[j2].zwarning NE 0) then begin
                     outdat[j2].specprimary = 0
                  endif else if (outdat[j1].zwarning NE 0 $
                   AND outdat[j2].zwarning EQ 0) then begin
                     outdat[j1].specprimary = 0
                  endif else if (outdat[j1].platesn2 GE outdat[j2].platesn2) then begin
                     outdat[j2].specprimary = 0
                  endif else begin
                     outdat[j1].specprimary = 0
                  endelse
               endif
            endfor
         endif
      endfor
   endfor

   splog, 'Time to assign primaries = ', systime(1)-t2, ' sec'

   ;----------
   ; Write the output FITS file, in chunks of 20 plates

   mwrfits_chunks, outdat, outroot[0]+'.fits', /create, chunksize=640*20

   ;----------
   ; Create the structure for ASCII output

   adat = create_struct( $
    'plate'      ,  0L, $
    'mjd'        ,  0L, $
    'fiberid'    ,  0L, $
    'class'      ,  '', $
    'subclass'   ,  '', $
    'z'          , 0.0, $
    'z_err'      , 0.0, $
    'zwarning'   ,  0L, $
    'rchi2'      , 0.0, $
    'ra'         , 0.0d, $
    'dec'        , 0.0d, $
    'platesn2'   ,  0.0, $
    'modelflux'  ,  fltarr(5), $
    'objc_type'  ,  '', $
    'primtarget' ,  0L, $
    'sectarget'  ,  0L, $
    'progname',     '', $
    'specprimary',  0L, $
    'objtype'    ,  '' )
   sz2 = n_tags(adat, /length)
   splog, 'Size of one ASCII structure = ', sz2, ' bytes'
   splog, 'Number of objects = ', nout
   splog, 'Total size of ASCII structure = ', float(sz2)*nout/1.e6, ' Mbyte'
   adat = replicate(adat, nout)
   struct_assign, outdat, adat

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
outdat = 0 ; Free memory

   struct_print, adat, filename=outroot[0]+'.dat'

   ;----------
   ; Create the merged line data

   plate = adat.plate
   mjd = adat.mjd
   fiberid = adat.fiberid
adat = 0 ; Free memory

   readspec, plate, mjd=mjd, fiberid, zline=linedat
   nobj = n_elements(plate)
   nper = n_elements(linedat) / nobj
   sxaddpar, linehdr, 'DIMS0', nper, ' Number of emission lines'
   sxaddpar, linehdr, 'DIMS1', nobj, ' Number of objects'

   mwrfits, 0, outroot[1]+'.fits', linehdr, /create
   mwrfits_chunks, linedat, outroot[1]+'.fits', chunksize=640*nper*20

   splog, 'Total time = ', systime(1)-t1, ' sec'

   return
end
;------------------------------------------------------------------------------
