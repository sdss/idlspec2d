;+
; NAME:
;   platelist
;
; PURPOSE:
;   Make list of reduced plates
;
; CALLING SEQUENCE:
;   platelist, [fullplatefile, /create, plist= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   infile      - Either a list of combine-plan files or a list of plate files.
;                 If not specified, then search for all plan files matching
;                   '$SPECTRO_DATA/*/spPlancomb-*.par'.
;                 If no such plan files are found, then search for all plate
;                 files matching
;                   '$SPECTRO_DATA/*/spPlate-*.fits'.
;   create      - If set, then re-generate the "platelist.fits" file;
;                 if not set, then simply read this file from a previous call.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   plist       - Output structure with information for each plate.
;
; COMMENTS:
;   Two files are generated: 'platelist.txt' and 'platelist.fits'.
;
;   If INFILE is a list of plan files, i.e.
;     spPlancomb-0306-51690.par
;   then look for the following files for the 2D reductions:
;     spPlancomb-0306-51690.par
;     spDiagcomb-0306-51690.log
;     spPlan2d-0306-51690.par (as specified by 'planfile2d' in spPlancomb)
;     spDiag2d-0306-51690.log
;     spPlate-0306-51690.fits (as specified by 'combinefile' in spPlancomb)
;   and look for the following files for the 1D reductions:
;     spPlan1d-0306-51690.par
;     spZbest-0306-51690.fits
;     spDiag1d-0306-51690.log
;
;   Decide which plates constitute unique tiles with the required S/N,
;   then set QSURVEY=1.  Require (S/N)^2 > 13 for G1,I1,G2,I2.
;   Also require that the target version is not "special" or "devel".
;
; EXAMPLES:
;
; BUGS:
;   Spawns the Unix command 'tail' to get the last line of log files.
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/spPlateList.par
;   $SPECTRO_DATA/platelist.fits
;   $SPECTRO_DATA/platelist.txt
;
; PROCEDURES CALLED:
;   chunkinfo()
;   djs_filepath()
;   fileandpath()
;   headfits()
;   mrdfits()
;   repstr()
;   splog
;   sxpar()
;   tai2airmass()
;   yanny_free
;   yanny_par()
;   yanny_read
;
; REVISION HISTORY:
;   29-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platelist, infile, plist=plist, create=create

   minsn2 = 13.0
   fitsfile = djs_filepath('platelist.fits', root_dir=getenv('SPECTRO_DATA'))
   ascfile = djs_filepath('platelist.txt', root_dir=getenv('SPECTRO_DATA'))

   ;----------
   ; If the /CREATE flag is not set, and the platelist file already exists
   ; on disk, then simply return the info in that file.

   if (NOT keyword_set(create)) then begin
      thisfile = (findfile(fitsfile))[0]
      if (keyword_set(thisfile)) then begin
         plist = mrdfits(thisfile,1)
         return
      endif
   endif

   ;----------
   ; Generate the list of plan files or plate files if not specified

   if (NOT keyword_set(infile)) then $
    infile = djs_filepath('spPlancomb-*.par', $
     root_dir=getenv('SPECTRO_DATA'), subdirectory='*')
   fullfile = findfile(infile, count=nfile)
   if (nfile EQ 0) then begin
      infile = djs_filepath('spPlate-*.fits', $
       root_dir=getenv('SPECTRO_DATA'), subdirectory='*')
      fullfile = findfile(infile, count=nfile)
   endif
   if (nfile EQ 0) then return

   fullfile = fullfile[sort(fullfile)] ; Sort these files

   ;----------
   ; Create output structure

   plist = create_struct( $
    'plate'    , 0L, $
    'tile'     , 0L, $
    'mjd'      , 0L, $
    'ra'       , 0.0, $
    'dec'      , 0.0, $
    'tai'      , 0.0D, $
    'airmass'  , 0.0, $
    'mapname'  , ' ', $
    'vers2d'   , ' ', $
    'verscomb' , ' ', $
    'vers1d'   , ' ', $
    'verstarg' , ' ', $
    'chunkname', ' ', $
    'mjdlist'  , ' ', $
    'sn2_g1'   , 0.0, $
    'sn2_r1'   , 0.0, $
    'sn2_i1'   , 0.0, $
    'sn2_g2'   , 0.0, $
    'sn2_r2'   , 0.0, $
    'sn2_i2'   , 0.0, $
    'n_galaxy' , 0L,  $
    'n_qso'    , 0L,  $
    'n_star'   , 0L,  $
    'n_unknown', 0L,  $
    'n_sky'    , 0L, $
    'status2d' , 'Missing', $
    'status1d' , 'Missing', $
    'qsurvey'  , 0L, $
    'public'  , '' )
   plist = replicate(plist, nfile)

   ;----------
   ; Read the data file with the public plate information

   publicfile = filepath('spPlateList.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   yanny_read, publicfile, pdata
   publicdata = *pdata[0]
   yanny_free, pdata

   ;---------------------------------------------------------------------------
   ; Loop through all files
   ;---------------------------------------------------------------------------

   for ifile=0, nfile-1 do begin

      splog, 'Looking at ' + fullfile[ifile]

      junk = fileandpath(fullfile[ifile], path=path)  ; Determine PATH

      ;----------
      ; Test if INFILE specifies Yanny param files for spPlancomb.

      if (strmid(fullfile[ifile],strlen(fullfile[ifile])-4) EQ '.par') $
       then begin
         combparfile = fullfile
         yanny_read, fullfile[ifile], hdr=hdrp
         platefile = $
          djs_filepath(yanny_par(hdrp, 'combinefile'), root_dir=path)
      endif else begin
         platefile = fullfile[ifile]
         combparfile = repstr(platefile, 'spPlate', 'spPlancomb')
         combparfile = repstr(combparfile, '.fits', '.par')
      endelse

      ;----------
      ; Determine names of associated files

      comblogfile = repstr(combparfile, '.par', '.log')
      platemjd = strmid(fileandpath(platefile), 8, 10)
      zbestfile = 'spZbest-' + platemjd + '.fits'

      fullzfile = djs_filepath(zbestfile, root_dir=path)

      ;----------
      ; Read plate file - get status of 2D

      hdr1 = headfits(platefile)
      if (size(hdr1, /tname) EQ 'STRING') then begin
;         plist[ifile].plate = sxpar(hdr1, 'PLATEID')
;         plist[ifile].mjd = sxpar(hdr1, 'MJD')
         plist[ifile].mjdlist = sxpar(hdr1, 'MJDLIST')
         plist[ifile].tile = sxpar(hdr1, 'TILEID')
         plist[ifile].ra = sxpar(hdr1, 'RA')
         plist[ifile].dec = sxpar(hdr1, 'DEC')
         plist[ifile].tai = sxpar(hdr1, 'TAI')
         plist[ifile].airmass = tai2airmass(plist[ifile].ra, $
          plist[ifile].dec, tai=plist[ifile].tai)
         plist[ifile].sn2_g1 = sxpar(hdr1, 'SPEC1_G')
         plist[ifile].sn2_r1 = sxpar(hdr1, 'SPEC1_R')
         plist[ifile].sn2_i1 = sxpar(hdr1, 'SPEC1_I')
         plist[ifile].sn2_g2 = sxpar(hdr1, 'SPEC2_G')
         plist[ifile].sn2_r2 = sxpar(hdr1, 'SPEC2_R')
         plist[ifile].sn2_i2 = sxpar(hdr1, 'SPEC2_I')
         plist[ifile].mapname = strtrim(sxpar(hdr1, 'NAME'))
         plist[ifile].vers2d = strtrim(sxpar(hdr1, 'VERS2D'))
         plist[ifile].verscomb = strtrim(sxpar(hdr1, 'VERSCOMB'))
         plist[ifile].status2d = 'Done'
      endif else begin
         ; Get the mapping name from the combine par file
         yanny_read, combparfile[ifile], pp, hdr=hdrcomb
         plist[ifile].mapname = (*pp[0])[0].mapname
         yanny_free, pp

         ;----------
         ; Find the state of the 2D reductions -- spPlate file is missing

         statusdone = 0
         statusrun = 0
         statusmissing = 0

         ; Check status of individual 2D runs
         planlist = yanny_par(hdrcomb, 'planfile2d') ; Assume we find this
         for iplan=0, n_elements(planlist)-1 do begin
            yanny_read, djs_filepath(planlist[iplan], root_dir=path), hdr=hdr2d
            plist[ifile].mjdlist = strtrim(plist[ifile].mjdlist $
             + ' ' + yanny_par(hdr2d, 'MJD'),2)
            thislogfile = djs_filepath(yanny_par(hdr2d, 'logfile'), root_dir=path)
            thislogfile = (findfile(thislogfile))[0]
            if (keyword_set(thislogfile)) then begin
               spawn, 'tail -1 '+thislogfile, lastline
               if (strmatch(lastline[0], '*Successful completion*')) then begin
                  ; Case where this 2D log file completed
                  statusdone = statusdone + 1
               endif else begin
                  ; Case where this 2D log file isn't completed
                  statusrun = statusrun + 1
               endelse
            endif else begin
               ; Case where this 2D log file missing
               statusmissing = statusmissing + 1
            endelse
         endfor

         ; Check status of 2D combining
         thislogfile = djs_filepath(yanny_par(hdrcomb, 'logfile'), root_dir=path)
         thislogfile = (findfile(thislogfile))[0]
         if (keyword_set(thislogfile)) then begin
            spawn, 'tail -1 '+thislogfile, lastline
            if (strmatch(lastline[0], '*Successful completion*')) then begin
               ; Case where this combine log file completed
               statusdone = statusdone + 1
            endif else begin
               ; Case where this combine log file isn't completed
               statusrun = statusrun + 1
            endelse
         endif else begin
            ; Case where this combine log file missing
            statusmissing = statusmissing + 1
         endelse

         if (statusmissing EQ 0 AND statusrun EQ 0) then $
          plist[ifile].status2d = 'FAILED' $ ; Should have found spPlate file
         else if (statusrun EQ 0 AND statusdone EQ 0) then $
          plist[ifile].status2d = 'Pending' $ ; No log files created
         else $
          plist[ifile].status2d = 'RUNNING' ; Some log files created

      endelse

      ;----------
      ; Get the following from the file names, since sometimes they
      ; are wrong in the file headers!!

      plist[ifile].plate = long( strmid(fileandpath(platefile), 8, 4) )
      plist[ifile].mjd = long( strmid(fileandpath(platefile), 13, 5) )

      ;----------
      ; Determine the chunk name and the version of target used

      cinfo = chunkinfo(plist[ifile].plate)
      plist[ifile].chunkname = cinfo.chunkname
      plist[ifile].verstarg = cinfo.verstarg

      ;----------
      ; Determine which public data release has this plate+MJD

      j = where(plist[ifile].plate EQ publicdata.plate $
       AND plist[ifile].mjd EQ publicdata.mjd)
      if (j[0] NE -1) then $
       plist[ifile].public = publicdata[j[0]].public

      ;----------
      ; The RA,DEC in the header is sometimes wrong, so try to derive
      ; the field center from the plug-map information.  Choose the
      ; coordinates of the object closest to the center of the plate
      ; as defined by XFOCAL=YFOCAL=0.

      plug = mrdfits(platefile, 5, /silent)
      if (keyword_set(plug)) then begin
         iobj = where(strtrim(plug.holetype,2) EQ 'OBJECT')
         junk = min( plug[iobj].xfocal^2 + plug[iobj].yfocal^2, imin)
         plist[ifile].ra = plug[iobj[imin]].ra
         plist[ifile].dec = plug[iobj[imin]].dec
         plist[ifile].airmass = tai2airmass(plist[ifile].ra, $
          plist[ifile].dec, tai=plist[ifile].tai)
      endif

      ;----------
      ; Read Zbest file - get status of 1D

      hdr2 = headfits(fullzfile)
      if (size(hdr2, /tname) EQ 'STRING') then begin
         zans = mrdfits(fullzfile, 1, /silent)
         class = strtrim(zans.class,2)
         ; Use the ZWARNING flag if it exists to identify SKY or UNKNOWN.
         if ((where(tag_names(zans) EQ 'ZWARNING'))[0] NE -1) then $
          zwarning = zans.zwarning $
         else $
          zwarning = bytarr(n_elements(zans))
         qsky = (zwarning AND 1) NE 0
         plist[ifile].n_galaxy = total(class EQ 'GALAXY' AND zwarning EQ 0)
         plist[ifile].n_qso = total(class EQ 'QSO' AND zwarning EQ 0)
         plist[ifile].n_star = total(class EQ 'STAR' AND zwarning EQ 0)
         plist[ifile].n_unknown = total(class EQ 'UNKNOWN' $
          OR (zwarning NE 0 AND qsky EQ 0))
         plist[ifile].n_sky = total(class EQ 'SKY' OR qsky EQ 1)
         plist[ifile].vers1d = strtrim(sxpar(hdr2, 'VERS1D'))
         plist[ifile].status1d = 'Done'
      endif else begin
         ;----------
         ; Find the state of the 1D reductions -- spZbest file is missing

         thislogfile = repstr(fileandpath(platefile), 'spPlate', 'spDiag1d')
         thislogfile = repstr(thislogfile, '.fits', '.log')
         thislogfile = djs_filepath(thislogfile, root_dir=path)
         thislogfile = (findfile(thislogfile))[0]
         if (keyword_set(thislogfile)) then begin
            spawn, 'tail -1 '+thislogfile, lastline
            if (strmatch(lastline[0], '*Successful completion*')) then begin
               ; Case where this 1D log file completed
               plist[ifile].status1d = 'FAILED'; Should have found spZbest file
            endif else begin
               ; Case where this 1D log file isn't completed
               plist[ifile].status1d = 'RUNNING'
            endelse
         endif else begin
            ; Case where this 1D log file missing
            plist[ifile].status1d = 'Pending'
         endelse
      endelse

   endfor

   ;----------
   ; Remove from the plate list earlier reductions of the same plugging
   ; (keeping only the most recent MJD of each plugging).

   qkeep = bytarr(nfile)

   ; First get the unique list of MAPNAME, then mark the most recent MJD
   ; of each as the good one.
   isort = sort(plist.mapname)
   isort = isort[ uniq(plist[isort].mapname) ]
   maplist = plist[isort].mapname
   for imap=0, n_elements(maplist)-1 do begin
      indx = where(plist.mapname EQ maplist[imap])
      junk = max(plist[indx].mjd, imax)
      qkeep[indx[imax]] = 1
   endfor

   ; Don't discard any whre MAPNAME isn't set
   indx = where(strtrim(plist.mapname) EQ '')
   if (indx[0] NE -1) then qkeep[indx] = 1

   ; List partially-combined plates that we're discarding from the list
   for ifile=0, nfile-1 do $
    if (qkeep[ifile] NE 1) then $
     splog, 'Discard partially-combined ' + fullfile[ifile]

   ; Trim the plate list, and update NFILE to this trimmed number
   plist = plist[where(qkeep, nfile)]

   ;----------
   ; Make a list of one S/N for each plate which is the minimum of
   ; G1, I1, G2, I2.

   snvec = fltarr(nfile)
   for ifile=0, nfile-1 do $
    snvec[ifile] = min([plist[ifile].sn2_g1, plist[ifile].sn2_i1, $
     plist[ifile].sn2_g2, plist[ifile].sn2_i2])

   ;----------
   ; Decide which plates constitute unique tiles with the required S/N,
   ; then set QSURVEY=1.
   ; Also insist that the target version isn't "special" or "devel".

   ; First get the unique list of TILE
   isort = sort(plist.tile)
   isort = isort[ uniq(plist[isort].tile) ]
   tilelist = plist[isort].tile

   for itile=0, n_elements(tilelist)-1 do begin
      indx = where(plist.tile EQ tilelist[itile] $
       AND strtrim(plist.verstarg,2) NE 'special' $
       AND strtrim(plist.verstarg,2) NE 'devel')
      if (indx[0] NE -1) then begin
         snbest = max(snvec[indx], ibest)
         if (snbest GE minsn2) then plist[indx[ibest]].qsurvey = 1
      endif
   endfor

   ;---------------------------------------------------------------------------
   ; Open output file

   if (keyword_set(ascfile)) then $
    openw, olun, ascfile, /get_lun $
   else $
    olun = -1L

   printf, olun, 'PLATE  MJD   RA    DEC   SN2_G1 SN2_I1 SN2_G2 SN2_I2 ' $
    + 'Ngal Nqso Nsta Nunk Nsky Stat2D  Stat1D  Vers2D    Vers1D    ? Pub'
   printf, olun, '-----  ----- ----- ----- ------ ------ ------ ------ ' $
    + '---- ---- ---- ---- ---- ------- ------- --------- --------- - ---'

   ;----------
   ; Loop through all files

   for ifile=0, nfile-1 do begin
      printf, olun, plist[ifile].plate, plist[ifile].mjd, $
       plist[ifile].ra, plist[ifile].dec, plist[ifile].sn2_g1, $
       plist[ifile].sn2_i1, plist[ifile].sn2_g2, plist[ifile].sn2_i2, $
       plist[ifile].n_galaxy, plist[ifile].n_qso, plist[ifile].n_star, $
       plist[ifile].n_unknown, plist[ifile].n_sky, $
       plist[ifile].status2d, plist[ifile].status1d, $
       plist[ifile].vers2d, plist[ifile].vers1d, $
       plist[ifile].qsurvey, plist[ifile].public, $
       format='(i5,i7,f6.1,f6.1,4f7.1,5i5,2(1x,a7),2(1x,a9),i2,a4)'
   endfor

   ;----------
   ; Close output file

   if (keyword_set(ascfile)) then begin
      close, olun
      free_lun, olun
   endif

   ;----------
   ; Write the FITS binary table

   if (keyword_set(fitsfile)) then $
    mwrfits, plist, fitsfile, /create

end
;------------------------------------------------------------------------------
