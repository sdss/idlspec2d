;+
; NAME:
;   platelist
;
; PURPOSE:
;   Make list of reduced plates
;
; CALLING SEQUENCE:
;   platelist, [infile, /create, /purge2d, /purge1d, /killpartial, plist= ]
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
;   purge2d     - If set, then delete all log files for plates that are
;                 considered to be 'RUNNING', but not those that are 'Done',
;                 'Pending' or 'FAILED'.  Those plates are then listed as
;                 'Pending'.  Setting /PURGE2D also sets /CREATE.
;                 Deleting these log files will cause the next invocation
;                 of BATCH2D to re-reduce those plates.
;   purge1d     - If set, then delete all log files for plates that are
;                 considered to be 'RUNNING', but not those that are 'Done',
;                 'Pending' or 'FAILED'.  Those plates are then listed as
;                 'Pending'.  Setting /PURGE1D also sets /CREATE.
;                 Deleting these log files will cause the next invocation
;                 of BATCH1D to re-reduce those plates.
;   killpartial - If set, then delete all files associated with a combine
;                 of only some nights of a multi-night plate.  Such files
;                 can be produced by the Spectro-Robot when it fully reduces
;                 data from one night, but then more data is obtained for
;                 that plugging of the same plate on a later date.  This
;                 deletes spPlate and spZ files and their logs files.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   plist       - Output structure with information for each plate.
;
; COMMENTS:
;   Three files are generated:
;     $SPECTRO_DATA/platelist.fits
;     $SPECTRO_DATA/platelist.txt
;     $SPECTRO_DATA/platelist.html
;     $SPECTRO_DATA/platequal.txt
;     $SPECTRO_DATA/platequal.html
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
;   PLATESN2 is set to the minimum of the 4 cameras.
;   PLATEQUALITY defaults to 'good'.
;   PLATEQUALITY is set to 'marginal' if MINSN2 < 15.0
;                          'bad'      if MINSN2 < 13.0
;   PLATEQUALITY is set to 'marginal' if FBADPIX > 0.05
;                          'bad'      if FBADPIX > 0.10
;
;   Decide which plates constitute unique tiles with the required S/N,
;   then set QSURVEY=1.  Require PLATEQUALITY='good' or 'marginal'.
;   Also require PROGNAME='main'.
;
; EXAMPLES:
;
; BUGS:
;   Spawns the Unix command 'tail' to get the last line of log files.
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/spChunkList.par
;   $IDLSPEC2D_DIR/etc/spPlateList.par
;
; PROCEDURES CALLED:
;   apo_checklimits()
;   chunkinfo()
;   copy_struct_inx
;   djs_filepath()
;   fileandpath()
;   headfits()
;   mrdfits()
;   repstr()
;   rmfile
;   splog
;   sxpar()
;   tai2airmass()
;   yanny_free
;   yanny_par()
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;   platelist_write
;
; REVISION HISTORY:
;   29-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platelist_write, plist, trimtags=trimtags, alias=alias, $
 fileprefix=fileprefix, title=title

   ascfile = djs_filepath(fileprefix+'.txt', root_dir=getenv('SPECTRO_DATA'))
   htmlfile = djs_filepath(fileprefix+'.html', root_dir=getenv('SPECTRO_DATA'))

   trimdat  = struct_trimtags(plist, select_tags=trimtags[0,*])
   trimstring = struct_trimtags(plist, select_tags=trimtags[0,*], $
    format=trimtags[1,*])
   struct_print, trimstring, filename=ascfile, fdigit=3, alias=alias

   for itag=0, n_tags(trimdat)-1 do begin
      for iarr=0, n_elements(trimdat[0].(itag))-1 do begin
         for irow=0, n_elements(trimdat)-1 do begin
            markstring = apo_checklimits('SUMMARY', $
             strupcase(trimtags[0,itag]), '', $
             trimdat[irow].(itag)[iarr], /html)
            trimstring[irow].(itag)[iarr] = markstring $
             + trimstring[irow].(itag)[iarr]
         endfor
      endfor
   endfor
   struct_print, trimstring, /html, alias=alias, tarray=tarray
   openw, lun, htmlfile, /get_lun
   printf, lun, '<HTML>'
   printf, lun, '<HEAD>'
   printf, lun, '<TITLE>' + title + '</TITLE>'
   printf, lun, '<H1>' + title + '</H1>'
   printf, lun, '<BODY>'
   printf, lun, '<A HREF="' + fileandpath(ascfile) + '">ASCII version</A>'
   printf, lun, '<A HREF="platelist.fits">FITS version</A>'
   for iline=0, n_elements(tarray)-1 do $
    printf, lun, tarray[iline]
   printf, lun, '</BODY>'
   close, lun

   return
end
;------------------------------------------------------------------------------
pro platelist, infile, plist=plist, create=create, $
 purge2d=purge2d, purge1d=purge1d, killpartial=killpartial

   fitsfile = djs_filepath('platelist.fits', root_dir=getenv('SPECTRO_DATA'))

   ;----------
   ; If the /CREATE flag is not set, and the platelist file already exists
   ; on disk, then simply return the info in that file.

   if (NOT keyword_set(create) AND NOT keyword_set(purge2d) $
    AND NOT keyword_set(purge1d)) then begin
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
    'cartid'   , 0L, $
    'tai'      , 0.0D, $
    'tai_beg'  , 0.0D, $
    'tai_end'  , 0.0D, $
    'airmass'  , 0.0, $
    'exptime'  , 0.0, $
    'mapname'  , ' ', $
    'vers2d'   , ' ', $
    'verscomb' , ' ', $
    'vers1d'   , ' ', $
    'progname' , ' ', $
    'chunkname', ' ', $
    'platequality' , ' ', $
    'platesn2' , 0.0, $
    'qsurvey'  , 0L,  $
    'mjdlist'  , ' ', $
    'nexp'     , 0L,  $
    'expt_b1'  , 0.0, $
    'expt_b2'  , 0.0, $
    'expt_r1'  , 0.0, $
    'expt_r2'  , 0.0, $
    'sn2_g1'   , 0.0, $
    'sn2_r1'   , 0.0, $
    'sn2_i1'   , 0.0, $
    'sn2_g2'   , 0.0, $
    'sn2_r2'   , 0.0, $
    'sn2_i2'   , 0.0, $
    'smearuse' , ' ', $
    'nguide'   , 0L , $
    'seeing20' , 0.0, $
    'seeing50' , 0.0, $
    'seeing80' , 0.0, $
    'rmsoff20' , 0.0, $
    'rmsoff50' , 0.0, $
    'rmsoff80' , 0.0, $
    'airtemp'  , 0.0, $
    'xsigma'   , 0.0, $
    'xsigmin'  , 0.0, $
    'xsigmax'  , 0.0, $
    'wsigma'   , 0.0, $
    'wsigmin'  , 0.0, $
    'wsigmax'  , 0.0, $
    'xchi2'    , 0.0, $
    'xchi2min' , 0.0, $
    'xchi2max' , 0.0, $
    'skychi2'  , 0.0, $
    'schi2min' , 0.0, $
    'schi2max' , 0.0, $
    'fbadpix'  , 0.0, $
    'fbadpix1' , 0.0, $
    'fbadpix2' , 0.0, $
    'n_galaxy' , 0L,  $
    'n_qso'    , 0L,  $
    'n_star'   , 0L,  $
    'n_unknown', 0L,  $
    'n_sky'    , 0L, $
    'n_target_main',  0L, $
    'n_target_lrg' ,  0L, $
    'n_target_qso' ,  0L, $
    'success_main' , 0.0, $
    'success_lrg'  , 0.0, $
    'success_qso'  , 0.0, $
    'status2d' , 'Missing', $
    'statuscombine', 'Missing', $
    'status1d' , 'Missing', $
    'public'  , ' ', $
    'qualcomments'  , ' ' )
   trimtags1 = [ $
    ['plate'        ,   'i4'], $
    ['mjd'          ,   'i5'], $
    ['ra'           , 'f6.2'], $
    ['dec'          , 'f6.2'], $
    ['progname'     ,    'a'], $
    ['chunkname'    ,    'a'], $
    ['platequality' ,    'a'], $
    ['platesn2'     , 'f5.1'], $
    ['n_galaxy'     ,   'i3'], $
    ['n_qso'        ,   'i3'], $
    ['n_star'       ,   'i3'], $
    ['n_unknown'    ,   'i3'], $
    ['n_sky'        ,   'i3'], $
    ['public'       ,    'a']  ]
   trimtags2 = [ $
    ['plate'        ,   'i4'], $
    ['mjd'          ,   'i5'], $
    ['progname'     ,    'a'], $
    ['chunkname'    ,    'a'], $
    ['sn2_g1'       , 'f5.1'], $
    ['sn2_i1'       , 'f5.1'], $
    ['sn2_g2'       , 'f5.1'], $
    ['sn2_i2'       , 'f5.1'], $
    ['fbadpix'      , 'f5.3'], $
    ['success_main' , 'f5.1'], $
    ['success_lrg'  , 'f5.1'], $
    ['success_qso'  , 'f5.1'], $
    ['status2d'     ,    'a'], $
    ['statuscombine',    'a'], $
    ['status1d'     ,    'a'], $
    ['platequality' ,    'a'], $
    ['qualcomments' ,    'a']  ]
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

   platefile = strarr(nfile)
   combparfile = strarr(nfile)
   comblogfile = strarr(nfile)
   combpsfile = strarr(nfile)
   zlogfile = strarr(nfile)
   zbestfile = strarr(nfile)
   zallfile = strarr(nfile)

   for ifile=0, nfile-1 do begin

      splog, 'Looking at ' + fullfile[ifile]

      junk = fileandpath(fullfile[ifile], path=path)  ; Determine PATH

      ;----------
      ; Test if INFILE specifies Yanny param files for spPlancomb.

      if (strmid(fullfile[ifile],strlen(fullfile[ifile])-4) EQ '.par') $
       then begin
         combparfile[ifile] = fullfile[ifile]
         yanny_read, fullfile[ifile], hdr=hdrp
         platefile[ifile] = $
          djs_filepath(yanny_par(hdrp, 'combinefile'), root_dir=path)
      endif else begin
         platefile[ifile] = fullfile[ifile]
         combparfile[ifile] = repstr(platefile[ifile], 'spPlate', 'spPlancomb')
         combparfile[ifile] = repstr(combparfile[ifile], '.fits', '.par')
      endelse

      ;----------
      ; Determine names of associated files

      comblogfile[ifile] = repstr(combparfile[ifile], '.par', '.log')
      comblogfile[ifile] = repstr(comblogfile[ifile], 'spPlancomb', 'spDiagcomb')
      combpsfile[ifile] = repstr(comblogfile[ifile], '.log', '.ps')
      platemjd = strmid(fileandpath(platefile[ifile]), 8, 10)
      zbestfile[ifile] = djs_filepath('spZbest-' + platemjd + '.fits', $
       root_dir=path)
      zallfile[ifile] = djs_filepath('spZall-' + platemjd + '.fits', $
       root_dir=path)

      ; Read the combine plan file to get the list of all the 2D plan files
      ; from its Yanny header.
      ; Also get the mapping name from the combine par file in case we were
      ; unable to get it from the spPlate file.
      yanny_read, combparfile[ifile], pp, hdr=hdrcomb
      plist[ifile].mapname = (*pp[0])[0].mapname
      yanny_free, pp

      ;----------
      ; Find the state of the 2D reductions (not the combine step)

      statusdone = 0
      statusrun = 0
      statusmissing = 0

      ; Check status of individual 2D runs
      planlist = yanny_par(hdrcomb, 'planfile2d') ; Assume we find this
      planlist = djs_filepath(planlist, root_dir=path)
      logfile2d = '' ; List of 2D log files that exist
      for iplan=0, n_elements(planlist)-1 do begin
         yanny_read, planlist[iplan], hdr=hdr2d
         plist[ifile].mjdlist = strtrim(plist[ifile].mjdlist $
          + ' ' + yanny_par(hdr2d, 'MJD'),2)
         thislogfile = djs_filepath(yanny_par(hdr2d, 'logfile'), root_dir=path)
         thislogfile = (findfile(thislogfile))[0]
         if (keyword_set(thislogfile)) then begin
            if (NOT keyword_set(logfile2d)) then logfile2d = thislogfile $
             else logfile2d = [logfile2d, thislogfile]
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

      if (statusmissing GT 0 AND statusrun EQ 0) then begin
         plist[ifile].status2d = 'Pending'
      endif else if (statusmissing GT 0 OR statusrun GT 0) then begin
         plist[ifile].status2d = 'RUNNING'
      endif else begin
         plist[ifile].status2d = 'Done'
      endelse

      ;----------
      ; Read plate file - get status of Combine

      hdr1 = headfits(platefile[ifile])
      if (size(hdr1, /tname) EQ 'STRING') then begin
;         plist[ifile].plate = sxpar(hdr1, 'PLATEID')
;         plist[ifile].mjd = sxpar(hdr1, 'MJD')
         plist[ifile].mjdlist = sxpar(hdr1, 'MJDLIST')
         plist[ifile].tile = sxpar(hdr1, 'TILEID')
         plist[ifile].ra = sxpar(hdr1, 'RA')
         plist[ifile].dec = sxpar(hdr1, 'DEC')
         plist[ifile].cartid = sxpar(hdr1, 'CARTID')
         plist[ifile].tai = sxpar(hdr1, 'TAI')
         plist[ifile].tai_beg = sxpar(hdr1, 'TAI-BEG')
         plist[ifile].tai_end = sxpar(hdr1, 'TAI-END')
         plist[ifile].airmass = sxpar(hdr1, 'AIRMASS')
;         plist[ifile].airmass = tai2airmass(plist[ifile].ra, $
;          plist[ifile].dec, tai=plist[ifile].tai)
         plist[ifile].exptime = sxpar(hdr1, 'EXPTIME')
         plist[ifile].nexp = sxpar(hdr1, 'NEXP')
         plist[ifile].expt_b1 = sxpar(hdr1, 'EXPT_B1')
         plist[ifile].expt_b2 = sxpar(hdr1, 'EXPT_B2')
         plist[ifile].expt_r1 = sxpar(hdr1, 'EXPT_R1')
         plist[ifile].expt_r2 = sxpar(hdr1, 'EXPT_R2')
         plist[ifile].sn2_g1 = sxpar(hdr1, 'SPEC1_G')
         plist[ifile].sn2_r1 = sxpar(hdr1, 'SPEC1_R')
         plist[ifile].sn2_i1 = sxpar(hdr1, 'SPEC1_I')
         plist[ifile].sn2_g2 = sxpar(hdr1, 'SPEC2_G')
         plist[ifile].sn2_r2 = sxpar(hdr1, 'SPEC2_R')
         plist[ifile].sn2_i2 = sxpar(hdr1, 'SPEC2_I')
         plist[ifile].smearuse = sxpar(hdr1, 'SMEARUSE') ? 'T' : 'F'
         plist[ifile].nguide = sxpar(hdr1, 'NGUIDE')
         plist[ifile].seeing20 = sxpar(hdr1, 'SEEING20')
         plist[ifile].seeing50 = sxpar(hdr1, 'SEEING50')
         plist[ifile].seeing80 = sxpar(hdr1, 'SEEING80')
         plist[ifile].rmsoff20 = sxpar(hdr1, 'RMSOFF20')
         plist[ifile].rmsoff50 = sxpar(hdr1, 'RMSOFF50')
         plist[ifile].rmsoff80 = sxpar(hdr1, 'RMSOFF80')
         plist[ifile].airtemp = sxpar(hdr1, 'AIRTEMP')
         plist[ifile].xsigma = sxpar(hdr1, 'XSIGMA')
         plist[ifile].xsigmin = sxpar(hdr1, 'XSIGMIN')
         plist[ifile].xsigmax = sxpar(hdr1, 'XSIGMAX')
         plist[ifile].wsigma = sxpar(hdr1, 'WSIGMA')
         plist[ifile].wsigmin = sxpar(hdr1, 'WSIGMIN')
         plist[ifile].wsigmax = sxpar(hdr1, 'WSIGMAX')
         plist[ifile].xchi2 = sxpar(hdr1, 'XCHI2')
         plist[ifile].xchi2min = sxpar(hdr1, 'XCHI2MIN')
         plist[ifile].xchi2max = sxpar(hdr1, 'XCHI2MAX')
         plist[ifile].skychi2 = sxpar(hdr1, 'SKYCHI2')
         plist[ifile].schi2min = sxpar(hdr1, 'SCHI2MIN')
         plist[ifile].schi2max = sxpar(hdr1, 'SCHI2MAX')
         plist[ifile].fbadpix = sxpar(hdr1, 'FBADPIX')
         plist[ifile].fbadpix1 = sxpar(hdr1, 'FBADPIX1')
         plist[ifile].fbadpix2 = sxpar(hdr1, 'FBADPIX2')
         plist[ifile].mapname = strtrim(sxpar(hdr1, 'NAME'))
         plist[ifile].vers2d = strtrim(sxpar(hdr1, 'VERS2D'))
         plist[ifile].verscomb = strtrim(sxpar(hdr1, 'VERSCOMB'))
         plist[ifile].statuscombine = 'Done'
      endif else begin
         ; Case where no spPlate file exists
         thislogfile = djs_filepath(yanny_par(hdrcomb, 'logfile'), root_dir=path)
         thislogfile = (findfile(thislogfile))[0]
         if (keyword_set(thislogfile)) then begin
            spawn, 'tail -1 '+thislogfile, lastline
            if (strmatch(lastline[0], '*Successful completion*')) then begin
               ; Case where this combine log file completed
               plist[ifile].statuscombine = 'FAILED'
            endif else begin
               ; Case where this combine log file isn't completed
               spawn, 'grep ABORT '+thislogfile, abortline
               if (keyword_set(abortline)) then begin
                  plist[ifile].statuscombine = 'FAILED' ; Combining step aborted
               endif else begin
                  plist[ifile].statuscombine = 'RUNNING'
               endelse
            endelse
         endif else begin
            ; Case where this combine log file missing
            plist[ifile].statuscombine = 'Pending'
         endelse

         if (keyword_set(purge2d) $
          AND plist[ifile].statuscombine NE 'Done') then begin
            splog, 'PURGE2D ', logfile2d
            rmfile, logfile2d
            splog, 'PURGE2D ', comblogfile[ifile]
            rmfile, comblogfile[ifile]
            splog, 'PURGE2D ', combpsfile[ifile]
            rmfile, combpsfile[ifile]
            plist[ifile].status2d = 'Pending'
            plist[ifile].statuscombine = 'Pending'
         endif
      endelse

      ;----------
      ; Get the following from the file names, since sometimes they
      ; are wrong in the file headers!!

      plist[ifile].plate = long( strmid(fileandpath(platefile[ifile]), 8, 4) )
      plist[ifile].mjd = long( strmid(fileandpath(platefile[ifile]), 13, 5) )

      ;----------
      ; Determine the chunk name and the version of target used

      cinfo = chunkinfo(plist[ifile].plate)
      plist[ifile].chunkname = cinfo.chunkname
      plist[ifile].progname = cinfo.progname

      ;----------
      ; Determine which public data release has this plate+MJD

      j = where(plist[ifile].plate EQ publicdata.plate $
       AND plist[ifile].mjd EQ publicdata.mjd)
      if (j[0] NE -1) then begin
         copy_struct_inx, publicdata[j[0]], plist, index_to=ifile
      endif

      ;----------
      ; The RA,DEC in the header is sometimes wrong, so try to derive
      ; the field center from the plug-map information.  Choose the
      ; coordinates of the object closest to the center of the plate
      ; as defined by XFOCAL=YFOCAL=0.

;      plug = mrdfits(platefile[ifile], 5, /silent)
;      if (keyword_set(plug)) then begin
;         iobj = where(strtrim(plug.holetype,2) EQ 'OBJECT')
;         junk = min( plug[iobj].xfocal^2 + plug[iobj].yfocal^2, imin)
;         plist[ifile].ra = plug[iobj[imin]].ra
;         plist[ifile].dec = plug[iobj[imin]].dec
;         plist[ifile].airmass = tai2airmass(plist[ifile].ra, $
;          plist[ifile].dec, tai=plist[ifile].tai)
;      endif

      ;----------
      ; Read Zbest file - get status of 1D

      hdr2 = headfits(zbestfile[ifile])
      if (size(hdr2, /tname) EQ 'STRING') then begin
         zans = mrdfits(zbestfile[ifile], 1, /silent)
         plug = mrdfits(platefile[ifile], 5, /silent)
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

         nobj = n_elements(zans)
         targets = strarr(nobj)
         for iobj=0, nobj-1 do $
          targets[iobj] = sdss_flagname('TARGET',plug[iobj].primtarget, /silent, /concat)+' '
         imain = where(strmatch(targets,'*GALAXY *') $
          OR strmatch(targets,'*GALAXY_BIG *') $
          OR strmatch(targets,'*GALAXY_BRIGHT_CORE *'), nmain)
         ilrg = where(strmatch(targets,'*GALAXY_RED *') $
          OR strmatch(targets,'*GALAXY_RED_II *'), nlrg)
         iqso = where(strmatch(targets,'*QSO_HIZ *') $
          OR strmatch(targets,'*QSO_CAP *') $
          OR strmatch(targets,'*QSO_SKIRT *'), nqso)
          OR strmatch(targets,'*QSO_FIRST_CAP *') $
          OR strmatch(targets,'*QSO_FIRST_SKIRT *'), nqso)
         plist[ifile].n_target_main = nmain
         plist[ifile].n_target_lrg = nlrg
         plist[ifile].n_target_qso = nqso
         if (nmain GT 0) then $
          plist[ifile].success_main = $
           100 * total(zans[imain].zwarning EQ 0 $
           AND (strmatch(zans[imain].class,'GALAXY*') $
            OR strmatch(zans[imain].class,'QSO*'))) / nmain
         if (nlrg GT 0) then $
          plist[ifile].success_lrg = $
           100 * total(zans[ilrg].zwarning EQ 0 $
           AND strmatch(zans[ilrg].class,'GALAXY*')) / nlrg
         if (nqso GT 0) then $
          plist[ifile].success_qso = $
           100 * total(zans[iqso].zwarning EQ 0 $
           AND strmatch(zans[iqso].class,'QSO*')) / nqso
      endif else begin
         ;----------
         ; Find the state of the 1D reductions -- spZbest file is missing

         zlogfile[ifile] = repstr(fileandpath(platefile[ifile]), 'spPlate', 'spDiag1d')
         zlogfile[ifile] = repstr(zlogfile[ifile], '.fits', '.log')
         zlogfile[ifile] = djs_filepath(zlogfile[ifile], root_dir=path)
         zlogfile[ifile] = (findfile(zlogfile[ifile]))[0]
         if (keyword_set(zlogfile[ifile])) then begin
            spawn, 'tail -1 '+zlogfile[ifile], lastline
            if (strmatch(lastline[0], '*Successful completion*')) then begin
               ; Case where this 1D log file completed
               plist[ifile].status1d = 'FAILED'; Should have found spZbest file
            endif else begin
               ; Case where this 1D log file isn't completed
               if (keyword_set(purge1d)) then begin
                  splog, 'PURGE1D ', zlogfile[ifile]
                  rmfile, zlogfile[ifile]
                  plist[ifile].status1d = 'Pending'
               endif else begin
                  plist[ifile].status1d = 'RUNNING'
               endelse
            endelse
         endif else begin
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

   ; Don't discard any where MAPNAME isn't set
   indx = where(strtrim(plist.mapname) EQ '')
   if (indx[0] NE -1) then qkeep[indx] = 1

   ; List partially-combined plates that we're discarding from the list
   for ifile=0, nfile-1 do begin
      if (qkeep[ifile] NE 1) then begin
         splog, 'Discard partially-combined ' + combparfile[ifile]
         if (keyword_set(killpartial)) then begin
; Also should be purging the spFluxcalib files and spSN2d files !!!???
            killfiles = [ combparfile[ifile], $
                          comblogfile[ifile], $
                          combpsfile[ifile], $
                          platefile[ifile], $
                          zbestfile[ifile], $
                          zallfile[ifile], $
                          zlogfile[ifile] ]
            for ikill=0, n_elements(killfiles)-1 do $
             splog, 'KILLPARTIAL ', killfiles[ikill]
            rmfile, killfiles
         endif
      endif
   endfor

   ; Trim the plate list, and update NFILE to this trimmed number
   plist = plist[where(qkeep, nfile)]

   ;----------
   ; Make a list of one S/N for each plate which is the minimum of
   ; G1, I1, G2, I2.
   ; Assign a plate quality, but do not over-write any plate quality
   ; from the manually-assigned one in the "spPlateList.par" file.

   qualstring = ['bad', 'marginal', 'good']
   for ifile=0, nfile-1 do begin
      if (strtrim(plist[ifile].statuscombine,2) EQ 'Done') then begin
         plist[ifile].platesn2 = min( $
          [plist[ifile].sn2_g1, plist[ifile].sn2_i1, $
          plist[ifile].sn2_g2, plist[ifile].sn2_i2])
         iqual = 2
         if (plist[ifile].platesn2 LT 13) then iqual = iqual < 0
         if (plist[ifile].platesn2 LT 15) then iqual = iqual < 1
         if (plist[ifile].fbadpix GT 0.10) then iqual = iqual < 0
         if (plist[ifile].fbadpix GT 0.05) then iqual = iqual < 1
         if (NOT keyword_set(strtrim(plist[ifile].platequality))) then $
          plist[ifile].platequality = qualstring[iqual]
      endif
   endfor

   ;----------
   ; Decide which plates constitute unique tiles with the required S/N,
   ; then set QSURVEY=1.
   ; Also insist that PROGNAME='main'.

   ; First get the unique list of TILE
   isort = sort(plist.tile)
   isort = isort[ uniq(plist[isort].tile) ]
   tilelist = plist[isort].tile

   for itile=0, n_elements(tilelist)-1 do begin
      indx = where(plist.tile EQ tilelist[itile] $
       AND (strtrim(plist.platequality,2) EQ 'good' $
         OR strtrim(plist.platequality,2) EQ 'marginal') $
       AND strtrim(plist.progname,2) EQ 'main')
      if (indx[0] NE -1) then begin
         snbest = max(plist[indx].platesn2, ibest)
         plist[indx[ibest]].qsurvey = 1
      endif
   endfor

   ;---------------------------------------------------------------------------
   ; Write ASCII + HTML output files

   alias = [['CHUNKNAME'    , 'CHUNK'   ], $
            ['PROGNAME'     , 'PROG'    ], $
            ['PLATESN2'     , 'SN^2'    ], $
            ['N_GALAXY'     , 'N_gal'   ], $
            ['N_QSO'        , 'N_QSO'   ], $
            ['N_STAR'       , 'N_star'  ], $
            ['N_UNKNOWN'    , 'N_unk'   ], $
            ['N_SKY'        , 'N_sky'   ], $
            ['FBADPIX'      , 'Badpix'  ], $
            ['SUCCESS_MAIN' , '%Main'   ], $
            ['SUCCESS_LRG'  , '%LRG'    ], $
            ['SUCCESS_QSO'  , '%QSO'    ], $
            ['STATUS2D'     , '2D'      ], $
            ['STATUSCOMBINE', 'Combine' ], $
            ['STATUS1D'     , '1D'      ], $
            ['PLATEQUALITY' , 'QUALITY' ] ]

    platelist_write, plist, trimtags=trimtags1, alias=alias, $
    fileprefix='platelist', title='SDSS Spectroscopy Plates Observed List'
    platelist_write, plist, trimtags=trimtags2, alias=alias, $
    fileprefix='platequality', title='SDSS Spectroscopy Plate Quality List'

   ;----------
   ; Write the FITS binary table

   mwrfits, plist, fitsfile, /create

end
;------------------------------------------------------------------------------
