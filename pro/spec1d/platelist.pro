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
;   infile      - Either a list of combine-plan files or a list of plate files,
;                 which can contain wildcards.
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
;   The following files are generated:
;     $SPECTRO_DATA/platelist.fits
;     $SPECTRO_DATA/platelist.txt
;     $SPECTRO_DATA/platelist.html
;     $SPECTRO_DATA/platelist-mjdsort.txt
;     $SPECTRO_DATA/platelist-mjdsort.html
;     $SPECTRO_DATA/platequality.txt
;     $SPECTRO_DATA/platequality.html
;     $SPECTRO_DATA/platequality-mjdsort.txt
;     $SPECTRO_DATA/platequalitymjdsort.html
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
;   PLATEQUALITY is set to 'marginal' if min(NEXP_*) < 3
;                          'bad'      if min(NEXP_*) < 2
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
;   $IDLSPEC2D_DIR/etc/spPlateList.par
;
; PROCEDURES CALLED:
;   apo_checklimits()
;   chunkinfo()
;   copy_struct_inx
;   djs_diff_angle()
;   djs_filepath()
;   fileandpath()
;   headfits()
;   mrdfits()
;   repstr()
;   rmfile
;   splog
;   struct_print
;   sxpar()
;   tai2airmass()
;   yanny_free
;   yanny_par()
;   yanny_read
;   yanny_readone()
;
; INTERNAL SUPPORT ROUTINES:
;   platelist_write
;
; REVISION HISTORY:
;   29-Oct-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
PRO platelist_write, plist, trimtags=trimtags, alias=alias, $
    fileprefix=fileprefix, title=title, toptext=toptext
    ;
    ;
    ;
    ascfile = djs_filepath(fileprefix+'.txt', root_dir=GETENV('SPECTRO_DATA'))
    htmlfile = djs_filepath(fileprefix+'.html', root_dir=GETENV('SPECTRO_DATA'))
    trimdat  = struct_trimtags(plist, select_tags=trimtags[0,*])
    trimstring = struct_trimtags(plist, select_tags=trimtags[0,*], $
        format=trimtags[1,*])
    struct_print, trimstring, filename=ascfile, fdigit=3, alias=alias
    FOR itag=0, N_TAGS(trimdat)-1 DO BEGIN
        FOR iarr=0, N_ELEMENTS(trimdat[0].(itag))-1 DO BEGIN
            FOR irow=0, N_ELEMENTS(trimdat)-1 DO BEGIN
                markstring = apo_checklimits('SUMMARY', $
                    STRUPCASE(trimtags[0,itag]), '', $
                    trimdat[irow].(itag)[iarr], /html)
                trimstring[irow].(itag)[iarr] = markstring $
                    + trimstring[irow].(itag)[iarr]
                IF STRMATCH(markstring,'<span*') THEN $
                    trimstring[irow].(itag)[iarr] = $
                    trimstring[irow].(itag)[iarr] + '</span>'
            ENDFOR
        ENDFOR
    ENDFOR
    struct_print, trimstring, /html, alias=alias, tarray=tarray, css=css
    OPENW, lun, htmlfile, /GET_LUN
    PRINTF, lun, '<?xml version="1.0" encoding="UTF-8"?>'
    PRINTF, lun, '<!DOCTYPE html'
    PRINTF, lun, '     PUBLIC "-//W3C//DTD XHTML 1.1//EN"'
    PRINTF, lun, '     "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">'
    PRINTF, lun, '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">'
    PRINTF, lun, '<head>'
    PRINTF, lun, '<title>' + title + '</title>'
    PRINTF, lun, '<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />'
    FOR c=0, N_ELEMENTS(css)-1 DO $
        PRINTF, lun, css[c]
    PRINTF, lun, '</head>'
    PRINTF, lun, '<body>'
    PRINTF, lun, '<h1>' + title + '</h1>'
    FOR i=0, N_ELEMENTS(toptext)-1 DO $
        PRINTF, lun, toptext[i]
    FOR iline=0, N_ELEMENTS(tarray)-1 DO $
        PRINTF, lun, tarray[iline]
    PRINTF, lun, '</body>'
    PRINTF, lun, '</html>'
    CLOSE, lun
    FREE_LUN, lun
    RETURN
END
;------------------------------------------------------------------------------
PRO platelist, infile, plist=plist, create=create, $
    purge2d=purge2d, purge1d=purge1d, killpartial=killpartial
    ;
    ;
    ;
    fitsfile = djs_filepath('platelist.fits', root_dir=GETENV('SPECTRO_DATA'))
    ;
    ; If the /CREATE flag is not set, and the platelist file already exists
    ; on disk, then simply return the info in that file.
    ;
    If (~KEYWORD_SET(create) && ~KEYWORD_SET(purge2d) $
        && ~KEYWORD_SET(purge1d)) THEN BEGIN
        thisfile = (FINDFILE(fitsfile))[0]
        IF (KEYWORD_SET(thisfile)) THEN BEGIN
            plist = mrdfits(thisfile,1)
            RETURN
        ENDIF
    ENDIF
    ;
    ; Generate the list of plan files or plate files if not specified
    ;
    IF (KEYWORD_SET(infile)) THEN BEGIN
        FOR i=0L, N_ELEMENTS(infile)-1L DO BEGIN
            thisfile = FINDFILE(infile[i], COUNT=ct)
            IF (ct GT 0) THEN $
                fullfile = KEYWORD_SET(fullfile) ? [fullfile, thisfile] : thisfile
        ENDFOR
    ENDIF ELSE BEGIN
        dirlist = get_mjd_dir( GETENV('SPECTRO_DATA'), mjstart=1, mjend=9999)
        FOR i=0L, N_ELEMENTS(dirlist)-1L DO BEGIN
            thisfile = FINDFILE(djs_filepath('spPlancomb-*.par', $
                root_dir=GETENV('SPECTRO_DATA'), subdir=dirlist[i]), COUNT=ct)
            IF (ct EQ 0) THEN $
                thisfile = FINDFILE(djs_filepath('spPlate-*.fits', $
                root_dir=GETENV('SPECTRO_DATA'), subdir=dirlist[i]), COUNT=ct)
            IF (ct GT 0) THEN $
                fullfile = KEYWORD_SET(fullfile) ? [fullfile, thisfile] : thisfile
        ENDFOR
    ENDELSE
    nfile = N_ELEMENTS(fullfile)
    IF (nfile EQ 0) THEN RETURN
    ;
    ; Sort these files
    ;
    fullfile = fullfile[SORT(fullfile)]
    ;
    ; Create output structure
    ;
    plist = CREATE_STRUCT( $
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
        'nexp_b1'  , 0L,  $
        'nexp_b2'  , 0L,  $
        'nexp_r1'  , 0L,  $
        'nexp_r2'  , 0L,  $
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
        'goffstd'  , 0., $
        'grmsstd'  , 0., $
        'roffstd'  , 0., $
        'rrmsstd'  , 0., $
        'ioffstd'  , 0., $
        'irmsstd'  , 0., $
        'groffstd' , 0., $
        'grrmsstd' , 0., $
        'rioffstd' , 0., $
        'rirmsstd' , 0., $
        'goffgal'  , 0., $
        'grmsgal'  , 0., $
        'roffgal'  , 0., $
        'rrmsgal'  , 0., $
        'ioffgal'  , 0., $
        'irmsgal'  , 0., $
        'groffgal' , 0., $
        'grrmsgal' , 0., $
        'rioffgal' , 0., $
        'rirmsgal' , 0., $
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
    plist = REPLICATE(plist, nfile)
    ;
    ; Read the data file with the public plate information
    ;
    publicfile = FILEPATH('spPlateList.par', $
        ROOT_DIR=GETENV('IDLSPEC2D_DIR'), SUBDIRECTORY='etc')
    publicdata = yanny_readone(publicfile, 'SPPLATELIST')
    ;
    ; Loop through all files
    ;
    platefile = STRARR(nfile)
    combparfile = STRARR(nfile)
    comblogfile = STRARR(nfile)
    combpsfile = STRARR(nfile)
    zlogfile = STRARR(nfile)
    zbestfile = STRARR(nfile)
    zallfile = STRARR(nfile)
    FOR ifile=0, nfile-1 DO BEGIN
        splog, 'Looking at ' + fullfile[ifile]
        ;
        ; Determine PATH
        ;
        junk = fileandpath(fullfile[ifile], path=path)
        ;
        ; Test if INFILE specifies Yanny param files for spPlancomb.
        ;
        IF (STRMID(fullfile[ifile],STRLEN(fullfile[ifile])-4) EQ '.par') $
            THEN BEGIN
            combparfile[ifile] = fullfile[ifile]
            yanny_read, fullfile[ifile], hdr=hdrp
            platefile[ifile] = $
                djs_filepath(yanny_par(hdrp, 'combinefile'), root_dir=path)
        ENDIF ELSE BEGIN
            platefile[ifile] = fullfile[ifile]
            combparfile[ifile] = repstr(platefile[ifile], 'spPlate', 'spPlancomb')
            combparfile[ifile] = repstr(combparfile[ifile], '.fits', '.par')
        ENDELSE
        ;
        ; Determine names of associated files
        ;
        comblogfile[ifile] = repstr(combparfile[ifile], '.par', '.log')
        comblogfile[ifile] = repstr(comblogfile[ifile], 'spPlancomb', 'spDiagcomb')
        combpsfile[ifile] = repstr(comblogfile[ifile], '.log', '.ps')
        platemjd = STRMID(fileandpath(platefile[ifile]), 8, 10)
        zbestfile[ifile] = djs_filepath('spZbest-' + platemjd + '.fits', $
            root_dir=path)
        zallfile[ifile] = djs_filepath('spZall-' + platemjd + '.fits', $
            root_dir=path)
        ;
        ; Read the combine plan file to get the list of all the 2D plan files
        ; from its Yanny header.
        ; Also get the mapping name from the combine par file in case we were
        ; unable to get it from the spPlate file.
        ;
        plist[ifile].mapname = (yanny_readone(combparfile[ifile], 'SPEXP', $
            hdr=hdrcomb))[0].mapname
        ;
        ; Find the state of the 2D reductions (not the combine step)
        ;
        statusdone = 0
        statusrun = 0
        statusmissing = 0
        ;
        ; Check status of individual 2D runs
        ;
        planlist = yanny_par(hdrcomb, 'planfile2d') ; Assume we find this
        planlist = djs_filepath(planlist, root_dir=path)
        logfile2d = '' ; List of 2D log files that exist
        FOR iplan=0, N_ELEMENTS(planlist)-1 DO BEGIN
            yanny_read, planlist[iplan], hdr=hdr2d
            plist[ifile].mjdlist = STRTRIM(plist[ifile].mjdlist $
                + ' ' + yanny_par(hdr2d, 'MJD'),2)
            thislogfile = djs_filepath(yanny_par(hdr2d, 'logfile'), root_dir=path)
            thislogfile = (FINDFILE(thislogfile))[0]
            IF KEYWORD_SET(thislogfile) THEN BEGIN
                IF ~KEYWORD_SET(logfile2d) THEN logfile2d = thislogfile $
                ELSE logfile2d = [logfile2d, thislogfile]
                SPAWN, 'tail -1 '+thislogfile, lastline
                IF (STRMATCH(lastline[0], '*Successful completion*')) THEN BEGIN
                    ; Case where this 2D log file completed
                    statusdone = statusdone + 1
                ENDIF ELSE BEGIN
                    ; Case where this 2D log file isn't completed
                    statusrun = statusrun + 1
                ENDELSE
            ENDIF ELSE BEGIN
                ; Case where this 2D log file missing
                statusmissing = statusmissing + 1
            ENDELSE
        ENDFOR
        IF (statusmissing GT 0 AND statusrun EQ 0) THEN BEGIN
            plist[ifile].status2d = 'Pending'
        ENDIF ELSE IF (statusmissing GT 0 OR statusrun GT 0) THEN BEGIN
            plist[ifile].status2d = 'RUNNING'
        ENDIF ELSE BEGIN
             plist[ifile].status2d = 'Done'
        ENDELSE
        ;
        ; Read plate file - get status of Combine
        ;
        hdr1 = headfits(platefile[ifile], /silent, errmsg=errmsg)
        IF (SIZE(hdr1, /TNAME) EQ 'STRING') THEN BEGIN
            ; plist[ifile].plate = sxpar(hdr1, 'PLATEID')
            ; plist[ifile].mjd = sxpar(hdr1, 'MJD')
            plist[ifile].mjdlist = sxpar(hdr1, 'MJDLIST')
            plist[ifile].tile = sxpar(hdr1, 'TILEID')
            plist[ifile].ra = sxpar(hdr1, 'RA')
            plist[ifile].dec = sxpar(hdr1, 'DEC')
            plist[ifile].cartid = sxpar(hdr1, 'CARTID')
            plist[ifile].tai = sxpar(hdr1, 'TAI')
            plist[ifile].tai_beg = sxpar(hdr1, 'TAI-BEG')
            plist[ifile].tai_end = sxpar(hdr1, 'TAI-END')
            plist[ifile].airmass = sxpar(hdr1, 'AIRMASS')
            ; plist[ifile].airmass = tai2airmass(plist[ifile].ra, $
            ;     plist[ifile].dec, tai=plist[ifile].tai)
            plist[ifile].exptime = sxpar(hdr1, 'EXPTIME')
            plist[ifile].nexp = sxpar(hdr1, 'NEXP')
            plist[ifile].nexp_b1 = sxpar(hdr1, 'NEXP_B1')
            plist[ifile].nexp_b2 = sxpar(hdr1, 'NEXP_B2')
            plist[ifile].nexp_r1 = sxpar(hdr1, 'NEXP_R1')
            plist[ifile].nexp_r2 = sxpar(hdr1, 'NEXP_R2')
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
            plist[ifile].goffstd = sxpar(hdr1, 'GOFFSTD')
            plist[ifile].grmsstd = sxpar(hdr1, 'GRMSSTD')
            plist[ifile].roffstd = sxpar(hdr1, 'ROFFSTD')
            plist[ifile].rrmsstd = sxpar(hdr1, 'RRMSSTD')
            plist[ifile].ioffstd = sxpar(hdr1, 'IOFFSTD')
            plist[ifile].irmsstd = sxpar(hdr1, 'IRMSSTD')
            plist[ifile].groffstd = sxpar(hdr1, 'GROFFSTD')
            plist[ifile].grrmsstd = sxpar(hdr1, 'GRRMSSTD')
            plist[ifile].rioffstd = sxpar(hdr1, 'RIOFFSTD')
            plist[ifile].rirmsstd = sxpar(hdr1, 'RIRMSSTD')
            plist[ifile].goffgal = sxpar(hdr1, 'GOFFGAL')
            plist[ifile].grmsgal = sxpar(hdr1, 'GRMSGAL')
            plist[ifile].roffgal = sxpar(hdr1, 'ROFFGAL')
            plist[ifile].rrmsgal = sxpar(hdr1, 'RRMSGAL')
            plist[ifile].ioffgal = sxpar(hdr1, 'IOFFGAL')
            plist[ifile].irmsgal = sxpar(hdr1, 'IRMSGAL')
            plist[ifile].groffgal = sxpar(hdr1, 'GROFFGAL')
            plist[ifile].grrmsgal = sxpar(hdr1, 'GRRMSGAL')
            plist[ifile].rioffgal = sxpar(hdr1, 'RIOFFGAL')
            plist[ifile].rirmsgal = sxpar(hdr1, 'RIRMSGAL')
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
            plist[ifile].mapname = STRTRIM(sxpar(hdr1, 'NAME'))
            plist[ifile].vers2d = STRTRIM(sxpar(hdr1, 'VERS2D'))
            plist[ifile].verscomb = STRTRIM(sxpar(hdr1, 'VERSCOMB'))
            plist[ifile].statuscombine = 'Done'
        ENDIF ELSE BEGIN
            ;
            ; Case where no spPlate file exists
            ;
            thislogfile = djs_filepath(yanny_par(hdrcomb, 'logfile'), root_dir=path)
            thislogfile = (FINDFILE(thislogfile))[0]
            IF KEYWORD_SET(thislogfile) THEN BEGIN
                SPAWN, 'tail -1 '+thislogfile, lastline
                IF STRMATCH(lastline[0], '*Successful completion*') THEN BEGIN
                    ; Case where this combine log file completed
                    plist[ifile].statuscombine = 'FAILED'
                ENDIF ELSE BEGIN
                    ; Case where this combine log file isn't completed
                    SPAWN, 'grep ABORT '+thislogfile, abortline
                    IF KEYWORD_SET(abortline) THEN BEGIN
                        plist[ifile].statuscombine = 'FAILED' ; Combining step aborted
                    ENDIF ELSE BEGIN
                        plist[ifile].statuscombine = 'RUNNING'
                    ENDELSE
                ENDELSE
            ENDIF ELSE BEGIN
                ; Case where this combine log file missing
                plist[ifile].statuscombine = 'Pending'
            ENDELSE
            IF (KEYWORD_SET(purge2d) $
                AND plist[ifile].statuscombine NE 'Done') THEN BEGIN
                splog, 'PURGE2D ', logfile2d
                rmfile, logfile2d
                splog, 'PURGE2D ', comblogfile[ifile]
                rmfile, comblogfile[ifile]
                splog, 'PURGE2D ', combpsfile[ifile]
                rmfile, combpsfile[ifile]
                plist[ifile].status2d = 'Pending'
                plist[ifile].statuscombine = 'Pending'
            ENDIF
        ENDELSE
        ;
        ; Get the following from the file names, since sometimes they
        ; are wrong in the file headers!!
        ;
        plist[ifile].plate = LONG( STRMID(fileandpath(platefile[ifile]), 8, 4) )
        plist[ifile].mjd = LONG( STRMID(fileandpath(platefile[ifile]), 13, 5) )
        ;
        ; Determine the chunk name and the version of target used
        ;
        cinfo = chunkinfo(plist[ifile].plate)
        plist[ifile].chunkname = cinfo.chunkname
        plist[ifile].progname = cinfo.progname
        ;
        ; Determine which public data release has this plate+MJD
        ;
        j = WHERE(plist[ifile].plate EQ publicdata.plate $
             AND plist[ifile].mjd EQ publicdata.mjd)
        IF (j[0] NE -1) THEN BEGIN
            copy_struct_inx, publicdata[j[0]], plist, index_to=ifile
        ENDIF
        ;
        ; The RA,DEC in the header is sometimes wrong, so try to derive
        ; the field center from the plug-map information.  Choose the
        ; coordinates of the object closest to the center of the plate
        ; as defined by XFOCAL=YFOCAL=0.  Replace the plate center with
        ; the position of that object in the cases where they disagree
        ; by more than 1.5 degrees.
        ;
        IF KEYWORD_SET(FINDFILE(platefile[ifile])) THEN $
            plug = mrdfits(platefile[ifile], 5, /silent) $
        ELSE $
            plug = 0
        IF KEYWORD_SET(plug) THEN BEGIN
            iobj = WHERE(STRTRIM(plug.holetype,2) EQ 'OBJECT')
            junk = MIN( plug[iobj].xfocal^2 + plug[iobj].yfocal^2, imin)
            thisra = plug[iobj[imin]].ra
            thisdec = plug[iobj[imin]].dec
            adist = djs_diff_angle(plist[ifile].ra, plist[ifile].dec, $
                thisra, thisdec)
            IF (adist GT 1.5) THEN BEGIN
                splog, 'WARNING: Replacing plate center for plate', $
                    plist[ifile].plate
                plist[ifile].ra = thisra
                plist[ifile].dec = thisdec
                plist[ifile].airmass = tai2airmass(plist[ifile].ra, $
                    plist[ifile].dec, tai=plist[ifile].tai)
            ENDIF
        ENDIF
        ;
        ; Read Zbest file - get status of 1D
        ;
        IF KEYWORD_SET(zbestfile[ifile]) THEN $
            hdr2 = headfits(zbestfile[ifile], /silent, errmsg=errmsg) $
        ELSE $
            hdr2 = 0
        IF (SIZE(hdr2, /TNAME) EQ 'STRING') THEN BEGIN
            zans = mrdfits(zbestfile[ifile], 1, /silent)
            plug = mrdfits(platefile[ifile], 5, /silent)
            class = STRTRIM(zans.class,2)
            ; Use the ZWARNING flag if it exists to identify SKY or UNKNOWN.
            IF ((WHERE(TAG_NAMES(zans) EQ 'ZWARNING'))[0] NE -1) THEN $
                zwarning = zans.zwarning $
            ELSE $
                zwarning = BYTARR(N_ELEMENTS(zans))
            qsky = (zwarning AND 1) NE 0
            plist[ifile].n_galaxy = TOTAL(class EQ 'GALAXY' AND zwarning EQ 0)
            plist[ifile].n_qso = TOTAL(class EQ 'QSO' AND zwarning EQ 0)
            plist[ifile].n_star = TOTAL(class EQ 'STAR' AND zwarning EQ 0)
            plist[ifile].n_unknown = TOTAL(class EQ 'UNKNOWN' $
                OR (zwarning NE 0 AND qsky EQ 0))
            plist[ifile].n_sky = TOTAL(class EQ 'SKY' OR qsky EQ 1)
            plist[ifile].vers1d = STRTRIM(sxpar(hdr2, 'VERS1D'))
            plist[ifile].status1d = 'Done'
            nobj = N_ELEMENTS(zans)
            targets = STRARR(nobj)
            FOR iobj=0, nobj-1 DO $
                targets[iobj] = sdss_flagname('TARGET',plug[iobj].primtarget, $
                    /silent, /concat)+' '
            imain = WHERE(STRMATCH(targets,'*GALAXY *') $
                OR STRMATCH(targets,'*GALAXY_BIG *') $
                OR STRMATCH(targets,'*GALAXY_BRIGHT_CORE *'), nmain)
            ilrg = WHERE(STRMATCH(targets,'*GALAXY_RED *') $
                OR STRMATCH(targets,'*GALAXY_RED_II *'), nlrg)
            iqso = WHERE(STRMATCH(targets,'*QSO_HIZ *') $
                OR STRMATCH(targets,'*QSO_CAP *') $
                OR STRMATCH(targets,'*QSO_SKIRT *') $
                OR STRMATCH(targets,'*QSO_FIRST_CAP *') $
                OR STRMATCH(targets,'*QSO_FIRST_SKIRT *'), nqso)
            plist[ifile].n_target_main = nmain
            plist[ifile].n_target_lrg = nlrg
            plist[ifile].n_target_qso = nqso
            IF (nmain GT 0) THEN $
                plist[ifile].success_main = $
                    100 * TOTAL(zans[imain].zwarning EQ 0 $
                    AND (STRMATCH(zans[imain].class,'GALAXY*') $
                    OR STRMATCH(zans[imain].class,'QSO*'))) / nmain
            IF (nlrg GT 0) THEN $
                plist[ifile].success_lrg = $
                    100 * TOTAL(zans[ilrg].zwarning EQ 0 $
                    AND STRMATCH(zans[ilrg].class,'GALAXY*')) / nlrg
            IF (nqso GT 0) THEN $
                plist[ifile].success_qso = $
                    100 * TOTAL(zans[iqso].zwarning EQ 0 $
                    AND STRMATCH(zans[iqso].class,'QSO*')) / nqso
        ENDIF ELSE BEGIN
            ;
            ; Find the state of the 1D reductions -- spZbest file is missing
            ;
            zlogfile[ifile] = repstr(fileandpath(platefile[ifile]), 'spPlate', 'spDiag1d')
            zlogfile[ifile] = repstr(zlogfile[ifile], '.fits', '.log')
            zlogfile[ifile] = djs_filepath(zlogfile[ifile], root_dir=path)
            zlogfile[ifile] = (FINDFILE(zlogfile[ifile]))[0]
            IF KEYWORD_SET(zlogfile[ifile]) THEN BEGIN
                SPAWN, 'tail -1 '+zlogfile[ifile], lastline
                IF STRMATCH(lastline[0], '*Successful completion*') THEN BEGIN
                    ; Case where this 1D log file completed
                    plist[ifile].status1d = 'FAILED'; Should have found spZbest file
                ENDIF ELSE BEGIN
                    ; Case where this 1D log file isn't completed
                    IF KEYWORD_SET(purge1d) THEN BEGIN
                        splog, 'PURGE1D ', zlogfile[ifile]
                        rmfile, zlogfile[ifile]
                        plist[ifile].status1d = 'Pending'
                    ENDIF ELSE BEGIN
                        plist[ifile].status1d = 'RUNNING'
                    ENDELSE
                ENDELSE
            ENDIF ELSE BEGIN
                plist[ifile].status1d = 'Pending'
            ENDELSE
        ENDELSE
    ENDFOR
    ;
    ; Remove from the plate list earlier reductions of the same plugging
    ; (keeping only the most recent MJD of each plugging).
    ;
    qkeep = BYTARR(nfile)
    ;
    ; First get the unique list of MAPNAME, then mark the most recent MJD
    ; of each as the good one.
    ;
    isort = SORT(plist.mapname)
    isort = isort[ UNIQ(plist[isort].mapname) ]
    maplist = plist[isort].mapname
    FOR imp=0, N_ELEMENTS(maplist)-1 DO BEGIN
        indx = WHERE(plist.mapname EQ maplist[imp])
        junk = MAX(plist[indx].mjd, imax)
        qkeep[indx[imax]] = 1
    ENDFOR
    ;
    ; Don't discard any where MAPNAME isn't set
    ;
    indx = WHERE(STRTRIM(plist.mapname) EQ '')
    IF (indx[0] NE -1) THEN qkeep[indx] = 1
    ;
    ; List partially-combined plates that we're discarding from the list
    ;
    FOR ifile=0, nfile-1 DO BEGIN
        IF (qkeep[ifile] NE 1) THEN BEGIN
            splog, 'Discard partially-combined ' + combparfile[ifile]
            IF (KEYWORD_SET(killpartial)) THEN BEGIN
                ; Also should be purging the spFluxcalib files and spSN2d files !!!???
                killfiles = [ combparfile[ifile], $
                    comblogfile[ifile], $
                    combpsfile[ifile], $
                    platefile[ifile], $
                    zbestfile[ifile], $
                    zallfile[ifile], $
                    zlogfile[ifile] ]
                FOR ikill=0, N_ELEMENTS(killfiles)-1 DO $
                    splog, 'KILLPARTIAL ', killfiles[ikill]
                rmfile, killfiles
            ENDIF
        ENDIF
    ENDFOR
    ;
    ; Trim the plate list, and update NFILE to this trimmed number
    ;
    plist = plist[WHERE(qkeep, nfile)]
    ;
    ; Make a list of one S/N for each plate which is the minimum of
    ; G1, I1, G2, I2.
    ; Assign a plate quality, but do not over-write any plate quality
    ; from the manually-assigned one in the "spPlateList.par" file.
    ;
    qualstring = ['bad', 'marginal', 'good']
    FOR ifile=0, nfile-1 DO BEGIN
        IF (STRTRIM(plist[ifile].statuscombine,2) EQ 'Done') THEN BEGIN
            nexp_min = MIN( $
                [plist[ifile].nexp_b1, plist[ifile].nexp_b2, $
                plist[ifile].nexp_r1, plist[ifile].nexp_r2], MAX=nexp_max)
            plist[ifile].platesn2 = MIN( $
                [plist[ifile].sn2_g1, plist[ifile].sn2_i1, $
                plist[ifile].sn2_g2, plist[ifile].sn2_i2])
            iqual = 2
            IF (plist[ifile].platesn2 LT 13) THEN iqual = iqual < 0
            IF (plist[ifile].platesn2 LT 15) THEN iqual = iqual < 1
            IF (plist[ifile].fbadpix GT 0.10) THEN iqual = iqual < 0
            IF (plist[ifile].fbadpix GT 0.05) THEN iqual = iqual < 1
            ; For reductions before v5_1, NEXP_MIN and NEXP_MAX are always zero
            IF (nexp_max GT 0) THEN BEGIN
                IF (nexp_min LT 2) THEN iqual = iqual < 0
                IF (nexp_min LT 3) THEN iqual = iqual < 1
            ENDIF
            IF ~KEYWORD_SET(STRTRIM(plist[ifile].platequality)) THEN $
                plist[ifile].platequality = qualstring[iqual]
        ENDIF
    ENDFOR
    ;
    ; Decide which plates constitute unique tiles with the required S/N,
    ; then set QSURVEY=1.
    ; Also insist that PROGNAME='main'.
    ;
    ; First get the unique list of TILE
    ;
    isort = SORT(plist.tile)
    isort = isort[ UNIQ(plist[isort].tile) ]
    tilelist = plist[isort].tile
    FOR itile=0, N_ELEMENTS(tilelist)-1 DO BEGIN
        indx = WHERE(plist.tile EQ tilelist[itile] $
            AND (STRTRIM(plist.platequality,2) EQ 'good' $
            OR STRTRIM(plist.platequality,2) EQ 'marginal') $
            AND STRTRIM(plist.progname,2) EQ 'main')
        IF (indx[0] NE -1) THEN BEGIN
            snbest = max(plist[indx].platesn2, ibest)
            plist[indx[ibest]].qsurvey = 1
        ENDIF
    ENDFOR
    ;
    ; Write ASCII + HTML output files
    ;
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
    isort = REVERSE(SORT(plist.mjd))
    toptext = [ '<p>Last Update: '+ SYSTIME()+'</p>', $
        '<ul>', $
        '<li><a href="http://spectro.princeton.edu/">HOME</a></li>', $
        '<li>Sorted by plate: <a href="platelist.html">HTML</a>' $
        + ' <a href="platelist.txt">ASCII</a></li>', $
        '<li>Sorted by MJD: <a href="platelist-mjdsort.html">HTML</a>' $
        + ' <a href="platelist-mjdsort.txt">ASCII</a></li>', $
        '<li><a href="platelist.fits">FITS version</a></li>', $
        '</ul>']
    platelist_write, plist, trimtags=trimtags1, alias=alias, $
        fileprefix='platelist', toptext=toptext, $
        title='SDSS Spectroscopy Plates Observed List'
    platelist_write, plist[isort], trimtags=trimtags1, alias=alias, $
        fileprefix='platelist-mjdsort', toptext=toptext, $
        title='SDSS Spectroscopy Plates Observed List'
    toptext = [ '<p>Last Update: '+ SYSTIME()+'</p>', $
        '<ul>', $
        '<li><a href="http://spectro.princeton.edu/">HOME</a></li>', $
        '<li>Sorted by plate: <a href="platequality.html">HTML</a>' $
        + ' <a href="platequality.txt">ASCII</a></li>', $
        '<li>Sorted by MJD: <a href="platequality-mjdsort.html">HTML</a>' $
        + ' <a href="platequality-mjdsort.txt">ASCII</a></li>', $
        '<li><a href="platelist.fits">FITS version</a></li>', $
        '</ul>']
    platelist_write, plist, trimtags=trimtags2, alias=alias, $
        fileprefix='platequality', toptext=toptext, $
        title='SDSS Spectroscopy Plate Quality List'
    platelist_write, plist[isort], trimtags=trimtags2, alias=alias, $
        fileprefix='platequality-mjdsort', toptext=toptext, $
        title='SDSS Spectroscopy Plate Quality List'
    ;
    ; Write the FITS binary table
    ;
    mwrfits, plist, fitsfile, /create
    RETURN
END
;------------------------------------------------------------------------------
