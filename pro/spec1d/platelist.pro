;+
; NAME:
;   platelist
;
; PURPOSE:
;   Make list of reduced plates
;
; CALLING SEQUENCE:
;   platelist, [ /create, topdir=, outdir=, run2d=, run1d=, $
;    /purge2d, /purge1d, /killpartial, skipcart=, /rawsn2, plist= ]
;
; INPUTS
;
; OPTIONAL INPUTS:
;   create      - If set, then re-generate the "platelist.fits" file;
;                 if not set, then simply read this file from a previous call.
;   topdir      - Optional override value for the environment
;                 variable $BOSS_SPECTRO_REDUX.
;   outdir      - Optional override of both topdir and $BOSS_SPECTRO_REDUX
;                 for output directory location
;   run2d       - Optional RUN2D subdirectories to include in outputs;
;                 set to '' to not search subdirectories;
;                 set to '*' to search all subdirs; default to $RUN2D
;   run1d       - Optional RUN1D subdirectories to include in outputs
;                 set to '' to not search subdirectories;
;                 set to '*' to search all subdirs; default to $RUN1D
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
;   skipcart    - cart number or list of cart numbers to drop from platelist
;   rawsn2      - If set, output original raw SN2 numbers in html and text
;                 files; otherwise use dereddened (dust extinction corrected)
;                 values.  Both are always written to the FITS file output.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   plist       - Output structure with information for each plate.
;
; COMMENTS:
;   RUN2D: If option isn't set, use $RUN2D env var
;   RUN2D/RUN1D = '*' means seach all subdirs
;   Output directory:
;     If OUTDIR is set, use that.
;     Otherwise, if a single RUN2D is given, use $BOSS_SPECTRO_REDUX/[RUN2D]/
;        --> Note: RUN2D option, not $RUN2D environment variable
;     Otherwise (none or multiple RUN2D), use $BOSS_SPECTRO_REDUX/
;     Option TOPDIR can override the value of $BOSS_SPECTRO_REDUX
;        in the above, but it won't override OUTDIR
;           
;   TOPDIR is used to override $BOSS_SPECTRO_REDUX for both input and
;   output.  OUTDIR can be used to override TOPDIR for output if you want
;   to get input from one dir while writing output to a different dir.
;        
;   The following files are generated in the output directory
;     platelist.fits
;     platelist.txt
;     platelist.html
;     platelist-mjdsort.txt
;     platelist-mjdsort.html
;     platequality.txt
;     platequality.html
;     platequality-mjdsort.txt
;     platequalitymjdsort.html
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
;   DEREDSN2 like PLATESN2, but with dereddened SN2 values
;   PLATEQUALITY defaults to 'good'.
;   PLATEQUALITY is set to 'bad'      if MINSN2(B) < 10.0 or MINSN2(R) < 22.0
;                                     --> previously if MINSN2 < 13.0
;   PLATEQUALITY is set to 'bad'      if FBADPIX > 0.10
;   PLATEQUALITY is set to 'bad'      if min(NEXP_*) < 3
;
;   Decide which plates constitute unique tiles with the required S/N,
;   then set QSURVEY=1.  Require PLATEQUALITY='good'.
;   Also require PROGNAME='main'.
;
; EXAMPLES:
;
; BUGS:
;   Spawns the Unix commands 'tail' and 'grep', which is very slow.
;   If the spPlate file is missing, it doesn't try to figure out which
;        cart it really is and thus ignores skipcart
;
; DATA FILES:
;   $PLATELIST_DIR/platePlans.par
;   $SPECLOG_DIR/opfiles/spPlateList.par
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
;   11-Jan-2011  Stephen Bailey, LBNL
;     * Updated (S/N)^2 thresholds for "bad"
;     * Added get_lastline to be faster than spawn tail -1
;   21-Aug-2012  Stephen Bailey, LBNL
;     * Changed default output directory to RUN2D
;------------------------------------------------------------------------------

;----------
; get lastline of a file faster than spawning 'tail -1'
function get_lastline, filename
   ;;; spawn, 'tail -1 '+filename, lastline
   
   openr, ilun, filename, /get_lun
   lastline = ''
   while (NOT eof(ilun)) do begin
      readf, ilun, lastline
   endwhile
   close, ilun
   free_lun, ilun
   
   return, lastline
end

pro platelist_write, plist, trimtags=trimtags, alias=alias, $
 fileprefix=fileprefix, title=title, toptext=toptext, outdir=outdir1
 
   if (keyword_set(outdir1)) then begin
       outdir = outdir1
   endif else begin
       outdir = getenv('BOSS_SPECTRO_REDUX')
   endelse

   ascfile = djs_filepath(fileprefix+'.txt', root_dir=outdir)
   htmlfile = djs_filepath(fileprefix+'.html', root_dir=outdir)

   trimdat  = struct_trimtags(plist, select_tags=trimtags[0,*])
   trimstring = struct_trimtags(plist, select_tags=trimtags[0,*], $
    format=trimtags[1,*])
   struct_print, trimstring, filename=ascfile, fdigit=3, alias=alias

   for itag=0L, n_tags(trimdat)-1L do begin
      for iarr=0L, n_elements(trimdat[0].(itag))-1L do begin
         for irow=0L, n_elements(trimdat)-1L do begin
            markstring = apo_checklimits('SUMMARY', $
             strupcase(trimtags[0,itag]), '', $
             trimdat[irow].(itag)[iarr], /html)
            trimstring[irow].(itag)[iarr] = markstring $
             + trimstring[irow].(itag)[iarr]
            if strmatch(markstring,'<span*') then $
             trimstring[irow].(itag)[iarr] = $
             trimstring[irow].(itag)[iarr] + '</span>'
         endfor
      endfor
   endfor
   struct_print, trimstring, /html, alias=alias, tarray=tarray, css=css
   openw, lun, htmlfile, /get_lun
   printf, lun, '<?xml version="1.0" encoding="UTF-8"?>'
   printf, lun, '<!DOCTYPE html'
   printf, lun, '     PUBLIC "-//W3C//DTD XHTML 1.1//EN"'
   printf, lun, '     "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">'
   printf, lun, '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">'
   printf, lun, '<head>'
   printf, lun, '<title>' + title + '</title>'
   printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />'
   for c=0L, n_elements(css)-1L do $
    printf, lun, css[c]
   printf, lun, '</head>'
   printf, lun, '<body>'
   printf, lun, '<h1>' + title + '</h1>'
   for i=0L, n_elements(toptext)-1L do $
    printf, lun, toptext[i]

;   for iline=0, n_elements(tarray)-1 do $
;    printf, lun, tarray[iline]
   ifirst = where(plist.run2d NE shift(plist.run2d,1))
   for i=0, 2 do printf, lun, tarray[i]
   for iline=0L, n_elements(plist)-1L do begin
      ; Make a new html table for each new RUN2D
      if (iline NE 0 AND total(iline EQ ifirst) NE 0) then $
       printf, lun, tarray[1]
      printf, lun, tarray[3+iline]
   endfor
   printf, lun, '</table>'

   printf, lun, '</body>'
   printf, lun, '</html>'
   close, lun
   free_lun, lun

   return
end
;------------------------------------------------------------------------------
pro platelist, plist=plist, create=create, $
 topdir=topdir1, outdir=outdir1, $
 run2d=run2d1, run1d=run1d1, $
 purge2d=purge2d, purge1d=purge1d, killpartial=killpartial, $
 skipcart=skipcart, rawsn2=rawsn2

   if (n_elements(run2d1) GT 0) then run2d = strtrim(run2d1,2) $
    else run2d = getenv('RUN2D')
    
   if (n_elements(run1d1) GT 0) then begin
      run1d = strtrim(run1d1,2)
   endif else begin
      if (n_elements(run2d1) EQ 1) then begin
         run1d = run2d
      endif else begin
         if (n_elements(run2d) GT 1) then begin
            run1d = '*'
         endif else begin
            run1d = getenv('RUN1D')
         endelse
      endelse
   endelse

   if (keyword_set(topdir1)) then topdir = topdir1 $
    else topdir = getenv('BOSS_SPECTRO_REDUX')

   if (keyword_set(outdir1)) then begin
       outdir = outdir1
       if strmid(outdir, 0, 1) NE '/' then begin
           cd, current=current_dir
           outdir = current_dir + '/' + outdir
       endif
   endif else begin
      ; If exactly one run2d is given, use that subdir
      if (n_elements(run2d) EQ 1) && (run2d NE '*') then $
       outdir = topdir + '/' + run2d $
      else $
       outdir = topdir
   endelse
   
   if keyword_set(skipcart) then begin
      splog, 'WARNING: Dropping plates from carts', skipcart
   endif

   ;----------
   ; If the /CREATE flag is not set, and the platelist file already exists
   ; on disk, then simply return the info in that file.

   fitsfile = djs_filepath('platelist.fits', root_dir=outdir)

   if (NOT keyword_set(create) AND NOT keyword_set(purge2d) $
    AND NOT keyword_set(purge1d)) then begin
      thisfile = (findfile(fitsfile))[0]
      if (keyword_set(thisfile)) then begin
         plist = mrdfits(thisfile, 1, /silent)
         return
      endif
   endif

   ;----------
   ; Generate the list of plan files or plate files if not specified

   fitsfile = djs_filepath('platelist.fits', root_dir=outdir)

   if (keyword_set(run2d)) then $
    run2dlist = get_mjd_dir(topdir, mjd=run2d, /alldirs) $
   else $
    run2dlist = ''

   for j=0L, n_elements(run2dlist)-1L do begin
      if (keyword_set(run2dlist[j])) then $
       thisdir = djs_filepath('', root_dir=topdir, subdir=run2dlist[j]) $
      else $
       thisdir = topdir
      dirlist = get_mjd_dir(thisdir, mjstart=1, mjend=9999)
      if (keyword_set(dirlist)) then begin
         for i=0L, n_elements(dirlist)-1L do begin
            ; Select only those files matching the plate of the directory
            thisfile = findfile(djs_filepath('spPlancomb-'+dirlist[i]+'*.par', $
             root_dir=thisdir, subdir=dirlist[i]), count=ct)
            if (ct EQ 0) then $
             thisfile = findfile(djs_filepath('spPlate-'+dirlist[i]+'*.fits', $
              root_dir=thisdir, subdir=dirlist[i]), count=ct)
            if (ct GT 0) then begin
               if (keyword_set(fullfile)) then begin
                  fullfile = [fullfile, thisfile]
                  fullrun2d = [fullrun2d, replicate(run2dlist[j],ct)]
               endif else begin
                  fullfile = thisfile
                  fullrun2d = replicate(run2dlist[j],ct)
               endelse
            endif
         endfor
      endif
   endfor
   nfile = n_elements(fullfile)
   if (nfile EQ 0) then return

   fullfile = fullfile[sort(fullfile)] ; Sort these files

   ;----------
   ; Create output structure

   plist = create_struct( $
    'plate'        , 0L, $
    'tileid'       , 0L, $
    'mjd'          , 0L, $
    'run2d'        , '', $
    'run1d'        , '', $
    'racen'        , 0.0, $
    'deccen'       , 0.0, $
    'epoch'        , 0.0, $
    'cartid'       , 0L, $
    'tai'          , 0.0D, $
    'tai_beg'      , 0.0D, $
    'tai_end'      , 0.0D, $
    'airmass'      , 0.0, $
    'exptime'      , 0.0, $
    'mapname'      , ' ', $
    'survey'       , ' ', $
    'programname'  , ' ', $
    'chunk'        , ' ', $
    'platequality' , ' ', $
    'platesn2'     , 0.0, $
    'deredsn2'     , 0.0, $
    'qsurvey'      , 0L,  $
    'mjdlist'      , ' ', $
    'nexp'         , 0L,  $
    'nexp_b1'      , 0L,  $
    'nexp_b2'      , 0L,  $
    'nexp_r1'      , 0L,  $
    'nexp_r2'      , 0L,  $
    'expt_b1'      , 0.0, $
    'expt_b2'      , 0.0, $
    'expt_r1'      , 0.0, $
    'expt_r2'      , 0.0, $
    'sn2_g1'       , 0.0, $
    'sn2_r1'       , 0.0, $
    'sn2_i1'       , 0.0, $
    'sn2_g2'       , 0.0, $
    'sn2_r2'       , 0.0, $
    'sn2_i2'       , 0.0, $
    'dered_sn2_g1' , 0.0, $
    'dered_sn2_r1' , 0.0, $
    'dered_sn2_i1' , 0.0, $
    'dered_sn2_g2' , 0.0, $
    'dered_sn2_r2' , 0.0, $
    'dered_sn2_i2' , 0.0, $
    'goffstd'      , 0., $
    'grmsstd'      , 0., $
    'roffstd'      , 0., $
    'rrmsstd'      , 0., $
    'ioffstd'      , 0., $
    'irmsstd'      , 0., $
    'groffstd'     , 0., $
    'grrmsstd'     , 0., $
    'rioffstd'     , 0., $
    'rirmsstd'     , 0., $
    'goffgal'      , 0., $
    'grmsgal'      , 0., $
    'roffgal'      , 0., $
    'rrmsgal'      , 0., $
    'ioffgal'      , 0., $
    'irmsgal'      , 0., $
    'groffgal'     , 0., $
    'grrmsgal'     , 0., $
    'rioffgal'     , 0., $
    'rirmsgal'     , 0., $
    'nguide'       , 0L , $
    'seeing20'     , 0.0, $
    'seeing50'     , 0.0, $
    'seeing80'     , 0.0, $
    'rmsoff20'     , 0.0, $
    'rmsoff50'     , 0.0, $
    'rmsoff80'     , 0.0, $
    'airtemp'      , 0.0, $
    'xsigma'       , 0.0, $
    'xsigmin'      , 0.0, $
    'xsigmax'      , 0.0, $
    'wsigma'       , 0.0, $
    'wsigmin'      , 0.0, $
    'wsigmax'      , 0.0, $
    'xchi2'        , 0.0, $
    'xchi2min'     , 0.0, $
    'xchi2max'     , 0.0, $
    'skychi2'      , 0.0, $
    'schi2min'     , 0.0, $
    'schi2max'     , 0.0, $
    'fbadpix'      , 0.0, $
    'fbadpix1'     , 0.0, $
    'fbadpix2'     , 0.0, $
    'n_total'      , 0L,  $
    'n_galaxy'     , 0L,  $
    'n_qso'        , 0L,  $
    'n_star'       , 0L,  $
    'n_unknown'    , 0L,  $
    'n_sky'        , 0L, $
    'n_target_main',  0L, $
    'n_target_lrg1',  0L, $
    'n_target_lrg2',  0L, $
    'n_target_qso' ,  0L, $
    'success_main' , 0.0, $
    'success_lrg1' , 0.0, $
    'success_lrg2' , 0.0, $
    'success_qso'  , 0.0, $
    'status2d'     , 'Missing', $
    'statuscombine', 'Missing', $
    'status1d'     , 'Missing', $
    'public'       , ' ', $
    'qualcomments' , ' ' )
    
   if keyword_set(rawsn2) then begin
      sn2tag = 'platesn2'
      dereddened_sn2 = 0
   endif else begin
      sn2tag = 'deredsn2'
      dereddened_sn2 = 1
   endelse

   ;; Which columns to keep for HTML and ASCII files
   ;; For platelist
   trimtags1 = [ $
    ['plate'        ,   'i4'], $
    ['mjd'          ,   'i5'], $
    ['racen'        , 'f6.2'], $
    ['deccen'       , 'f6.2'], $
    ['run2d'        ,    'a'], $
    ['run1d'        ,    'a'], $
    ['platequality' ,    'a'], $
    [sn2tag         , 'f5.1'], $
    ['n_galaxy'     ,   'i3'], $
    ['n_qso'        ,   'i3'], $
    ['n_star'       ,   'i3'], $
    ['n_unknown'    ,   'i3'], $
    ['n_sky'        ,   'i3'], $
    ['survey'       ,    'a'], $
    ['programname'  ,    'a'], $
    ['chunk'        ,    'a'], $
    ['tileid'       ,    'i'], $
    ['public'       ,    'a']  ]

   ;; For platequality
   trimtags2 = [ $
    ['plate'        ,   'i4'], $
    ['mjd'          ,   'i5'], $
    ['run2d'        ,    'a'], $
    ['run1d'        ,    'a'], $
    ['sn2_g1'       , 'f5.1'], $
    ['sn2_i1'       , 'f5.1'], $
    ['sn2_g2'       , 'f5.1'], $
    ['sn2_i2' ,       'f5.1'], $
    ['fbadpix'      , 'f5.3'], $
;    ['success_main' , 'f5.1'], $
    ['success_lrg1' , 'f5.1'], $
    ['success_lrg2' , 'f5.1'], $
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
    root_dir=getenv('SPECLOG_DIR'), subdirectory='opfiles')
   publicdata = yanny_readone(publicfile, 'SPPLATELIST')
   if (NOT keyword_set(publicdata)) then $
    message, 'Missing spPlateList.par file'

   ;---------------------------------------------------------------------------
   ; Loop through all files
   ;---------------------------------------------------------------------------

   platefile = strarr(nfile)
   combparfile = strarr(nfile)
   comblogfile = strarr(nfile)
   combpsfile = strarr(nfile)
   zlogfile = ptrarr(nfile)
   zbestfile = ptrarr(nfile)
   zbestrun1d = ptrarr(nfile)

   for ifile=0L, nfile-1L do begin

      splog, 'Looking at ' + fullfile[ifile]

      junk = fileandpath(fullfile[ifile], path=path)  ; Determine PATH

      ;----------
      ; Test if INFILE specifies Yanny param files for spPlancomb.

      if (strmid(fullfile[ifile],strlen(fullfile[ifile])-4) EQ '.par') $
       then begin
         combparfile[ifile] = fullfile[ifile]
         yanny_read, fullfile[ifile], hdr=hdrp
         platefile[ifile] = $
          djs_filepath('spPlate-' $
           +string(yanny_par(hdrp,'plateid'),format='(i4.4)') $
           +'-'+string(yanny_par(hdrp,'MJD'),format='(i5.5)'), root_dir=path) $
           +'.fits'
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
      ;
      ; run1d is ALWAYS set.  If it is not input as a keyword, it has
      ; the value '*'
      ;
      if (keyword_set(run1d[0])) then begin
       allsubdirs = get_mjd_dir(path, mjd=run1d, /alldirs)
       zbestrun1d[ifile] = ptr_new( allsubdirs[where(~strcmp(allsubdirs,'runz_',5))] )
      endif else $
       zbestrun1d[ifile] = ptr_new('')
      nsubdir = n_elements(*zbestrun1d[ifile])
      zbestfile[ifile] = ptr_new(strarr(nsubdir))
      zlogfile[ifile] = ptr_new(strarr(nsubdir))
      for i=0L, nsubdir-1L do $
       (*zbestfile[ifile])[i] = djs_filepath('spZbest-'+platemjd+'.fits', $
        root_dir=path, subdir=((*zbestrun1d[ifile])[i]))
      for i=0L, nsubdir-1L do $
       (*zlogfile[ifile])[i] = djs_filepath('spDiag1d-'+platemjd+'.log', $
        root_dir=path, subdir=((*zbestrun1d[ifile]))[i])

      ; Read the combine plan file to get the list of all the 2D plan files
      ; from its Yanny header.
      ; Also get the mapping name from the combine par file in case we were
      ; unable to get it from the spPlate file.
      plist[ifile].mapname = (yanny_readone(combparfile[ifile], 'SPEXP', $
        hdr=hdrcomb))[0].mapname

      ;----------
      ; Find the state of the 2D reductions (not the combine step)

      statusdone = 0
      statusrun = 0
      statusmissing = 0

      ; Check status of individual 2D runs
      planlist = yanny_par(hdrcomb, 'planfile2d') ; Assume we find this
      logfile2d = '' ; List of 2D log files that exist
      for iplan=0L, n_elements(planlist)-1L do begin
         yanny_read, djs_filepath(planlist[iplan],root_dir=path), hdr=hdr2d
         plist[ifile].mjdlist += ' ' + strtrim(yanny_par(hdr2d, 'MJD'),2)
         thislogfile = repstr(planlist[iplan],'spPlan2d','spDiag2d')
         thislogfile = repstr(thislogfile,'.par','.log')
         thislogfile = (findfile(djs_filepath(thislogfile, root_dir=path)))[0]
         if (keyword_set(thislogfile)) then begin
            if (NOT keyword_set(logfile2d)) then logfile2d = thislogfile $
             else logfile2d = [logfile2d, thislogfile]
            ;;; spawn, 'tail -1 '+thislogfile, lastline
            lastline = get_lastline(thislogfile)
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

      ; Get RUN2D from the directory name...
      plist[ifile].run2d = fullrun2d[ifile]
      hdr1 = headfits(platefile[ifile], /silent, errmsg=errmsg)
      if (size(hdr1, /tname) EQ 'STRING') then begin
         plist[ifile].n_total = sxpar(hdr1, 'NAXIS2')
;         plist[ifile].plate = sxpar(hdr1, 'PLATEID')
;         plist[ifile].mjd = sxpar(hdr1, 'MJD')
         plist[ifile].mjdlist = sxpar(hdr1, 'MJDLIST')
;         thisrun2d = sxpar(hdr1, 'RUN2D', count=ct)
;         plist[ifile].run2d = (ct GT 0) ? thisrun2d : ''
;         plist[ifile].tileid = sxpar(hdr1, 'TILEID') ; Get from platePlans
;         plist[ifile].racen = sxpar(hdr1, 'RA') ; Get from platePlans
;         plist[ifile].deccen = sxpar(hdr1, 'DEC') ; Get from platePlans
         plist[ifile].cartid = sxpar(hdr1, 'CARTID')
         plist[ifile].tai = sxpar(hdr1, 'TAI')
         plist[ifile].tai_beg = sxpar(hdr1, 'TAI-BEG')
         plist[ifile].tai_end = sxpar(hdr1, 'TAI-END')
         plist[ifile].airmass = sxpar(hdr1, 'AIRMASS')
;         plist[ifile].airmass = tai2airmass(plist[ifile].ra, $
;          plist[ifile].dec, tai=plist[ifile].tai)
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
         
         ; If these keywords don't exist, these will just get 0
         plist[ifile].dered_sn2_g1 = sxpar(hdr1, 'SN2EXT1G')
         plist[ifile].dered_sn2_r1 = sxpar(hdr1, 'SN2EXT1R')
         plist[ifile].dered_sn2_i1 = sxpar(hdr1, 'SN2EXT1I')
         plist[ifile].dered_sn2_g2 = sxpar(hdr1, 'SN2EXT2G')
         plist[ifile].dered_sn2_r2 = sxpar(hdr1, 'SN2EXT2R')
         plist[ifile].dered_sn2_i2 = sxpar(hdr1, 'SN2EXT2I')

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
         plist[ifile].mapname = strtrim(sxpar(hdr1, 'NAME'))
         plist[ifile].statuscombine = 'Done'
      endif else begin
         ; Case where no spPlate file exists
         thislogfile = repstr(fileandpath(platefile[ifile]), $
          'spPlate','spDiagcomb')
         thislogfile = repstr(thislogfile,'.fits','.log')
         thislogfile = (findfile(djs_filepath(thislogfile, root_dir=path)))[0]

         if (keyword_set(thislogfile)) then begin
            ;;; spawn, 'tail -1 '+thislogfile, lastline
            lastline = get_lastline(thislogfile)
            if (strmatch(lastline[0], '*Successful completion*')) then begin
               ; Case where this combine log file completed,
               ; (but we're still missing the spPlate file, so must have failed)
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
      plist[ifile].survey = cinfo.survey
      plist[ifile].programname = cinfo.programname
      plist[ifile].chunk = cinfo.chunk
      plist[ifile].tileid = cinfo.tileid
      plist[ifile].racen = cinfo.racen
      plist[ifile].deccen = cinfo.deccen
      plist[ifile].epoch = cinfo.epoch

      ;----------
      ; Determine which public data release has this plate+MJD

      j = where(plist[ifile].plate EQ publicdata.plate $
       AND plist[ifile].mjd EQ publicdata.mjd)
      if (j[0] NE -1) then begin
         copy_struct_inx, publicdata[j[0]], plist, index_to=ifile
      endif

   endfor

   ;----------
   ; Remove from the plate list earlier reductions of the same plugging
   ; (keeping only the most recent MJD of each plugging).

   qkeep = bytarr(nfile)

   ; First get the unique list of MAPNAME+RUN2D, then mark the most recent MJD
   ; of each as the good one.
   isort = sort(plist.mapname+plist.run2d)
   isort = isort[ uniq(plist[isort].mapname+plist[isort].run2d) ]
   maplist = plist[isort].mapname+plist[isort].run2d
   for imap=0L, n_elements(maplist)-1L do begin
      indx = where(plist.mapname+plist.run2d EQ maplist[imap])
      junk = max(plist[indx].mjd, imax)
      qkeep[indx[imax]] = 1
   endfor

   ; Don't discard any where MAPNAME isn't set
   indx = where(strtrim(plist.mapname) EQ '')
   if (indx[0] NE -1) then qkeep[indx] = 1

   ; List partially-combined plates that we're discarding from the list
   for ifile=0L, nfile-1L do begin
      if (qkeep[ifile] NE 1) then begin
         splog, 'Discard partially-combined ' + combparfile[ifile]
         if (keyword_set(killpartial)) then begin
            killfiles = [ combparfile[ifile], $
                          comblogfile[ifile], $
                          combpsfile[ifile], $
                          platefile[ifile], $
                          *zbestfile[ifile], $
                          *zlogfile[ifile] ]
            for ikill=0L, n_elements(killfiles)-1L do $
             splog, 'KILLPARTIAL ', killfiles[ikill]
            rmfile, killfiles
         endif
      
      ; Drop plates from carts in skipcart
      endif else begin
         if keyword_set(skipcart) then begin
            tmp = where(skipcart EQ plist[ifile].cartid, N)
            if N GT 0 then begin
               cartid = STRTRIM(STRING(plist[ifile].cartid), 1)
               splog, 'Dropping ' + combparfile[ifile] + ' (cart ' + cartid + ')'
               qkeep[ifile] = 0
            endif
         endif
      endelse
      
   endfor

   ;-----
   ; Trim the plate list, and update NFILE to this trimmed number

   ikeep = where(qkeep, nfile)
   idelete = where(qkeep EQ 0, ndelete)
   ; Free pointers for partial reductions
   if (ndelete GT 0) then begin
      for i=0L, ndelete-1L do begin
         ptr_free, zlogfile[idelete[i]]
         ptr_free, zbestfile[idelete[i]]
         ptr_free, zbestrun1d[idelete[i]]
      endfor
   endif

   plist = plist[ikeep]
   platefile = platefile[ikeep]
   combparfile = combparfile[ikeep]
   comblogfile = comblogfile[ikeep]
   combpsfile = combpsfile[ikeep]
   zlogfile = zlogfile[ikeep]
   zbestfile = zbestfile[ikeep]
   zbestrun1d = zbestrun1d[ikeep]

   ;----------
   ; Make a list of one S/N for each plate which is the minimum of
   ; G1, I1, G2, I2.
   ; Assign a plate quality, but do not over-write any plate quality
   ; from the manually-assigned one in the "spPlateList.par" file.

   qualstring = ['bad', 'marginal', 'good']
   for ifile=0L, nfile-1L do begin
      if (strtrim(plist[ifile].statuscombine,2) EQ 'Done') then begin
         nexp_min = min( $
          [plist[ifile].nexp_b1, plist[ifile].nexp_b2, $
          plist[ifile].nexp_r1, plist[ifile].nexp_r2], max=nexp_max)
         plist[ifile].platesn2 = min( $
          [plist[ifile].sn2_g1, plist[ifile].sn2_i1, $
           plist[ifile].sn2_g2, plist[ifile].sn2_i2])
         plist[ifile].deredsn2 = min( $
          [plist[ifile].dered_sn2_g1, plist[ifile].dered_sn2_i1, $
           plist[ifile].dered_sn2_g2, plist[ifile].dered_sn2_i2])
         if keyword_set(rawsn2) then begin
            min_sn2_b = min([plist[ifile].sn2_g1, plist[ifile].sn2_g2])
            min_sn2_r = min([plist[ifile].sn2_i1, plist[ifile].sn2_i2])
         endif else begin
            min_sn2_b = min([plist[ifile].dered_sn2_g1, $
                             plist[ifile].dered_sn2_g2])
            min_sn2_r = min([plist[ifile].dered_sn2_i1, $
                             plist[ifile].dered_sn2_i2])
         endelse
         iqual = 2
         if (min_sn2_b LT 10.0) OR (min_sn2_r LT 22.0) then $
            iqual = iqual < 0
;         if (plist[ifile].platesn2 LT 13) then iqual = iqual < 0
;         if (plist[ifile].platesn2 LT 15) then iqual = iqual < 1
         if (plist[ifile].fbadpix GT 0.10) then iqual = iqual < 0
;         if (plist[ifile].fbadpix GT 0.05) then iqual = iqual < 1
         ; For reductions before v5_1, NEXP_MIN and NEXP_MAX are always zero
         if (nexp_max GT 0) then begin
            if (nexp_min LT 3) then iqual = iqual < 0
;            if (nexp_min LT 3) then iqual = iqual < 1
         endif
         if (NOT keyword_set(strtrim(plist[ifile].platequality))) then $
          plist[ifile].platequality = qualstring[iqual]
      endif
   endfor

   ;----------
   ; Decide which plates constitute unique tiles with the required S/N,
   ; then set QSURVEY=1.
   ; Also insist that PROGNAME='main'.

   ; First get the unique list of TILE
   isort = sort(plist.tileid)
   isort = isort[ uniq(plist[isort].tileid) ]
   tilelist = plist[isort].tileid

   for itile=0L, n_elements(tilelist)-1L do begin
      indx = where(plist.tileid EQ tilelist[itile] $
       AND (strtrim(plist.platequality,2) EQ 'good' $
         OR strtrim(plist.platequality,2) EQ 'marginal') $
       AND strtrim(plist.survey,2) EQ 'boss', ct)
      if (ct GT 0) then begin
         snbest = max(plist[indx].platesn2, ibest)
         plist[indx[ibest]].qsurvey = 1
      endif
   endfor

   ;---------------------------------------------------------------------------
   ; Read the Spectro-1D files

   for i=0L, nfile-1L do begin

      for j=0L, n_elements((*zbestfile[i]))-1L do begin

         if (j EQ 0) then begin
            ifile = i ; First instance of z-file for this plate
         endif else begin
            ifile = n_elements(plist) ; Append another entry to plate list
            plist = [plist, plist[i]]
         endelse

         plist[ifile].run1d = (*zbestrun1d[i])[j]

         ;----------
         ; Read Zbest file - get status of 1D

         if (keyword_set((*zbestfile[i])[j])) then $
          hdr2 = headfits((*zbestfile[i])[j], /silent, errmsg=errmsg) $
         else $
          hdr2 = 0
         if (size(hdr2, /tname) EQ 'STRING') then begin
            zans = mrdfits((*zbestfile[i])[j], 1, /silent)
            plug = mrdfits(platefile[i], 5, /silent)
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
            plist[ifile].status1d = 'Done'

            nobj = n_elements(zans)
            targets = strarr(nobj)

            for iobj=0L, nobj-1L do $
             targets[iobj] = sdss_flagname('TARGET',plug[iobj].primtarget, $
              /silent, /concat)+' '
            if (tag_exist(plug,'BOSS_TARGET1')) then $
             for iobj=0L, nobj-1L do $
              targets[iobj] = sdss_flagname('BOSS_TARGET1', $
               plug[iobj].boss_target1, /silent, /concat)+' '

            ;; Objects which shouldn't count against the success statistics
            bad_fiber = sdss_flagval('ZWARNING', 'LITTLE_COVERAGE') OR $
                        sdss_flagval('ZWARNING', 'UNPLUGGED') OR $
                        sdss_flagval('ZWARNING', 'BAD_TARGET')

            imain = where((strmatch(targets,'*GALAXY *') $
             OR strmatch(targets,'*GALAXY_BIG *') $
             OR strmatch(targets,'*GALAXY_BRIGHT_CORE *')) $
             AND (zans.zwarning_noqso AND bad_fiber) EQ 0, nmain)
            ilrg1 = where((strmatch(targets,'*GALAXY_RED *') $
             OR strmatch(targets,'*GALAXY_RED_II *') $
             OR strmatch(targets,'*GAL_LOZ *')) $
             AND (zans.zwarning_noqso AND bad_fiber) EQ 0,  nlrg1)
            ilrg2 = where((strmatch(targets,'*GAL_HIZ *') $
             OR strmatch(targets,'*GAL_CMASS*')) $
             AND (zans.zwarning_noqso AND bad_fiber) EQ 0, nlrg2)
;            iqso = where(strmatch(targets,'*QSO_HIZ *') $
;             OR strmatch(targets,'*QSO_CAP *') $
;             OR strmatch(targets,'*QSO_SKIRT *') $
;             OR strmatch(targets,'*QSO_FIRST_CAP *') $
;             OR strmatch(targets,'*QSO_FIRST_SKIRT *'), nqso)
            iqso = where(strmatch(targets,'*QSO_*') $
                AND (zans.zwarning AND bad_fiber) EQ 0, nqso)

            plist[ifile].n_target_main = nmain
            plist[ifile].n_target_lrg1 = nlrg1
            plist[ifile].n_target_lrg2 = nlrg2
            plist[ifile].n_target_qso = nqso
            if (nmain GT 0) then $
             plist[ifile].success_main = $
              100 * total(zans[imain].zwarning EQ 0 $
              AND (strmatch(zans[imain].class,'GALAXY*') $
               OR strmatch(zans[imain].class,'QSO*'))) / nmain
            if (nlrg1 GT 0) then $
             plist[ifile].success_lrg1 = $
              100 * total(zans[ilrg1].zwarning EQ 0 $
              AND strmatch(zans[ilrg1].class,'GALAXY*')) / nlrg1
            if (nlrg2 GT 0) then $
             plist[ifile].success_lrg2 = $
              100 * total(zans[ilrg2].zwarning EQ 0 $
              AND strmatch(zans[ilrg2].class,'GALAXY*')) / nlrg2
            if (nqso GT 0) then $
             plist[ifile].success_qso = $
              100 * total(zans[iqso].zwarning EQ 0 $
              AND strmatch(zans[iqso].class,'QSO*')) / nqso
         endif else begin
            ;----------
            ; Find the state of the 1D reductions -- spZbest file is missing

            print, run1d, ' ', (*zlogfile[i])[j]

            if (file_search((*zlogfile[i])[j])) then begin
               ;;; spawn, 'tail -1 '+(*zlogfile[i])[j], lastline
               lastline = get_lastline((*zlogfile[i])[j])
               if (strmatch(lastline[0], '*Successful completion*')) then begin
                  ; Case where this 1D log file completed, which is not
                  ; a case that should ever occur.
                  plist[ifile].status1d = 'FAILED'; Should have found spZbest
               endif else begin
                  ; Case where this 1D log file isn't completed
                  if (keyword_set(purge1d)) then begin
                     splog, 'PURGE1D ', (*zlogfile[i])[j]
                     rmfile, (*zlogfile[i])[j]
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
   endfor

   ;---------------------------------------------------------------------------
   ; Write ASCII + HTML output files

   alias = [['PROGRAMNAME'  , 'PROG'    ], $
            ['PLATESN2'     , 'SN^2'    ], $
            ['DEREDSN2'     , 'SN^2'    ], $
            ['N_GALAXY'     , 'N_gal'   ], $
            ['N_QSO'        , 'N_QSO'   ], $
            ['N_STAR'       , 'N_star'  ], $
            ['N_UNKNOWN'    , 'N_unk'   ], $
            ['N_SKY'        , 'N_sky'   ], $
            ['FBADPIX'      , 'Badpix'  ], $
;            ['SUCCESS_MAIN' , '%Main'   ], $
            ['SUCCESS_LRG1' , '%LRG1'   ], $
            ['SUCCESS_LRG2' , '%LRG2'   ], $
            ['SUCCESS_QSO'  , '%QSO'    ], $
            ['STATUS2D'     , '2D'      ], $
            ['STATUSCOMBINE', 'Combine' ], $
            ['STATUS1D'     , '1D'      ], $
            ['DERED_SN2_G1' , 'SN2_G1'  ], $
            ['DERED_SN2_G2' , 'SN2_G2'  ], $
            ['DERED_SN2_I1' , 'SN2_G1'  ], $
            ['DERED_SN2_I2' , 'SN2_I2'  ], $
            ['PLATEQUALITY' , 'QUALITY' ] ]

   isort1 = reverse(sort(strtrim(strcompress(plist.run2d+' ' $
    +string(99999-plist.plate)),2)))
   isort2 = reverse(sort(strtrim(strcompress(plist.run2d+' ' $
    +string(plist.mjd)),2)))

   toptext = [ $
    '<p>Last Update: '+ systime()+'</p>', '<ul>', $
    '<li><a href="http://data.sdss3.org/sas/bossredux/">HOME</a></li>', $
    '<li>Plate list sorted by <a href="platelist.html">plate</a>,' $
     + ' <a href="platelist-mjdsort.html">MJD</a></li>', $
    '<li>Plate quality sorted by <a href="platequality.html">plate</a>,' $
     + ' <a href="platequality-mjdsort.html">MJD</a></li>', $
    '<li>Plate list as <a href="platelist.fits">FITS</a></li>','</ul>']

   if dereddened_sn2 NE 0 then begin
      toptext = [toptext, '<p>(S/N)^2 values are corrected for galactic dust reddening</p>']
   endif else begin
      toptext = [toptext, '<p>(S/N)^2 values are <b>not</b> corrected for galactic dust reddening</p>']
   endelse

   platelist_write, plist[isort1], trimtags=trimtags1, alias=alias, $
    fileprefix='platelist', toptext=toptext, outdir=outdir, $
    title='SDSS Spectroscopy Plates Observed List'
   platelist_write, plist[isort2], trimtags=trimtags1, alias=alias, $
    fileprefix='platelist-mjdsort', toptext=toptext, outdir=outdir, $
    title='SDSS Spectroscopy Plates Observed List'

   platelist_write, plist[isort1], trimtags=trimtags2, alias=alias, $
    fileprefix='platequality', toptext=toptext, outdir=outdir, $
    title='SDSS Spectroscopy Plate Quality List'
   platelist_write, plist[isort2], trimtags=trimtags2, alias=alias, $
    fileprefix='platequality-mjdsort', toptext=toptext, outdir=outdir, $
    title='SDSS Spectroscopy Plate Quality List'

   ;----------
   ; Write the FITS binary table

   mwrfits, plist, fitsfile, /create

   for i=0L, n_elements(zlogfile)-1L do ptr_free, zlogfile[i]
   for i=0L, n_elements(zbestfile)-1L do ptr_free, zbestfile[i]
   for i=0L, n_elements(zbestrun1d)-1L do ptr_free, zbestrun1d[i]

   return
end
;------------------------------------------------------------------------------
