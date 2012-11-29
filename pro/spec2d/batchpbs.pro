;+
; NAME:
;   batchpbs
;
; PURPOSE:
;   Batch process Spectro-2D and Spectro-1D reductions based upon
;   already-built plan files.
;
; CALLING SEQUENCE:
;   batchpbs, [ platenums, topdir=, run2d=, run1d=, platestart=, plateend=, $
;    mjd=, mjstart=, mjend=, upsvers2d=, upsvers1d=, /zcode, queue=, /skip2d, $
;     /clobber, /nosubmit ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   platenums  - Plate numbers to reduce; default to '*'
;   topdir     - Optional override value for the environment variable
;                $BOSS_SPECTRO_REDUX.  Used as input directory for generating
;                batch scripts and written into the scripts themselves.
;   run2d      - Optional override value for the environment variable $RUN2D
;   run1d      - Optional override value for the environment variable $RUN1D
;   platestart - Starting plate number.
;   plateend   - Ending plate number.
;   mjd        - MJD dates to reduce; default to all.
;                Select based upon the MJD of the combine plan file, and
;                reduce data from all nights needed for that combined plate+MJD.
;   mjstart    - Starting MJD dates to reduce.
;   mjend      - Ending MJD dates to reduce.
;   upsvers2d  - If set, then do a "setup idlspec2d $UPSVERS2D" on the
;                remote machine before executing Spectro-2D.  This allows
;                you to batch jobs using a version other than that which
;                is declared current under UPS.
;   upsvers1d  - If set, then do a "setup idlspec2d $UPSVERS2D" on the
;                remote machine before executing Spectro-1D.  This allows
;                you to batch jobs using a version other than that which
;                is declared current under UPS.
;   zcode      - If set, run Zcode in auto mode.
;   queue      - If set, sets the submit queue.
;   skip2d     - If set, then skip the Spectro-2D reductions.
;   skip1d     - If set, then skip the Spectro-1D reductions.
;   clobber    - If set, then reduce all specified plates.  The default is
;                to not reduce plates where the script file already exists.
;   nosubmit   - If set, generate script file but don't submit to queue
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   This is currently written to batch all Spectro-2D and Spectro-1D
;   reductions, without the option of doing only one or the other.
;
;   Script files are generated for partial reductions of a plate,
;   although those jobs are not submitted.  It would be less confusing
;   to not make those script files.
;
; REVISION HISTORY:
;   17-Jan-2010  Written by D. Schlegel, LBL
;   03-Oct-2012  Added /skip1d option, Stephen Bailey LBL
;-
;------------------------------------------------------------------------------
pro batchpbs, platenums1, topdir=topdir1, run2d=run2d1, run1d=run1d1, $
 platestart=platestart, plateend=plateend, $
 mjd=mjd, mjstart=mjstart, mjend=mjend, $
 upsvers2d=upsvers2d, upsvers1d=upsvers1d, $
 rawdata_dir=rawdata_dir, $
 boss_spectro_redux=boss_spectro_redux, $
 zcode=zcode, $
 queue=queue, skip2d=skip2d, skip1d=skip1d, clobber=clobber, nosubmit=nosubmit

   if (size(platenums1,/tname) EQ 'STRING') then platenums = platenums1 $
    else if (keyword_set(platenums1)) then $
      platenums = string(platenums1,format='(i4.4)') $
    else platenums = '*'

   ;----------
   ; Determine the top-level of the output directory tree

   if (keyword_set(topdir1)) then topdir = topdir1 $
    else topdir = getenv('BOSS_SPECTRO_REDUX')
   splog, 'Setting TOPDIR=', topdir
   if (keyword_set(run2d1)) then run2d = strtrim(run2d1,2) $
    else run2d = getenv('RUN2D')
   splog, 'Setting RUN2D=', run2d
   if (keyword_set(run1d1)) then run1d = strtrim(run1d1,2) $
    else run1d = getenv('RUN1D')
   splog, 'Setting RUN1D=', run1d

   topdir2d = djs_filepath('', root_dir=topdir, subdir=run2d)
   if (keyword_set(run1d)) then run1dstr = ',run1d="'+run1d+'"' $
    else run1dstr = ''
   if (keyword_set(run2d)) then run2dstr = ',run2d="'+run2d+'"' $
    else run2dstr = ''

   ;----------
   ; Create list of plate directories
   ; Limit the list to only those specified by PLATENUMS,PLATESTART,PLATEEND

   platedirs = get_mjd_dir(topdir2d, mjd=platenums, mjstart=platestart, $
    mjend=plateend)

   if (NOT keyword_set(platedirs[0])) then begin
      splog, 'No directories found'
      return
   endif
   ndir = n_elements(platedirs)

   ;----------
   ; In each plate directory, find all 'spPlancomb*.par' files

   for idir=0L, ndir-1L do begin
      planfile = findfile( $
       djs_filepath('spPlancomb*.par', root_dir=topdir2d, $
        subdir=platedirs[idir]), count=nfile)

      for ifile=0, nfile-1 do begin
         yanny_read, planfile[ifile], hdr=hdr
         thismjd = long(yanny_par(hdr, 'MJD'))

         ; Decide if THISMJD is within the bounds specified by MJD,MJSTART,MJEND
         if (mjd_match(thismjd, mjd=mjd, mjstart=mjstart, mjend=mjend)) $
          then begin
            if (keyword_set(platelist)) then begin
               platelist = [platelist, platedirs[idir]]
               planlist = [planlist, planfile[ifile]]
            endif else begin
               platelist = platedirs[idir]
               planlist = planfile[ifile]
            endelse
         endif
      endfor
   endfor

   nplate = n_elements(planlist)
   if (nplate EQ 0) then begin
      splog, 'No plan files found'
      return
   endif

; Do not use spPlancomb files that only have a subset of the MJDs !!!???

   ;----------
   ; For each combine plan file, generate the IDL script files

   fullscriptfile = strarr(nplate)
   plateid = lonarr(nplate)
   mjd_beg = lonarr(nplate)
   mjd_end = lonarr(nplate)

   fq = "'"
   qbatch = bytarr(nplate) + 1B ; default to reduce all plates
   for iplate=0, nplate-1 do begin
      ; Find all relevant 2D plan files
      yanny_read, planlist[iplate], hdr=hdr
      planfile2d = yanny_par(hdr, 'planfile2d')
      plateid[iplate] = yanny_par(hdr, 'plateid')
      mjd = yanny_par(hdr, 'MJD')
      platemjd = string(plateid[iplate],format='(i4.4)') + '-' $
       + string(mjd,format='(i5.5)')
      platefile = 'spPlate-'+platemjd+'.fits'

      ; Track the beginning and ending MJD going into this plate
      mjd_beg[iplate] = min(strmid(planfile2d,14,5))
      mjd_end[iplate] = max(strmid(planfile2d,14,5))

      ; Split the combine plan file name into a directory and file name
      planfilecomb = fileandpath(planlist[iplate], path=pathcomb)

      ; Construct the name of the batch file
      ; "script" -> "redux" to make names fit within qstat column widths
      ; fullscriptfile[iplate] = djs_filepath('script-'+platemjd, $
      fullscriptfile[iplate] = djs_filepath('redux-'+platemjd, $
       root_dir=pathcomb)
      if (keyword_set(skip2d)) then fullscriptfile[iplate] += '-' + run1d

      ; Write the batch file
      if (keyword_set(clobber) EQ 0) then $
       qbatch[iplate] = file_test(fullscriptfile[iplate]) EQ 0

      if (qbatch[iplate]) then begin
         openw, olun, fullscriptfile[iplate], /get_lun
         printf, olun, '# Auto-generated batch file '+systime()
         printf, olun, '#PBS -l nodes=1'
         printf, olun, '#PBS -l walltime=48:00:00'
         printf, olun, '#PBS -W umask=0022'
         printf, olun, '#PBS -V'
         printf, olun, '#PBS -j oe'
         ; set queue if asked
         if (keyword_set(queue)) then $
             printf, olun, '#PBS -q ' + queue
         printf, olun, 'cd $PBS_O_WORKDIR'

         ; Override environment variables if requested
         if (keyword_set(rawdata_dir)) then begin
             printf, olun, 'export BOSS_SPECTRO_DATA='+rawdata_dir
         endif
         if (keyword_set(boss_spectro_redux)) then begin
             printf, olun, 'export BOSS_SPECTRO_REDUX='+boss_spectro_redux
         endif

         printf, olun, ''
         printf, olun, '#- Echo commands to make debugging easier'
         printf, olun, 'set -o verbose'

         ; printf, olun, ''
         ; printf, olun, '#- Dump job environment for debugging'
         ; envlog = 'env-'+platemjd+'.txt'
         ; printf, olun, 'printenv > '+envlog
         
         if topdir ne getenv('BOSS_SPECTRO_REDUX') then begin
            printf, olun, ''
            printf, olun, '#- Override default BOSS_SPECTRO_REDUX'
            printf, olun, 'export BOSS_SPECTRO_REDUX='+topdir
         endif
         
         printf, olun, ''
         printf, olun, '#- The real work'
         if (keyword_set(skip2d) EQ 0) then begin
            ; Set up requested code version
            if (keyword_set(upsvers2d)) then $
             printf, olun, 'setup idlspec2d '+upsvers2d

            ; Create sorted photoPlate files
            for i=0, n_elements(planfile2d)-1 do $
             printf, olun, 'echo '+fq+'sdss_plate_sort,"'+planfile2d[i]+'"'+fq+' | idl'

            ; Run Spectro-2D
            for i=0, n_elements(planfile2d)-1 do $
             printf, olun, 'echo '+fq+'spreduce2d,"'+planfile2d[i]+'"'+fq+' | idl'
            printf, olun, 'echo '+fq+'spcombine_v5,"'+planfilecomb+'"'+fq+' | idl'
         endif

         ; Run Spectro-1D
         if (keyword_set(skip1d) EQ 0) then begin
            if (keyword_set(upsvers1d)) then $
             printf, olun, 'setup idlspec2d '+upsvers1d
            printf, olun, 'echo '+fq+'spreduce1d,"'+platefile+'"'+run1dstr+fq+' | idl'
         endif

         ; Run Zcode
         if (keyword_set(zcode)) then begin
            printf, olun, ''
            printf, olun, 'setup runz'
            printf, olun, 'runz_BOSS.sh ' + platefile +' -a'
            printf, olun, 'runz_BOSS.sh ' + platefile +' -a -G -t GAL'
         endif
         
         ; splog, "run1d is ", run1d
         ; splog, "run2d is ", run2d
         
         ; Make pretty pictures
         ;- post-DR9, no longer supported; use spectrawebapp or plotspec instead
         ; idlcmd  = "plate_spec_image, " + string(plateid[iplate],format='(i4.4)') 
         ; idlcmd += ", mjd=" + string(mjd,format='(i5.5)')
         ; idlcmd += ", run1d='" + run1d + "'"
         ; idlcmd += ", run2d='" + run2d + "'"
         ; idlcmd += ", /silent"
         ; printf, olun, ''
         ; printf, olun, '#- Make pretty pictures'
         ; printf, olun, 'idl -e "' + idlcmd + '"'
         
         close, olun
         free_lun, olun
      endif
   endfor

   ;----------
   ; Do not reduce any plan files that are only partial reductions
   ; (If we tried, then we would get multiple instances of SPREDUCE2D
   ; running on the same data.)
   ; Make this decision regardless of the values of CLOBBER.

   for i=0, nplate-1 do begin
      indx = where(plateid EQ plateid[i] AND mjd_end GT mjd_end[i] $
       AND mjd_beg LE mjd_end[i], ct)
      if (ct GT 0) then qbatch[i] = 0B
   endfor

   ;----------
   ; Trim the plate list to only those needing reductions

   ibatch = where(qbatch, nbatch)
   if (nbatch EQ 0) then begin
      splog, 'All plates have been reduced'
      return
   endif

   ;----------
   ; Submit jobs to the PBS queue

   for i=0L, nbatch-1L do begin
      thisfile = fileandpath(fullscriptfile[ibatch[i]], path=thispath)
      if (keyword_set(thispath)) then cd, thispath, current=origdir
      if keyword_set(nosubmit) then begin
          splog, 'Generated '+thisfile+' but not submitting to queue'
      endif else begin
          splog, 'Submitting '+thisfile
          spawn, 'qsub '+thisfile
      endelse
   endfor
   cd, origdir

   return
end
;------------------------------------------------------------------------------
