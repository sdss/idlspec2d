;+
; NAME:
;   remove2redo
;
; PURPOSE:
;   Removed specified exposures from the SOS log files and re-reduce them.
;
; CALLING SEQUENCE:
;   remove2redo, [ mjd=, plate=, expnum=, /not_sos ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - MJD number; if not set, then select the most recent MJD
;                in the $SPECTROLOG_DIR directory.
;   plate      - If specified, then search for all sdR files for this plate
;                that were not successfully reduced.  They are assumed to
;                have been successfully reduced if they have an entry in
;                HDU #1, 2, or 3 in the SOS log file.
;   expnum     - If specified, then force the re-reduction of this exposure
;                number(s) even if they have already been successfully reduced.
;                Do not specify both PLATE and EXPNUM.
;   not_sos    - This keyword can be set to run this proc on a machine
;                that is not named "sos".  This would only be done for
;                testing purposes, or if Son-of-Spectro has been moved
;                to another machine.
;
; OUTPUT:
;
; COMMENTS:
;   The way this procedure works is to actually remove the sdR*.fit files
;   from /data/spectro.  The Son-of-Spectro cron jobs will then re-copy
;   and re-reduce those files.  It does not look for or remove sdR*.fit.gz
;   files, which SOS currently makes too.  If the raw data no longer exists
;   on sdsshost.apo, then the data will not be re-reduced.
;
;   If $SPECTROLOG_DIR is not set, then it is assumed to be
;     /data/spectro/spectrologs
;   When using the PLATE keyword, never re-reduce the very last exposure
;   number on disk.  If the last exposure is not for the specified plate,
;   then we can try re-reducing any exposures for that plate.
;
; EXAMPLES:
;
; BUGS:
;   Currently, either EXPNUM or PLATE must be specified.  Probably nothing
;   happens if neither is set.
;
; PROCEDURES CALLED:
;   apo_log2html
;   djs_filename()
;   djs_lockfile()
;   djs_modfits
;   djs_unlockfile
;   sdsshead()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   03-April-2001  Written by Scott Burles, FNAL
;-
;------------------------------------------------------------------------------
pro remove2redo, mjd=mjd, plate=plate, expnum=expnum, not_sos=not_sos

   quiet = !quiet
   !quiet = 1

   if (keyword_set(plate) AND keyword_set(expnum)) then begin
      splog, 'Cannot specify both MJD and PLATE -- quitting!'
      !quiet = quiet
      return
   endif

   if (NOT keyword_set(getenv('SPECLOG_DIR'))) then $
    setenv,'SPECLOG_DIR=/data/spectro/astrolog'

   ;----------
   ; Insist that this proc only run on machines named "sos"

   if (NOT keyword_set(apo_uname)) then begin
      spawn, 'uname -n', uname_string
      apo_uname = (strsplit(uname_string[0], '.', /extract))[0]
   endif

   if (apo_uname NE 'sos' AND keyword_set(not_sos) EQ 0) then begin
      print, 'This procedure can only be run on the machine sos.apo.nmsu.edu'
      print, 'If Son-of-Spectro is now running on a different machine, you'
      print, 'can call this routine with /not_sos'
      return
   endif

   ;----------
   ; If MJD is not specified, then find the most recent MJD for output files

   spectrolog_dir = getenv('SPECTROLOG_DIR')
   if (NOT keyword_set(spectrolog_dir)) then $
    spectrolog_dir = '/data/spectro/spectrologs'

   if (NOT keyword_set(mjd)) then begin
      mjdlist = get_mjd_dir(spectrolog_dir, mjstart=1, mjend=99999, mjd='?????')
      mjd = (reverse(mjdlist[sort(mjdlist)]))[0]
   endif
   splog, 'Selecting MJD=', mjd

   ;----------
   ; Find the log file

   mjdstr = string(mjd, format='(i5.5)')
   mjddir = concat_dir(spectrolog_dir, mjdstr)

   logfile = 'logfile-' + mjdstr + '.fits'
   logfile = filepath(logfile, root_dir=mjddir)
   if (NOT keyword_set(findfile(logfile))) then begin
      splog, 'Unable to find logfile '+logfile
      !quiet = quiet
      return
   endif

   if (keyword_set(plate) OR keyword_set(expnum)) then begin
      rawdata_dir = getenv('BOSS_SPECTRO_DATA')
      if (NOT keyword_set(rawdata_dir)) then $
       rawdata_dir = '/data/spectro'
   endif

   ;----------
   ; If PLATE is specified, then read the headers of all the files
   ; to find out which files correspond to this plate.
   ; All those files are put in the array FULLNAME.

   if (keyword_set(plate)) then begin
      filename = 'sdR-*-*.fit' ; Only look at ".fit" files, not ".fit.gz"
      fullname = findfile( djs_filepath(filename, root_dir=rawdata_dir, $
       subdirectory=mjdstr), count=nfile)
      if (nfile EQ 0) then begin
         splog, 'No sdR files for MJD=', mjd
         !quiet = quiet
         return
      endif
      qthisplate = bytarr(nfile)
      allexpnum = lonarr(nfile)
      splog, 'Reading headers for all ', nfile, ' files on MJD=', mjd
      for i=0, nfile-1 do begin
         ; Print something since this might take a while to read all the
         ; FITS headers...
         print, format='(".",$)'

         hdr = sdsshead(fullname[i], /do_lock)
         allexpnum[i] = sxpar(hdr, 'EXPOSURE')
         if (sxpar(hdr, 'PLATEID') EQ plate) then qthisplate[i] = 1B
      endfor

      ifile = where(qthisplate, nfile)
      if (nfile EQ 0) then begin
         splog, 'No sdR files for MJD=', mjd, ' PLATE=', plate
         !quiet = quiet
         return
      endif
      fullname = fullname[ifile]

      shortname = fileandpath(fullname)
      qdone = bytarr(nfile) ; Start as all 0's

      ; Always mark the last exposure number as 'done', since SOS may
      ; still be trying to reduce it.
      idone = where(allexpnum[ifile] EQ max(allexpnum))
      if (idone[0] NE -1) then $
       qdone[idone] = 1B
   endif

   ;----------
   ; Lock the file to do this - otherwise we might read/write to a partially
   ; written file.

   splog, 'Trying to lock the logfile  ' + logfile
   while(djs_lockfile(logfile) EQ 0) do wait, 1

   ;----------
   ; Loop through each HDU in the log file.
   ; If EXPNUM is specified, then remove all entries corresponding to this
   ; exposure number.
   ; If PLATE is specified, then determine which of our list of files for
   ; that plate have been succesfully reduced, and we will try re-reducing
   ; the others.

   splog, 'Reading the logfile ' + logfile
   for thishdu=1, 5 do begin
      rstruct = mrdfits(logfile, thishdu, /silent)
      nstruct = n_elements(rstruct)

      if (keyword_set(expnum) AND keyword_set(rstruct)) then begin
         qkeep = bytarr(nstruct)
         for i=0, nstruct-1 do $
          qkeep[i] = total(rstruct[i].expnum EQ expnum) EQ 0
         ikeep = where(qkeep, nkeep)
         if (nkeep LT nstruct) then begin
            if (nkeep EQ 0) then $
             djs_modfits, logfile, 0, exten_no=thishdu, /delete_data $
            else $
             djs_modfits, logfile, rstruct[ikeep], exten_no=thishdu
         endif
      endif

      if (keyword_set(plate) AND keyword_set(rstruct)) then begin
         if (thishdu LE 4) then begin
            ; Bias, flat, arc, or science/smear -- note any file that's been
            ; successfully reduced.
            for i=0, nstruct-1 do $
             qdone = qdone OR (rstruct[i].filename EQ shortname)
         endif else begin
            ; Warning or aborts -- remove any for files that haven't been
            ; successfully reduced and that we'll re-reduce.
            qkeep = bytarr(nstruct)
            for i=0, nstruct-1 do $
             qkeep[i] = total(rstruct[i].filename EQ shortname $
              AND qdone EQ 0) EQ 0
            ikeep = where(qkeep, nkeep)
            if (nkeep LT nstruct) then begin
               splog, 'Removing ', nstruct-nkeep, ' warning and abort messages'
               if (nkeep EQ 0) then $
                djs_modfits, logfile, 0, exten_no=thishdu, /delete_data $
               else $
                djs_modfits, logfile, rstruct[ikeep], exten_no=thishdu
            endif
         endelse
      endif
   endfor

   ;----------
   ; Now unlock the log file. 

   djs_unlockfile, logfile

   ;----------
   ; Construct the list of files to re-reduce.

   if (keyword_set(expnum)) then begin
      camnames = ['b1', 'b2', 'r1', 'r2']
      redofiles = ''
      for iexp=0, n_elements(expnum)-1 do $
       for icam=0, n_elements(camnames)-1 do $
        redofiles = [redofiles, $
         findfile(filepath('sdR-' + camnames[icam] + '-' $
         + string(expnum[iexp],format='(i8.8)') + '.fit', $
         root_dir=rawdata_dir, subdirectory=mjdstr))]
      iredo = where(redofiles NE '')
      if (iredo[0] NE -1) then redofiles = redofiles[iredo] $
       else redofiles = 0
   endif

   if (keyword_set(plate)) then begin
      iredo = where(qdone EQ 0)
      if (iredo[0] NE -1) then redofiles = fullname[iredo]
   endif

   for iredo=0, n_elements(redofiles)-1 do begin
      splog, 'Marking file to re-reduce: ' + redofiles[iredo]
      spawn, 'rm -f ' + redofiles[iredo]
   endfor

   if (NOT keyword_set(redofiles)) then begin
      splog, 'No files to re-reduce'
   endif else begin
      ; Update the HTML file so that the files to be re-reduced
      ; are removed from the web page until then.
      apo_log2html, logfile
   endelse

   !quiet = quiet
   return
end
;------------------------------------------------------------------------------
