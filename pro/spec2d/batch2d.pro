;+
; NAME:
;   batch2d
;
; PURPOSE:
;   Batch process Spectro-2D reductions based upon already-built plan files.
;
; CALLING SEQUENCE:
;   batch2d, [ platenums, topdir=, platestart=, plateend=, $
;    mjd=, mjstart=, mjend=, upsversion=, nice=, /clobber ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   platenums  - Plate numbers to reduce; default to '*'
;   topdir     - Top directory for reductions; default to current directory.
;   platestart - Starting plate number.
;   plateend   - Ending plate number.
;   mjd        - MJD dates to reduce; default to all.
;                Select based upon the MJD of the combine plan file, and
;                reduce data from all nights needed for that combined plate+MJD.
;   mjstart    - Starting MJD dates to reduce.
;   mjend      - Ending MJD dates to reduce.
;   upsversion - If set, then do a "setup idlspec2d $UPSVERSION" on the
;                remote machine before executing the IDL job.  This allows
;                you to batch jobs using a version other than that which
;                is declared current under UPS.
;   nice       - Unix nice-ness for spawned jobs; default to 19.
;   clobber    - If set, then reduce all specified plates, overwriting
;                any previous reductions.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The list of hosts and protocols should be in the Yanny parameter file
;   specified in the file TOPDIR/batch2d.par if it exists, or the default
;   file "$IDLSPEC2D_DIR/examples/batch2d.par" is used.
;
;   If using machines in Peyton, set
;     topdir='/peyton/scr/spectro0/data/2d_v4'
;   A plate is considered not reduced if any of the "spPlan2d*.par" files
;   do not have a corresponding "spDiag2d*.log" file.
;
;   The command is piped to the bash shell on the remote machine, so IDL
;   and the idlspec2d product must be present when running "bash --login".
;   Normally, your .bashrc file should set up all the necessary path info.
;   If the UPSVERSION keyword is used, then the UPS "setup" command must
;   also be set up in the .bashrc file.
;
;   The command that is spawned will look something like (but all in one line):
;     ssh1 wire1.princeton.edu 'cd /u/dss/spectro;
;       echo "DISPLAY=:0.0; setup idlspec2d v4_9_6; /bin/nice -n 10
;       idl 0406/spPlancomb-0406-51817.batch" | bash --login >& /dev/null'
;
;   The $DISPLAY environment variable is always set to ":0.0" on the remote
;   machine to make certain that we only use one IDL license per machine.  
;   (Any IDL jobs that have the same the username, machine name, and $DISPLAY 
;   use the same license.)
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/batch2d.par
;
; PROCEDURES CALLED:
;   concat_dir()
;   djs_batch
;   djs_filepath()
;   fileandpath()
;   get_mjd_dir()
;   mjd_match()
;   repstr()
;   splog
;   yanny_free
;   yanny_read
;   yanny_par()
;
; INTERNAL SUPPORT ROUTINES:
;   batch2d_nolog()
;   batch2d_rawfiles()
;   batch2d_combfiles()
;
; REVISION HISTORY:
;   17-Oct-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
; For a list of files with names like 'path/spPlanXXX.par', 
; return the indexes of all that do **not** have a corresponding 
; log file of the form 'path/spDiagXXX.log'.

function batch2d_nolog, planfile

   retindx = -1L
   for ifile=0, n_elements(planfile)-1 do begin
      thisfile = fileandpath(planfile[ifile], path=thispath)
      kk = rstrpos(thisfile, '.')
      logfile = 'spDiag' + strmid(thisfile, 6, kk-6) + '.log'
      logfile = djs_filepath(logfile, root_dir=thispath)
      if (NOT keyword_set(findfile(logfile))) then retindx = [retindx, ifile]
   endfor

   nfound = n_elements(retindx)-1
   if (nfound EQ 0) then return, retindx $
    else return, retindx[1:nfound]
end

;------------------------------------------------------------------------------
function batch2d_rawfiles, planfile, outfile=outfile

   nplan = n_elements(planfile)
   if (nplan GT 1) then begin
      infiles = batch2d_rawfiles(planfile[0], outfile=outfile)
      for i=1, nplan-1 do begin
         infiles = [infiles, $
          batch2d_rawfiles(planfile[i], outfile=tmpout)]
         outfile = [outfile, tmpout]
      endfor
      return, infiles
   endif

   yanny_read, planfile, pp, hdr=hdr
   if (NOT keyword_set(pp)) then begin
      splog, 'WARNING: Could not find plan file ' + planfile
      return, ''
   endif
   extractdir = yanny_par(hdr, 'extractdir')
   thismjd = long(yanny_par(hdr, 'MJD'))
   logfile = yanny_par(hdr, 'logfile')
   plotfile = yanny_par(hdr, 'plotfile')
   for ii=0, n_elements(pp)-1 do begin
      sname = tag_names(*pp[ii], /structure_name)
      if (sname EQ 'SPEXP') then begin
         expfiles = ((*pp[ii]).name)
      endif
   endfor
   if (NOT keyword_set(expfiles)) then begin
      splog, 'WARNING: Could not find any exposures in plan file ' + planfile
      return, ''
   endif

   ; Add a wildcard to the end of the raw FITS files so that we
   ; find compressed (.gz) files.
   mjdstr = string(thismjd, format='(i05.5)')
   infiles = djs_filepath(expfiles[*]+'*', root_dir=mjdstr)

   ;----------
   ; Replace the prefix 'sdR' with 'spSky' in the science frames,
   ; with 'spArc' for the arc frames, and with 'spFlat' for the flat frames.
   ; Also, force the suffix from '.fit' to '.fits'.

   otherfiles = expfiles
   if (size(otherfiles,/n_dimen) EQ 1) then nloop=1 $
    else nloop = (size(otherfiles,/dimens))[1]
   for i=0, nloop-1 do begin
      case ((*pp[0]).flavor)[i] of
         'arc': otherfiles[*,i] = repstr(otherfiles[*,i], 'sdR', 'spArc')
         'flat': otherfiles[*,i] = repstr(otherfiles[*,i], 'sdR', 'spFlat')
         'science': otherfiles[*,i] = repstr(otherfiles[*,i], 'sdR', 'spSky')
         'smear': otherfiles[*,i] = repstr(otherfiles[*,i], 'sdR', 'spSky')
      endcase
      otherfiles[*,i] = repstr(otherfiles[*,i], '.fit', '.fits')
   endfor

   ;----------
   ; Now get the reduced frames.

   i = where((*pp[0]).flavor EQ 'science' OR (*pp[0]).flavor EQ 'smear')
   framefiles = expfiles[*,i]
   framefiles = repstr(framefiles, 'sdR', 'spFrame')
   framefiles = repstr(framefiles, '.fit', '.fits')

   outfile = [ djs_filepath(logfile, root_dir=extractdir), $
               djs_filepath(plotfile, root_dir=extractdir), $
               djs_filepath(framefiles[*], root_dir=extractdir), $
               djs_filepath(otherfiles[*], root_dir=extractdir) ]

   yanny_free, pp

   return, infiles
end

;------------------------------------------------------------------------------
function batch2d_combfiles, planfile, outfile=outfile

   nplan = n_elements(planfilecomb)
   if (nplan GT 1) then begin
      rawfiles = batch2d_combfiles(planfile[0], outfile=outfile)
      for i=1, nplan-1 do begin
         infiles = [infiles, $
          batch2d_rawfiles(planfile[i], outfile=tmpout)]
         outfile = [outfile, tmpout]
      endfor
      return, infiles
   endif

   yanny_read, planfile, pp, hdr=hdr
   if (NOT keyword_set(pp)) then begin
      splog, 'WARNING: Could not find plan file ' + planfile
      return, ''
   endif
   extractdir = yanny_par(hdr, 'extractdir')
   combinedir = yanny_par(hdr, 'combinedir')
   logfile = yanny_par(hdr, 'logfile')
   plotfile = yanny_par(hdr, 'plotfile')
   plotsnfile = yanny_par(hdr, 'plotsnfile')
   fcalibprefix = yanny_par(hdr, 'fcalibprefix')
   combinefile = yanny_par(hdr, 'combinefile')
   for ii=0, n_elements(pp)-1 do begin
      sname = tag_names(*pp[ii], /structure_name)
      if (sname EQ 'SPEXP') then begin
         expfiles = ((*pp[ii]).name)[*]
         expnums = strmid((*pp[ii]).name[0], 11, 8)
      endif
   endfor
   if (NOT keyword_set(expfiles)) then begin
      splog, 'WARNING: Could not find any exposures in plan file ' + planfile
      return, ''
   endif

   fcorrprefix = 'spFluxcorr-' + expnums

   infiles = [ djs_filepath(expfiles, root_dir=mjdstr) ]

   junk = fileandpath(planfile, path=thisdir)
   outfile = [ djs_filepath(logfile, root_dir=thisdir), $
               djs_filepath(plotfile, root_dir=thisdir), $
               djs_filepath(plotsnfile, root_dir=thisdir), $
               djs_filepath(fcalibprefix+'-*.fits', root_dir=thisdir), $
               djs_filepath(fcorrprefix+'-*.fits', root_dir=thisdir), $
               djs_filepath(combinefile, root_dir=thisdir) ]

   yanny_free, pp

   return, infiles
end

;------------------------------------------------------------------------------
pro batch2d, platenums, topdir=topdir, $
 platestart=platestart, plateend=plateend, $
 mjd=mjd, mjstart=mjstart, mjend=mjend, $
 upsversion=upsversion, nice=nice, clobber=clobber

   if (NOT keyword_set(platenums)) then platenums = '*'
   if (NOT keyword_set(topdir)) then begin
      cd, current=topdir
   endif
   cd, topdir
   if (n_elements(nice) EQ 0) then nice = 19

   splog, prelog='(2D)'

   ;----------
   ; Create symoblic link from current directory to raw data directory

   rawdata_dir = getenv('RAWDATA_DIR')
   if (NOT keyword_set(rawdata_dir)) then $
    message, 'Must set environment variable RAWDATA_DIR'
   junk = findfile('rawdata', count=ct)
   if (ct EQ 0) then $
    spawn, 'ln -s ' + rawdata_dir + ' rawdata'

   ;----------
   ; Create list of plate directories
   ; Limit the list to only those specified by PLATENUMS,PLATESTART,PLATEEND

;   if (size(platenums,/tname) EQ 'INT' OR size(platenums,/tname) EQ 'LONG') $
;    then platestr = string(platenums,format='(i4.4)') $
;   else platestr = platenums
;   spawn, 'ls -d ' + string(platestr+' ', $
;    format='(99(a," "))'), platedirs

   platedirs = get_mjd_dir(topdir, mjd=platenums, mjstart=platestart, $
    mjend=plateend)

   if (NOT keyword_set(platedirs[0])) then begin
      splog, 'No directories found'
      return
   endif
   ndir = n_elements(platedirs)

   ;----------
   ; In each plate directory, find all 'spPlancomb*.par' files

   for idir=0, ndir-1 do begin
      planfile = findfile( $
       djs_filepath('spPlancomb*.par', root_dir=platedirs[idir]), count=nfile)

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

   ;----------
   ; Create pointers for input and output files

   pinfile = ptrarr(nplate)
   poutfile = ptrarr(nplate)

   ;----------
   ; For each combine plan file, generate the IDL script files

   fq = "'"
   fullscriptfile = strarr(nplate)
   for iplate=0, nplate-1 do begin
      ; Find all relevant 2D plan files
      yanny_read, planlist[iplate], hdr=hdr
      planfile2d = yanny_par(hdr, 'planfile2d')

      ; Find which of these plan files do **not** have log files already.
      ; Presume those are the ones that need to be reduced.
      ; The exception is if /CLOBBER is set -- in that case, we reduce
      ; all specified plates, overwriting any previous reductions.
      junk = fileandpath(planlist[iplate], path=thispath)
      if (keyword_set(clobber)) then $
       ido2d = lindgen(n_elements(planfile2d)) $
      else $
       ido2d = batch2d_nolog(djs_filepath(planfile2d, root_dir=thispath))

      if (ido2d[0] NE -1) then begin
         ; Trim the list of 2D plan files to those not reduced yet.
         planfile2d = planfile2d[ido2d]

         ; Split the combine plan file name into a directory and file name
         planfilecomb = fileandpath(planlist[iplate], path=pathcomb)

         ; Construct the name of the batch file
         i = rstrpos(planfilecomb, '.')
         if (i EQ -1) then i = strlen(planfilecomb)
         fullscriptfile[iplate] = $
          djs_filepath(strmid(planfilecomb,0,i)+'.batch', root_dir=pathcomb)

         ; Write the batch file
         openw, olun, fullscriptfile[iplate], /get_lun
         printf, olun, '; Auto-generated batch file '+systime()
         printf, olun, 'cd, ' + fq+pathcomb+fq
;         printf, olun, 'setenv, ' + fq+'RAWDATA_DIR=../rawdata'+fq
         printf, olun, 'setenv, ' + fq+'RAWDATA_DIR=' $
          +concat_dir(topdir,'rawdata')+fq
         for i=0, n_elements(planfile2d)-1 do $
          printf, olun, 'spreduce2d, ' + fq+planfile2d[i]+fq
         printf, olun, 'spcombine, ' + fq+planfilecomb+fq
         printf, olun, 'spawn, '+fq+'rm -f spArc*.fits* spFlat*.fits* spSky*.fits*'+fq
         printf, olun, 'exit'
         close, olun
         free_lun, olun

         ; List of input files
         planfile2d = filepath(planfile2d, root_dir=pathcomb)
         rawfiles = batch2d_rawfiles(planfile2d, outfile=outfile2d)
         rawfiles = djs_filepath(rawfiles, root_dir='rawdata')
         outfile2d = djs_filepath(outfile2d, root_dir=pathcomb)
         junk = batch2d_combfiles(planlist[iplate], outfile=outfilecomb)

         pinfile[iplate] = ptr_new([ fullscriptfile[iplate], $
          planlist[iplate], planfile2d, rawfiles ])

         ; List of output files
         poutfile[iplate] = ptr_new([ outfile2d, outfilecomb ])
      endif

   endfor

   ;----------
   ; Trim the plate list to only those needing reductions.

   iplate = where(pinfile NE ptr_new(), nplate)
   if (iplate[0] EQ -1) then begin
      splog, 'All plates have been reduced'
      return
   endif

   platelist = platelist[iplate]
   pinfile = pinfile[iplate]
   poutfile = poutfile[iplate]
   fullscriptfile = fullscriptfile[iplate]

   ;----------
   ; Prioritize to do the lowest-numbered plates first

   priority = lonarr(nplate)
   isort = sort(platelist)
   priority[isort] = reverse(lindgen(nplate)) + 1

   ;----------
   ; Determine which computers to use for these reductions.
   ; Use TOPDIR/batch2d.par if it exists, otherwise
   ; use "$IDLSPEC2D/examples/batch2d.par".

   hostfile = djs_filepath('batch2d.par', root_dir=topdir)
   hostfile = (findfile(hostfile))[0]
   if (NOT keyword_set(hostfile)) then $
    hostfile = filepath('batch2d.par', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
   splog, 'Reading batch file ' + hostfile
   yanny_read, hostfile, pp
   if (NOT keyword_set(pp)) then begin
      splog, 'WARNING: Could not file batch file ' + hostfile
      return
   endif
   hostconfig = *pp[0]
   yanny_free, pp

   ;----------
   ; Begin the batch jobs.
   ; Force this to be sent to a bash shell locally, and pipe to bash remotely.
   ; Set the environment variable $RAWDATA_DIR.
   ; Redirect output to /dev/null; this redirection should be valid for
   ;   either bash or csh shells.
   ; The command will look something like (but all in one line):
   ;   cd /u/dss/spectro;
   ;     echo "DISPLAY=:0.0; setup idlspec2d v4_9_6; /bin/nice -n 10
   ;     idl 0406/spPlancomb-0406-51817.batch" | bash --login >& /dev/null'

   setenv, 'RAWDATA_DIR=../rawdata'
   setenv, 'SHELL=bash'
   precommand = 'echo "DISPLAY=:0.0; '
   if (keyword_set(upsversion)) then $
    precommand = precommand + 'setup idlspec2d ' + upsversion + '; '
   if (keyword_set(nice)) then $
    precommand = precommand + '/bin/nice -n ' + strtrim(string(nice),2)
   command = precommand + ' idl ' + fullscriptfile + '" | bash --login >& /dev/null'

   djs_batch, topdir, pinfile, poutfile, $
    hostconfig.protocol, hostconfig.remotehost, hostconfig.remotedir, $
    command, priority=priority

   ;----------
   ; Remove symbolic link to raw data

;   spawn, 'rm -f rawdata' ; But this could break another batch process!

   return
end
;------------------------------------------------------------------------------
