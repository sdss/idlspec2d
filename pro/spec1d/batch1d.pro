;+
; NAME:
;   batch1d
;
; PURPOSE:
;   Batch process Spectro-1D reductions based upon existing 2D plate files.
;
; CALLING SEQUENCE:
;   batch1d, [ fullplatefile, topdir=, upsversion=, nice=, /clobber ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   fullplatefile - Plate files to reduce; default to all files matching
;                   '*/spPlate*.fits' from the top-level directory.
;   topdir     - Top directory for reductions; default to current directory.
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
;   specified in the file TOPDIR/batch1d.par if it exists, or the default
;   file "$IDLSPEC2D_DIR/examples/batch1d.par" is used.
;
;   If using machines in Peyton, set
;     topdir='/peyton/scr/spectro0/data/2d_v4'
;   A plate is considered not reduced if any of the "spPlate*.par" files
;   do not have a corresponding "spDiag1d*.log" file.
;
;   The command is piped to the bash shell on the remote machine, so IDL
;   and the idlspec2d product must be present when running "bash --login".
;   Normally, your .bashrc file should set up all the necessary path info.
;   If the UPSVERSION keyword is used, then the UPS "setup" command must
;   also be set up in the .bashrc file.
;
;   The command that is spawned will look something like (but all in one line):
;     ssh1 wire1.princeton.edu 'cd /u/dss/spectro;
;       echo "DISPLAY=; setup idlspec2d v4_9_6; /bin/nice -n 10
;       idl 0406/spPlate-0406-51817.batch" | bash --login >& /dev/null'
;
;   The $DISPLAY environment variable is always set to "" on the remote
;   machine to make certain that we only use one IDL license per machine.
;   (Any IDL jobs that have the same the username, machine name, and $DISPLAY
;   use the same license.)
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/batch1d.par
;
; PROCEDURES CALLED:
;   djs_batch
;   djs_filepath()
;   fileandpath()
;   splog
;   yanny_free
;   yanny_read
;
; REVISION HISTORY:
;   17-Oct-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro batch1d, fullplatefile, topdir=topdir, upsversion=upsversion, $
 nice=nice, clobber=clobber

   if (NOT keyword_set(topdir)) then begin
      cd, current=topdir
   endif
   cd, topdir
   if (n_elements(nice) EQ 0) then nice = 19

   splog, prelog='(1D)'

   ;----------
   ; Create list of plate files

   if (NOT keyword_set(fullplatefile)) then $
    fullplatefile = findfile( filepath('spPlate-*.fits', root_dir='*') )
   if (NOT keyword_set(fullplatefile)) then begin
      splog, 'No plate files found'
      return
   endif
   platefile = fileandpath(fullplatefile, path=localpath)
   nplate = n_elements(platefile)

   ;----------
   ; Find which programs are already done by testing for the existing
   ; of a log file.

   platemjd = strmid(platefile, 8, 10)

   diaglog = strarr(nplate)
   qdone = bytarr(nplate)
   for iplate=0, nplate-1 do begin
      diaglog[iplate] = $
       djs_filepath('spDiag1d-' + platemjd[iplate] + '.log', $
        root_dir=localpath[iplate])
;      endstring = 'Successful completion of SPREDUCE1D'
;      spawn, 'tail -1 ' + diaglog[iplate], tailstring
;      qdone[iplate] = strpos(tailstring[0], endstring) NE -1
      qdone[iplate] = keyword_set(findfile(diaglog[iplate]))
      if (qdone[iplate]) then $
       splog, 'File ' + platefile[iplate] + ' already reduced' $
      else $
       splog, 'File ' + platefile[iplate] + ' not reduced'
   endfor

   ; If /CLOBBER is set, then reduce all plates, overwriting any
   ; previous reductions.
   if (keyword_set(clobber)) then qdone[*] = 0

   indx = where(qdone EQ 0, nplate)
   if (nplate EQ 0) then begin
      splog, 'All plates have been reduced already'
      return
   endif

   fullplatefile = fullplatefile[indx]
   platefile = platefile[indx]
   localpath = localpath[indx]
   platemjd = platemjd[indx]
   diaglog = diaglog[indx]

   ;----------
   ; Generate the IDL script files

   fq = "'"
   fullscriptfile = strarr(nplate)
   for iplate=0, nplate-1 do begin
      i = rstrpos(platefile[iplate], '.')
      if (i EQ -1) then i = strlen(platefile[iplate])
      fullscriptfile[iplate] = $
       djs_filepath(strmid(platefile[iplate],0,i)+'.batch', $
       root_dir=localpath[iplate])
      openw, olun, fullscriptfile[iplate], /get_lun
      printf, olun, 'cd, ' + fq+localpath[iplate]+fq
      printf, olun, 'spreduce1d, ' + fq+platefile[iplate]+fq
      printf, olun, 'exit'
      close, olun
      free_lun, olun
   endfor

   ;----------
   ; Create lists of input files

   infile = ptrarr(nplate)
   for i=0, nplate-1 do $
    infile[i] = ptr_new([ fullscriptfile[i],fullplatefile[i] ])

   ;----------
   ; Create lists of expected output files

   outfile = ptrarr(nplate)
   for i=0, nplate-1 do begin
      zallfile = djs_filepath('spZall-'+platemjd[i]+'.fits', $
       root_dir=localpath[i])
      zbestfile = djs_filepath('spZbest-'+platemjd[i]+'.fits', $
       root_dir=localpath[i])
      diagps = djs_filepath('spDiag1d-'+platemjd[i]+'.ps', $
       root_dir=localpath[i])
      outfile[i] = ptr_new([ diaglog[i], diagps, zallfile, zbestfile ])
   endfor

   ;----------
   ; Prioritize to do the lowest-numbered plates first

   priority = lonarr(nplate)
   isort = sort(platemjd)
   priority[isort] = reverse(lindgen(nplate)) + 1

   ;----------
   ; Determine which computers to use for these reductions.
   ; Use TOPDIR/batch1d.par if it exists, otherwise
   ; use "$IDLSPEC2D/examples/batch1d.par".

   hostfile = djs_filepath('batch1d.par', root_dir=topdir)
   hostfile = (findfile(hostfile))[0]
   if (NOT keyword_set(hostfile)) then $
    hostfile = filepath('batch1d.par', $
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
   ; Force this to be sent to a bash shell locally, and pipe to a bash shell remotely.
   ; Redirect output to /dev/null; this redirection should be valid for
   ;   either bash or csh shells.
   ; The command will look something like (but all in one line):
   ;   cd /u/dss/spectro;
   ;     echo "DISPLAY=; setup idlspec2d v4_9_6; /bin/nice -n 10
   ;     idl 0406/spPlate-0406-51817.batch" | bash --login >& /dev/null'

   setenv, 'SHELL=bash'
   precommand = 'echo "DISPLAY=; '
   if (keyword_set(upsversion)) then $
    precommand = precommand + 'setup idlspec2d ' + upsversion + '; '
   if (keyword_set(nice)) then $
    precommand = precommand + '/bin/nice -n ' + strtrim(string(nice),2)
   command = precommand + ' idl ' + fullscriptfile + '" | bash --login >& /dev/null'

   djs_batch, topdir, infile, outfile, $
    hostconfig.protocol, hostconfig.remotehost, hostconfig.remotedir, $
    command, priority=priority

   return
end
;------------------------------------------------------------------------------
