;+
; NAME:
;   batch1d
;
; PURPOSE:
;   Batch process Spectro-1D reductions based upon existing 2D plate files.
;
; CALLING SEQUENCE:
;   batch1d, [ fullplatefile, topdir=, nice= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   fullplatefile - Plate files to reduce; default to all files matching
;                   '*/spPlate*.fits' from the top-level directory.
;   topdir     - Top directory for reductions; default to current directory.
;   nice       - Unix nice-ness for spawned jobs; default to 10.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The list of hosts and protocols should be in the Yanny parameter file
;   specified by the environment variable BATCH2DFILE, or the default file
;   "$IDLSPEC2D_DIR/examples/batch1d.par" is used.
;
;   If using machines in Peyton, set
;     topdir='/peyton/scr/spectro0/data/2d_v4'
;   A plate is considered not reduced if any of the "spPlate*.par" files
;   do not have a corresponding "spDiag1d*.log" file.
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
pro batch1d, fullplatefile, topdir=topdir, nice=nice

   if (NOT keyword_set(topdir)) then begin
      cd, current=topdir
   endif
   cd, topdir
   if (NOT keyword_set(nice)) then nice = 10

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
;   endstring = 'Successful completion of SPREDUCE1D'
   for iplate=0, nplate-1 do begin
      diaglog[iplate] = $
       djs_filepath('spDiag1d-' + platemjd[iplate] + '.log', $
        root_dir=localpath[iplate])
;      spawn, 'tail -1 ' + diaglog[iplate], tailstring
;      qdone[iplate] = strpos(tailstring[0], endstring) NE -1
      qdone[iplate] = keyword_set(findfile(diaglog[iplate]))
      if (qdone[iplate]) then $
       splog, 'File ' + platefile[iplate] + ' already reduced' $
      else $
       splog, 'File ' + platefile[iplate] + ' not reduced'
   endfor

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
   ; Determine which computers to use for these reductions

   hostfile = getenv('BATCH1DFILE')
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
   ; Force this to be sent to a bash shell.
   ; Redirect output to /dev/null; this redirection should be valid for
   ;  either bash or csh shells.

   setenv, 'SHELL=bash'
   nicestr = '/bin/nice -n ' + strtrim(string(nice),2)
   command = nicestr + ' idl ' + fullscriptfile + ' >& /dev/null'

   djs_batch, topdir, infile, outfile, $
    hostconfig.protocol, hostconfig.remotehost, hostconfig.remotedir, $
    command, priority=priority

   return
end
;------------------------------------------------------------------------------
