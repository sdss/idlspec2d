;+
; NAME:
;   batch
;
; PURPOSE:
;   Batch
;
; CALLING SEQUENCE:
;   batch,
;
; INPUTS:
;   topdir     - Local top-level directory for input and output files.
;                Also use this directory for remote hosts where REMOTEDIR
;                is not specified.
;   localfile  - Array of pointers to input files on local machine [NPROGRAM]
;   outfile    - Array of pointers to output files created on remote machine
;                and copied to local machine upon completion [NPROGRAM]
;   remotehost - List of remote hosts [NHOST]
;   remotedir  - List of remote directories; scalar or [NHOST]
;   command    - Command to execute to begin a job; scalar or [NPROGAM]
;
; OPTIONAL KEYWORDS:
;   priority   - Priority for each job, where the jobs with the largest
;                value are done first [NPROGRAM]
;   wtime      - Sleep time between checking status of all jobs; default to
;                600 seconds.
;
; OUTPUTS:
;
; COMMENTS:
;   The file names will support wildcards.  For example, if you want to
;   return all files from the directory REMOTEDIR/abc on the remote machine,
;   then set OUTFILE = 'abc/*'.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;   batch_spawn
;   create_program_list()
;   create_host_list()
;   batch_if_done()
;   batch_assign_job
;   batch_finish_job
;
; REVISION HISTORY:
;   17-Oct-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro batch_spawn, command, retval, _EXTRA=EXTRA

   splog, command
   spawn, command, retval, _EXTRA=EXTRA

   return
end

;------------------------------------------------------------------------------
function create_program_list, localfile, outfile, command, $
 priority=priority

   nprog = n_elements(localfile)

   ftemp = create_struct( name='PROGLIST_STRUCT', $
    'PROGNAME', '', $
    'LOCALFILE', ptr_new(), $
    'OUTFILE', ptr_new(), $
    'COMMAND', '', $
    'PID', 0L, $
    'PRIORITY', 1L, $
    'STATUS', 'UNASSIGNED' )

   proglist = replicate(ftemp, nprog)

   for iprog=0, nprog-1 do begin
      proglist[iprog].localfile = localfile[iprog] ; (pointer)
      proglist[iprog].outfile = outfile[iprog] ; (pointer)
      proglist[iprog].progname = 'JOB#' + strtrim(string(iprog),2)
   endfor

   proglist.command = command
   if (keyword_set(priority)) then proglist.priority = priority

   return, proglist
end

;------------------------------------------------------------------------------
function create_host_list, protocol, remotehost, remotedir, topdir

   nhost = n_elements(remotehost)

   ftemp = create_struct( name='HOSTLIST_STRUCT', $
    'PROGNAME', '', $
    'REMOTEHOST', '', $
    'REMOTEDIR', '', $
    'PROTOCOL', '', $
    'CPSTRING', '', $
    'STATUS', 'IDLE' )

   hostlist = replicate(ftemp, nhost)

   hostlist[*].remotehost = remotehost
   hostlist[*].remotedir = remotedir
   hostlist[*].protocol = protocol

   for ihost=0, n_elements(hostlist)-1 do begin
      if (keyword_set(hostlist[ihost].remotedir)) then begin
         case hostlist[ihost].protocol of
            'ssh' : hostlist[ihost].cpstring = 'scp'
            'ssh1': hostlist[ihost].cpstring = 'scp1'
            'ssh2': hostlist[ihost].cpstring = 'scp2'
            'rsh' : hostlist[ihost].cpstring = 'rcp'
            ''    : hostlist[ihost].cpstring = ''
            else  : message, 'Invalid protocol: ' + hostlist[ihost].protocol
         endcase
      endif else begin
         ; Retain CPSTRING=''
         hostlist[ihost].remotedir = topdir
      endelse
   endfor

   return, hostlist
end

;------------------------------------------------------------------------------
function batch_if_done, remotehost, remotedir, protocol, pid

   sq = "'"

   if (keyword_set(protocol)) then begin
      prothost = protocol + ' ' + remotehost + ' '
      batch_spawn, prothost + sq+' ps -p '+string(pid)+sq, retstring
   endif else begin
      batch_spawn, 'ps -p '+string(pid), retstring
   endelse

   ; If there is only one line in the return string, then this process
   ; was not found with the 'ps' command, so the job must be done.
   if (n_elements(retstring) EQ 1) then retval = 'DONE' $
    else retval = 'NOTDONE'

   splog, 'Status of job on ' + remotehost + ' = ' + retval

   return, retval
end

;------------------------------------------------------------------------------
pro batch_assign_job, ihost, iprog

   common com_batch, hostlist, proglist

   sq = "'"

   if (hostlist[ihost].status NE 'IDLE') then $
    message, 'Host is not idle'

   splog, ''
   splog, 'Assigning job ' + proglist[iprog].progname $
    + ' to host ' + hostlist[ihost].remotehost

   cpstring = hostlist[ihost].cpstring
   prothost = hostlist[ihost].protocol + ' ' + hostlist[ihost].remotehost + ' '

   if (keyword_set(cpstring)) then begin

      ; Create directories on remote machine for input files.
      ; Create all remote directories at once, then copy files into one
      ; directory at a time.
      allinput = djs_filepath(*proglist[iprog].localfile, $
       root_dir=hostlist[ihost].remotedir)
      junk = fileandpath(allinput, path=newdir)

      iuniq = uniq(newdir, sort(newdir))
      batch_spawn, prothost + 'mkdir -p ' $
       + string(newdir[iuniq]+' ',format='(99a)')

      for i=0, n_elements(iuniq)-1 do begin
         newdir1 = newdir[iuniq[i]]
         indx = where(newdir EQ newdir1)

         tmp1 = string((*proglist[iprog].localfile)[indx]+' ',format='(99a )')
         batch_spawn, cpstring + ' ' + tmp1 + ' ' $
          + hostlist[ihost].remotehost + ':' + newdir1
      endfor

      ; Create directories on remote machine for output files
      alloutput = djs_filepath(*proglist[iprog].outfile, $
       root_dir=hostlist[ihost].remotedir)
      junk = fileandpath(alloutput, path=newdir)
      iuniq = uniq(newdir, sort(newdir))
      batch_spawn, prothost + 'mkdir -p ' $
       + string(newdir[iuniq]+' ',format='(99a)')

      ; Launch the command in the background, and pipe results to /dev/null
      batch_spawn, prothost + sq+'cd '+hostlist[ihost].remotedir+'; ' $
       +proglist[iprog].command + ' >& /dev/null &'+sq, retstring

   endif else begin

      ; Only need to create local directories for output files
      junk = fileandpath(*proglist[iprog].outfile, path=newdir)
      iuniq = uniq(newdir, sort(newdir))
      batch_spawn, 'mkdir -p ' $
       + string(newdir[iuniq]+' ',format='(99a)')

      ; Launch the command in the background, and pipe results to /dev/null
      batch_spawn, proglist[iprog].command + ' >& /dev/null &', retstring

   endelse

   ; Save the process ID number, which is the 2nd word of the return string
   retwords = str_sep(retstring, ' ')
   proglist[iprog].pid = long( retwords[1] )

   hostlist[ihost].status = 'BUSY'
   hostlist[ihost].progname = proglist[iprog].progname
   proglist[iprog].status = 'RUNNING'

   return
end

;------------------------------------------------------------------------------
pro batch_finish_job, ihost, iprog

   common com_batch, hostlist, proglist

   if (hostlist[ihost].status NE 'BUSY') then $
    message, 'Host is not busy'

   splog, 'Finishing job ' + proglist[iprog].progname $
    + ' on host ' + hostlist[ihost].remotehost

   cpstring = hostlist[ihost].cpstring
   prothost = hostlist[ihost].protocol + ' ' + hostlist[ihost].remotehost + ' '

   if (keyword_set(cpstring)) then begin

      alloutput = djs_filepath(*proglist[iprog].outfile, $
       root_dir=hostlist[ihost].remotedir)

      ; Create directories on local machine for output files.
      ; Create all local directories at once, then copy files into one
      ; directory at a time.
      junk = fileandpath(*proglist[iprog].outfile, path=newdir)

      iuniq = uniq(newdir, sort(newdir))
      batch_spawn, 'mkdir -p ' $
       + string(newdir[iuniq]+' ',format='(99a)')
      for i=0, n_elements(iuniq)-1 do begin
         newdir1 = newdir[iuniq[i]]
         indx = where(newdir EQ newdir1)

         ; Copy output files from remote machine to local
         tmp1 = string(hostlist[ihost].remotehost+':' $
          +alloutput[indx]+' ',format='(99a )')
         batch_spawn, cpstring + ' ' + tmp1 + ' ' + newdir1

      endfor

      ; Remove remote output files
      batch_spawn, prothost + ' rm -f ' $
       + string(alloutput+' ',format='(99a )')

      ; Remove remote input files
      allinput = djs_filepath(*proglist[iprog].localfile, $
       root_dir=hostlist[ihost].remotedir)
      batch_spawn, prothost + ' rm -f ' $
       + string(allinput+' ',format='(99a )')

   endif

   hostlist[ihost].status = 'IDLE'
   hostlist[ihost].progname = ''
   proglist[iprog].status = 'DONE'
   proglist[iprog].pid = 0L

   return
end

;------------------------------------------------------------------------------
pro batch, topdir, localfile, outfile, protocol, remotehost, remotedir, $
 command, priority=priority, wtime=wtime

   common com_batch, hostlist, proglist

   if (NOT keyword_set(wtime)) then wtime = 600 ; Default to wait 10 mins

   if (keyword_set(topdir)) then cd, topdir

   ;----------
   ; Create a list of programs to execute (and their status)

   proglist = create_program_list(localfile, outfile, command, $ 
    priority=priority)
   nprog = n_elements(proglist)
   splog, 'Number of batch programs = ', nprog

   ;----------
   ; Create a list of available remote hosts (and their status)

   hostlist = create_host_list(protocol, remotehost, remotedir, topdir)
   nhost = n_elements(hostlist)
   splog, 'Number of hosts = ', nhost

   ;----------
   ; Find which programs are already done by looking at local files

;   for iprog=0, nprog-1 do begin
;      qdone = batch_if_done('', '', '', (*proglist[iprog].outfile)[0], $
;       proglist[iprog].endstring)
;      if (qdone EQ 'DONE') then proglist[iprog].status = 'DONE' $
;       else proglist[iprog].status = 'UNASSIGNED'
;   endfor

   ;---------------------------------------------------------------------------
   ; MAIN LOOP
   ;---------------------------------------------------------------------------

   ndone = -1
   while (ndone LT nprog) do begin

      ;----------
      ; Find any jobs that may have completed

      for ihost=0, nhost-1 do begin
         if (hostlist[ihost].status EQ 'BUSY') then begin
            j = (where(proglist.progname EQ hostlist[ihost].progname, ct))[0]
            if (ct NE 1) then $
             message, 'No or multiple program names for this host'
            qdone = batch_if_done(hostlist[ihost].remotehost, $
             hostlist[ihost].remotedir, hostlist[ihost].protocol, $
             proglist[j].pid)
            if (qdone EQ 'DONE') then $
             batch_finish_job, ihost, j
         endif
      endfor

      ;----------

      iunassign = where(proglist.status EQ 'UNASSIGNED', nunassign)
      irun = where(proglist.status EQ 'RUNNING', nrunning)
      idone = where(proglist.status EQ 'DONE', ndone)

      iidle = where(hostlist.status EQ 'IDLE', nidle)

      splog, 'Number of UNASSIGNED jobs = ', nunassign
      splog, 'Number of RUNNING jobs = ', nrunning
      splog, 'Number of DONE jobs = ', ndone

      ;----------
      ; Assign jobs, doing the highest priority jobs first

      if (nidle GT 0 AND nunassign GT 0) then begin
         for j=0, nidle-1 do begin ; Loop over available hosts
            k = (where(proglist.status EQ 'UNASSIGNED'))
            if (k[0] NE -1) then begin
               junk = max(proglist[k].priority, kmax)
               kbest = k
               batch_assign_job, iidle[j], k[kmax]
            endif
         endfor
      endif

      ;----------
      ; Sleep

      splog, 'Sleeping for ', wtime, ' seconds'
      if (ndone LT nprog) then wait, wtime

   endwhile

   print, 'All jobs have completed.'

   return
end
;------------------------------------------------------------------------------
