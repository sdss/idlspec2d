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
;   localfile  - Array of input files on local machine, where there are
;                NFILE files per program [NINFILE,NPROGRAM]
;   outfile    - Array of output files created on remote machine and copied
;                to local machine upon completion [NOUTFILE,NPROGRAM]
;   remotehost - List of remote hosts [NHOST]
;   remotedir  - List of remote directories; scalar or [NHOST]
;   command    - Command to execute to begin a job; scalar or [NPROGAM]
;   endstring  - A job is considered done when the last line of the first
;                output file contains this string; scalar or [NPROGRAM]
;
; OPTIONAL KEYWORDS:
;   wtime      - Sleep time between checking status of all jobs; default to
;                600 seconds.
;
; OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;   create_program_list()
;
; REVISION HISTORY:
;   17-Oct-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro batch_spawn, command, retval

   splog, command
   spawn, command, retval
; ???
;retval = ''

   return
end
;------------------------------------------------------------------------------
function create_program_list, localfile, outfile, command, endstring

   nprog = n_elements(localfile[0,*])

   ftemp = create_struct( name='PROGLIST_STRUCT', $
    'PROGNAME', '', $
    'LOCALFILE', ptr_new(), $
    'OUTFILE', ptr_new(), $
    'COMMAND', '', $
    'ENDSTRING', '', $
    'PRIOIRTY', 1L, $
    'STATUS', 'NOTDONE' )

   proglist = replicate(ftemp, nprog)

   for iprog=0, nprog-1 do begin
      proglist[iprog].localfile = ptr_new(localfile[*,iprog])
      proglist[iprog].outfile = ptr_new(outfile[*,iprog])
      proglist[iprog].progname = 'JOB#' + strtrim(string(iprog),2)
   endfor

   proglist.command = command
   proglist.endstring = endstring

   return, proglist
end

;------------------------------------------------------------------------------
function create_host_list, protocol, remotehost, remotedir

   nhost = n_elements(remotehost)

   ftemp = create_struct( name='HOSTLIST_STRUCT', $
    'PROGNAME', '', $
    'REMOTEHOST', '', $
    'REMOTEDIR', '', $
    'PROTOCOL', '', $
    'STATUS', 'IDLE' )

   hostlist = replicate(ftemp, nhost)

   hostlist[*].remotehost = remotehost
   hostlist[*].remotedir = remotedir
   hostlist[*].protocol = protocol

   return, hostlist
end

;------------------------------------------------------------------------------
function batch_if_done, remotehost, remotedir, protocol, outfile0, endstring

   if (keyword_set(remotehost)) then begin
      batch_spawn, protocol + ' tail -1 ' + djs_filepath(outfile0, root_dir=remotedir)
   endif else begin
      batch_spawn, 'tail -1 ' + outfile0, tailstring
   endelse

   if (strpos(tailstring[0], endstring) NE -1) then retval = 'DONE' $
    else retval = 'NOTDONE'

   splog, 'Status of ' + outfile0 + ' = ' + retval

   return, retval
end

;------------------------------------------------------------------------------
pro batch_assign_job, ihost, iprog

   common com_batch, hostlist, proglist

   if (hostlist[ihost].status NE 'IDLE') then $
    message, 'Host is not idle'

   splog, 'Assigning job ' + proglist[iprog].progname $
    + ' to host ' + hostlist[ihost].remotehost

   if (hostlist[ihost].protocol EQ 'ssh') then cpstring = 'scp' $
    else cpstring = 'rcp'
   prothost = hostlist[ihost].protocol + ' ' + hostlist[ihost].remotehost + ' '

   if (keyword_set(hostlist[ihost].remotehost)) then begin

      ; Create directories on remote machine for input files.
      ; Create one remote directory at a time, copying files over it it.
      allinput = djs_filepath(*proglist[iprog].localfile, root_dir=hostlist[ihost].remotedir)
      junk = fileandpath(allinput, path=newdir)

      iuniq = uniq(newdir, sort(newdir))
      for i=0, n_elements(iuniq)-1 do begin
         newdir1 = newdir[iuniq[i]]
         indx = where(newdir EQ newdir1)

         batch_spawn, prothost + 'mkdir -p ' + newdir1
         tmp1 = string((*proglist[iprog].localfile)[indx]+' ',format='(99a )')
         batch_spawn, cpstring + ' ' + tmp1 + ' ' $
          + hostlist[ihost].remotehost + ':' + newdir1
      endfor

      ; Create directories on remote machine for output files
      alloutput = djs_filepath(*proglist[iprog].outfile, $
       root_dir=hostlist[ihost].remotedir)
      junk = fileandpath(alloutput, path=alloutdir)
      newdir = string(alloutdir+' ',format='(99a )')
      batch_spawn, prothost + 'mkdir -p ' + newdir

      ; Launch the command in the background
      batch_spawn, prothost + proglist[iprog].command + ' &'

   endif else begin

      ; Only need to create local directories for output files
      junk = fileandpath(*proglist[iprog].outfile, path=alloutdir)
      newdir = string(alloutdir+' ',format='(99a )')
      batch_spawn, 'mkdir -p ' + newdir

      ; Launch the command in the background
stop
      batch_spawn, proglist[iprog].command + ' &'

   endelse

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

   if (hostlist[ihost].protocol EQ 'ssh') then cpstring = 'scp' $
    else cpstring = 'rcp'
   prothost = hostlist[ihost].protocol + ' ' + hostlist[ihost].remotehost + ' '

   if (keyword_set(hostlist[ihost].remotehost)) then begin

      alloutput = djs_filepath(proglist[iprog].outfile, root_dir=hostlist[ihost].remotedir)

      ; Now create local directories for output files
      junk = fileandpath(proglist[iprog].outfile, path=newdir)
      newdir = string(alloutdir+' ',format='(99a )')
      batch_spawn, 'mkdir -p ' + newdir

      ; Create directories on local machine for output files.
      ; Create one local directory at a time, copying files over it it.
      junk = fileandpath(proglist[iprog].outfile, path=newdir)

      iuniq = uniq(newdir, sort(newdir))
      for i=0, n_elements(iuniq) do begin
         newdir1 = newdir[iuniq[i]]
         indx = where(newdir EQ newdir1)

         batch_spawn, 'mkdir -p ' + newdir1
         tmp1 = string(prothost+':'+alloutput[indx]+' ',format='(99a )')
         batch_spawn, cpstring + ' ' + tmp1 + ' ' + newdir1
      endfor

   endif

   hostlist[ihost].status = 'IDLE'
   hostlist[ihost].progname = ''
   proglist[iprog].status = 'DONE'

   return
end

;------------------------------------------------------------------------------
pro batch, localfile, outfile, protocol, remotehost, remotedir, $
 command, endstring, wtime=wtime

   common com_batch, hostlist, proglist

   if (NOT keyword_set(wtime)) then wtime = 600 ; Default to wait 10 mins

   ;----------
   ; Create a list of programs to execute (and their status)

   proglist = create_program_list(localfile, outfile, command, endstring)
   nprog = n_elements(proglist)
   splog, 'Number of batch programs = ', nprog

   ;----------
   ; Create a list of available remote hosts (and their status)

   hostlist = create_host_list(protocol, remotehost, remotedir)
   nhost = n_elements(hostlist)
   splog, 'Number of hosts = ', nhost

   ;----------
   ; Find which programs are already done by looking at local files

   for iprog=0, nprog-1 do begin
      qdone = batch_if_done('', '', '', (*proglist[iprog].outfile)[0], $
       proglist[iprog].endstring)
      if (qdone EQ 'DONE') then proglist[iprog].status = 'DONE' $
       else proglist[iprog].status = 'UNASSIGNED'
   endfor

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
             (*proglist[j].outfile)[0], proglist[j].endstring)
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
      ; Assign jobs

      if (nidle GT 0 AND nunassign GT 0) then begin
         for j=0, nidle-1 do begin ; Loop over available hosts
            k = (where(proglist.status EQ 'UNASSIGNED'))[0] ; Next program
            if (k NE -1) then begin
               batch_assign_job, iidle[j], k
            endif
         endfor
      endif

      ;----------
      ; Sleep

      wait, wtime

   endwhile

   return
end

