;+
; NAME:
;   nebular_batch
;
; PURPOSE:
;   Batch process the NEBULARSKY code
;
; CALLING SEQUENCE:
;   nebular_batch, [ plate, mjd=, topdir=, upsversion=, nice=, /onlysky ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   plate      - Plate number(s) to reduce; default to all non-bad plates
;                and all public plates.
;   mjd        - MJD(s) for each PLATE
;   topdir     - Top directory for reductions; default to current directory.
;   upsversion - If set, then do a "setup idlspec2d $UPSVERSION" on the 
;                remote machine before executing the IDL job.  This allows
;                you to batch jobs using a version other than that which
;                is declared current under UPS.
;   nice       - Unix nice-ness for spawned jobs; default to 19.
;   onlysky    - Keyword for NEBULARSKY
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The list of hosts and protocols should be in the Yanny parameter file
;   specified in the file TOPDIR/batch1d.par if it exists, or the default
;   file "$IDLSPEC2D_DIR/examples/batch1d.par" is used.
;   This batch processing only supports machines with cross-mounted disks.
;   It will not run if REMOTEDIR is set for any machine.
;
;   The command is piped to the bash shell on the remote machine, so IDL
;   and the idlspec2d product must be present when running "bash --login".
;   Normally, your .bashrc file should set up all the necessary path info.
;   If the UPSVERSION keyword is used, then the UPS "setup" command must
;   also be set up in the .bashrc file.
;
;   The $DISPLAY environment variable is always set to "" on the remote
;   machine to make certain that we only use one IDL license per machine.
;   (Any IDL jobs that have the same the username, machine name, and $DISPLAY
;   use the same license.)
;
;   Prioritize to do the highly-reddened plates first.
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
;   platelist
;   splog
;   yanny_readone
;
; REVISION HISTORY:
;   05-Dec-2006  Written by D. Schlegel, LBNL
;-
;------------------------------------------------------------------------------
pro nebular_batch, plate, mjd=mjd, topdir=topdir, upsversion=upsversion, $
 nice=nice, onlysky=onlysky

   if (NOT keyword_set(topdir)) then begin
      cd, current=topdir
   endif
   cd, topdir
   if (n_elements(nice) EQ 0) then nice = 19

   splog, prelog='(NEBULAR)'

   ;----------
   ; Create list of plate files

   if (keyword_set(plate)) then begin
      nplate = n_elements(plate)
      plist = replicate(create_struct('PLATE',0L,'MJD',0L), nplate)
      plist.plate = plate
      if (keyword_set(mjd1)) then begin
         if (n_elements(mjd1) NE nplate) then $
          message, 'Number of elements in PLATE and MJD do not agree'
         plist.mjd = mjd1
      endif
   endif else begin
      platelist, plist=plist
      if (keyword_set(plist)) then begin
         indx = where(strmatch(plist.status1d,'Done*') $
          AND strmatch(plist.platequality,'bad*') EQ 0 $
           OR (strtrim(plist.public) NE ''), nplate)
         if (nplate GT 0) then plist = plist[indx] $
          else plist = 0
      endif
   endelse
   if (nplate EQ 0) then begin
      splog, 'No plate files found'
      return
   endif

   ;----------
   ; Prioritize to do the highest-reddened plates first

   euler, plist.ra, plist.dec, ll, bb, 1
   priority = dust_getval(ll, bb,/ interp)

   ; Prioritize to do the lowest-numbered plates first
;   priority = lonarr(nplate)
;   isort = sort(plist.plate)
;   priority[isort] = reverse(lindgen(nplate)) + 1

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
   hostconfig = yanny_readone(hostfile)
   if (NOT keyword_set(hostconfig)) then begin
      splog, 'WARNING: Could not file batch file ' + hostfile
      return
   endif
   if (total(strtrim(hostconfig.remotedir) NE '') NE 0) then $
    message, 'This routine only supports cross-mounted machines!'

   ;----------
   ; Begin the batch jobs.
   ; Force this to be sent to a bash shell locally, and pipe to a bash shell remotely.
   ; Redirect output to /dev/null; this redirection should be valid for
   ;   either bash or csh shells.
   ; The command will look something like (but all in one line):
   ;   cd /u/dss/spectro;
   ;     echo "DISPLAY=; setup idlspec2d v4_9_6; /bin/nice -n 10
   ;     echo \"nebularsky,230,mjd=52251\" | idl " | bash --login >& /dev/null'

   platestr = strtrim(plist.plate,2)
   mjdstr = strtrim(plist.mjd,2)
   addstring = keyword_set(onlysky) ? ',/onlysky' : ''
   fq = '\"'
   setenv, 'SHELL=bash'
   precommand = 'echo "DISPLAY=; '
   if (keyword_set(upsversion)) then $
    precommand = precommand + 'setup idlspec2d ' + upsversion + '; '
   if (keyword_set(nice)) then $
    precommand = precommand + '/bin/nice -n ' + strtrim(string(nice),2)
   command = precommand + ' echo '+fq+'nebularsky,'+platestr+',mjd='+mjdstr+addstring+fq+' | idl ' + '" | bash --login >& /dev/null'

   djs_batch, topdir, 0, 0, $
    hostconfig.protocol, hostconfig.remotehost, hostconfig.remotedir, $
    command, priority=priority

   return
end
;------------------------------------------------------------------------------
