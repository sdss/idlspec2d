;+
; NAME:
;   nebular_batch
;
; PURPOSE:
;   Batch process the NEBULARSKY code
;
; CALLING SEQUENCE:
;   nebular_batch, [ plate, mjd=, topdir=, outdir=, upsversion=, nice=, $
;    /onlysky, outfile= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   plate      - Plate number(s) to reduce; default to all non-bad plates
;                and all public plates.
;   mjd        - MJD(s) for each PLATE
;   topdir     - Top directory for reductions; default $BOSS_SPECTRO_REDUX
;   outdir     - Top directory for outfile files; default to current directory
;   upsversion - If set, then do a "setup idlspec2d $UPSVERSION" on the 
;                remote machine before executing the IDL job.  This allows
;                you to batch jobs using a version other than that which
;                is declared current under UPS.
;   nice       - Unix nice-ness for spawned jobs; default to 19.
;   onlysky    - Keyword for NEBULARSKY
;   outfile    - Output file for all plates; default to 'nebular.fits'
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
;   Should pass OUTFILE to the batch commands, but there's a quote problem???
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
pro nebular_batch, plate, mjd=mjd, topdir=topdir, outdir=outdir, $
 upsversion=upsversion, nice=nice, onlysky=onlysky, outfile=outfile

   if (NOT keyword_set(topdir)) then topdir = getenv('BOSS_SPECTRO_REDUX')
   if (NOT keyword_set(outdir)) then begin
      cd, current=outdir
   endif
   cd, outdir
   if (n_elements(nice) EQ 0) then nice = 19

   if (NOT keyword_set(outfile)) then outfile = 'nebular.fits'

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
      if (keyword_set(topdir)) then setenv, 'BOSS_SPECTRO_REDUX=' + topdir
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
   ; Determine which plates are already reduced

   olddat = hogg_mrdfits(outfile, 1, columns=['PLATE','MJD'], $
    nchunk=10, /silent)
   qexist = bytarr(nplate)
   if (keyword_set(olddat)) then begin
      for i=0L, nplate-1L do begin
         if (total(olddat.plate EQ plist[i].plate $
          AND olddat.mjd EQ plist[i].mjd) GT 0) then $
           qexist[i] = 1B
      endfor
      olddat = 0 ; clear memory
   endif

   ikeep = where(qexist EQ 0, nkeep)
   if (nkeep EQ 0) then begin
      splog, 'All plates already reduced'
      return
   endif else begin
      splog, 'Skipping ', nplate-nkeep, ' reduced plates'
      plist = plist[ikeep]
      nplate = nkeep
   endelse

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
   ;     echo "DISPLAY=; setup idlspec2d v4_9_6;
   ;     echo \"nebularsky,230,mjd=52251\" |
   ;     /bin/nice -10 idl " | bash --login >& /dev/null'

   platestr = strtrim(plist.plate,2)
   mjdstr = strtrim(plist.mjd,2)
   addstring = keyword_set(onlysky) ? ',/onlysky' : ''
   setenv, 'SHELL=bash'
   precommand = 'echo "DISPLAY=; '
   fq = '\"'
   if (keyword_set(upsversion)) then $
    precommand = precommand + 'setup idlspec2d ' + upsversion + '; '
   if (keyword_set(topdir)) then $
    precommand = precommand + 'BOSS_SPECTRO_REDUX='+fq+topdir+fq+'; '
   if (keyword_set(nice)) then nicestr = '/bin/nice -' + strtrim(string(nice),2) $
    else nicestr = ''
   command = precommand + ' echo '+fq+'nebularsky,'+platestr+',mjd='+mjdstr+addstring+fq+' | ' + nicestr + ' idl ' + '" | bash --login >& /dev/null'

   djs_batch, outdir, 0, 0, $
    hostconfig.protocol, hostconfig.remotehost, hostconfig.remotedir, $
    command, priority=priority

   return
end
;------------------------------------------------------------------------------
