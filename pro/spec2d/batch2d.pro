;+
; NAME:
;   batch2d
;
; PURPOSE:
;   Batch process Spectro-2D reductions based upon already-built plan files.
;
; CALLING SEQUENCE:
;   batch2d, [ platenums, topdir=, mjd=, mjstart=, mjend=, nice= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   platenums  - Plate numbers to reduce.
;   topdir     - Top directory for reductions; default to current directory.
;   mjd        - MJD dates to reduce; default to all.
;   mjstart    - Starting MJD dates to reduce.
;   mjend      - Ending MJD dates to reduce.
;   nice       - Unix nice-ness for spawned jobs; default to 10.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The list of hosts and protocols should be in the Yanny parameter file
;   specified by the environment variable BATCH2DFILE, or the default file
;   "$IDLSPEC2D_DIR/examples/batch2d.par" is used.
;
;   If using machines in Peyton, set
;     topdir='/peyton/scr/spectro0/data/2d_v4'
;   A plate is considered not reduced if any of the "spPlan2d*.par" files
;   do not have a corresponding "spDiag2d*.log" file.
;
; EXAMPLES:
;
; BUGS:
;   Does not yet support MJD, MJSTART, MJEND ???
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/batch2d.par
;
; PROCEDURES CALLED:
;   djs_batch
;   djs_filepath()
;   fileandpath()
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

   i = where((*pp[0]).flavor EQ 'science')
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
pro batch2d, platenums, topdir=topdir, mjd=mjd, mjstart=mjstart, mjend=mjend, $
 nice=nice

   if (NOT keyword_set(platenums)) then platenums = '*'
   if (NOT keyword_set(topdir)) then begin
      cd, current=topdir
   endif
   cd, topdir
   if (NOT keyword_set(nice)) then nice = 10

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

   spawn, 'ls -d ' + string(platenums+' ', $
    format='(99(a," "))'), platedirs
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
; GET ONLY THE MJD's THAT WE WANT ???
         if (keyword_set(platelist)) then begin
            platelist = [platelist, platedirs[idir]]
            planlist = [planlist, planfile[ifile]]
         endif else begin
            platelist = platedirs[idir]
            planlist = planfile[ifile]
         endelse
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
      junk = fileandpath(planlist[iplate], path=thispath)
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
         printf, olun, 'setenv, ' + fq+'RAWDATA_DIR=../rawdata'+fq
         for i=0, n_elements(planfile2d)-1 do $
          printf, olun, 'spreduce2d, ' + fq+planfile2d[i]+fq
         printf, olun, 'spcombine, ' + fq+planfilecomb+fq
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
   ; Determine which computers to use for these reductions

   hostfile = getenv('BATCH2DFILE')
   if (NOT keyword_set(hostfile)) then $
    hostfile = filepath('batch2d.par', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
   yanny_read, hostfile, pp
   hostconfig = *pp[0]
   yanny_free, pp

   ;----------
   ; Begin the batch jobs.
   ; Force this to be sent to a bash shell.
   ; Set the environment variable $RAWDATA_DIR.
   ; Redirect output to /dev/null; this redirection should be valid for
   ;  either bash or csh shells.

   setenv, 'RAWDATA_DIR=../rawdata'
   setenv, 'SHELL=bash'
   nicestr = '/bin/nice -n ' + strtrim(string(nice),2)
   command = nicestr + ' idl ' + fullscriptfile + ' >& /dev/null'

   djs_batch, topdir, pinfile, poutfile, $
    hostconfig.protocol, hostconfig.remotehost, hostconfig.remotedir, $
    command, priority=priority

   ;----------
   ; Remove symbolic link to raw data

;   spawn, 'rm -f rawdata' ; But this could break another batch process!

   return
end
;------------------------------------------------------------------------------
