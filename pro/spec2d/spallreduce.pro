;+
; NAME:
;   spallreduce
;
; PURPOSE:
;   Calling script for SPREDUCE and COMBINE2DOUT that reduces a night
;   of data according to a plan file.
;
; CALLING SEQUENCE:
;   spallreduce, [ planfile, docams=, /nocombine, /xdisplay ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name(s) of output plan file; default to reducing all
;                plan files matching 'spPlan2d*.par'
;   docams     - Cameras to reduce; default to ['b1', 'b2', 'r1', 'r2'];
;                set to 0 or '' to disable running SPREDUCE and to only
;                combine files.
;   nocombine  - Only run SPREDUCE, not COMBINE2DOUT.
;   xdisplay   - Send plots to X display rather than to plot file
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Try to avoid using SPAWN.
;   How do I just combine files? and not re-run spreduce?
;
;   What happens when the gain has changed?
;   We need to list opECalib.par in spPlan file.
;   Use default opECalib if not found...
;
; PROCEDURES CALLED:
;   combine2dout
;   cpbackup
;   idlspec2d_version()
;   idlutils_version()
;   splog
;   spreduce
;   yanny_free
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------

pro spallreduce, planfile, docams=docams, nocombine=nocombine, $
 xdisplay=xdisplay

   if (NOT keyword_set(planfile)) then planfile = findfile('spPlan2d*.par')

   ;----------
   ; If multiple plan files exist, then call this script recursively
   ; for each such plan file.

   if (N_elements(planfile) GT 1) then begin
      for i=0, N_elements(planfile)-1 do $
       spallreduce, planfile[i], docams=docams, $
        nocombine=nocombine, xdisplay=xdisplay
      return
   endif

   if (N_elements(docams) EQ 0) then docams = ['b1', 'b2', 'r1', 'r2']
   if (keyword_set(docams)) then ndocam = N_elements(docams) $
    else ndocam = 0

   yanny_read, planfile[0], pdata, hdr=hdr

   ;----------
   ; Find the ONEEXP structure

   for i=0, N_elements(pdata)-1 do begin
      if (tag_names(*pdata[i], /structure_name) EQ 'PIXFLATS') then $
       pixflats = *pdata[i]
      if (tag_names(*pdata[i], /structure_name) EQ 'ONEPLUG') then $
       allplug = *pdata[i]
      if (tag_names(*pdata[i], /structure_name) EQ 'ONEEXP') then $
       allseq = *pdata[i]
   endfor

   yanny_free, pdata

   if (N_elements(allseq) EQ 0) then begin
      splog, 'ABORT: No ONEEXP structures in plan file ' + planfile
      return
   endif

   ;----------
   ; Find keywords from the header

   inputDir = yanny_par(hdr, 'inputDir')
   plugDir = yanny_par(hdr, 'plugDir')
   flatDir = yanny_par(hdr, 'flatDir')
   extractDir = yanny_par(hdr, 'extractDir')
   combineDir = yanny_par(hdr, 'combineDir')
   logfile = yanny_par(hdr, 'logfile')
   plotfile = yanny_par(hdr, 'plotfile')

   mjd = yanny_par(hdr, 'MJD')
   run = yanny_par(hdr, 'run')
   run = strtrim(string(run), 2)

   spawn, 'mkdir -p '+extractDir

   stime0 = systime(1)

   ;----------
   ; Open log files for output

   if (keyword_set(logfile)) then begin
      cpbackup, logfile
      splog, filename=logfile
      splog, 'Log file ', logfile, ' opened ', systime()
   endif
   if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then begin
      cpbackup, plotfile
      set_plot, 'ps'
      device, filename=plotfile, /color
      splog, 'Plot file ', plotfile, ' opened ', systime()
   endif
   splog, 'Plan file ', planfile
   splog, 'DOCAMS = ', docams

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()

   camnames = ['b1', 'b2', 'r1', 'r2']
   specnums = [1, 2, 1, 2] ; Spectrograph ID corresponding to CAMNAMES above
   ncam = N_elements(camnames)

   ;----------
   ; Find all the sequence IDs

   seqid = allseq[ sort(allseq.seqid) ].seqid
   seqid = allseq[ uniq(allseq.seqid) ].seqid

   for iseq=0, N_elements(seqid)-1 do begin

      ;----------
      ; Get the plate ID number from any (e.g., the first) exposure with
      ; this sequence ID number

      j = where(allseq.seqid EQ seqid[iseq])
      plateid = allseq[j[0]].plateid
      platestr = string(plateid, format='(i4.4)')

      stime1 = systime(1)
      splog, 'Begin plate ' + platestr + ' at ' + systime()

      ;----------
      ; Find the corresponding plug map file

      j = where(allplug.seqid EQ seqid[iseq] $
            AND allplug.plateid EQ plateid )
      if (j[0] NE -1) then begin
         plugfile = allplug[j[0]].name
      endif else begin
         splog, 'ABORT: No plug map file for SEQID= ' $
          + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' + platestr
         return
      endelse
      splog, 'Plug map file = ', plugfile

      for ido=0, ndocam-1 do begin

         icam = (where(camnames EQ docams[ido], camct))[0]
         splog, prelog=camnames[icam]
         if (camct NE 1) then message, 'Non-unique camera ID '

         ;----------
         ; Find the corresponding pixel flat

         pixflatname = pixflats.name[icam]
         if (pixflatname EQ 'UNKNOWN') then $
          message, 'No pixel flat for CAMERA= ' + camnames[icam]

         j = where(allseq.seqid EQ seqid[iseq] $
               AND allseq.flavor EQ 'science' $
               AND allseq.name[icam] NE 'UNKNOWN' )

         if (j[0] NE -1) then begin

            ; String array with all science exposures at this sequence + camera
            objname = allseq[j].name[icam]

            ;-----------
            ; Select **all** flat exposures at this sequence + camera

            j = where(allseq.seqid EQ seqid[iseq] $
                  AND allseq.flavor EQ 'flat' $
                  AND allseq.name[icam] NE 'UNKNOWN', nflat )
            if (nflat GT 0) then begin
               flatname = allseq[j].name[icam]
            endif else begin
               flatname = ''
               splog, 'ABORT: No flat for SEQID= ' $
                + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
                + platestr + ', CAMERA= ' + camnames[icam]
            endelse

            ;-----------
            ; Select **all** arc exposures at this sequence + camera

            j = where(allseq.seqid EQ seqid[iseq] $
                  AND allseq.flavor EQ 'arc' $
                  AND allseq.name[icam] NE 'UNKNOWN', narc )
            if (narc GT 0) then begin
               arcname = allseq[j].name[icam]
            endif else begin
               arcname = ''
               splog, 'ABORT: No arc for SEQID= ' $
                + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
                + platestr + ', CAMERA= ' + camnames[icam]
            endelse

            ;----------
            ; Get full name of pixel flat

            pixflatname = filepath(pixflatname, root_dir=flatDir)

            stime2 = systime(1)

            ;----------
            ; Reduce this set of frames (all objects w/same plate + camera)

            if (keyword_set(arcname) AND keyword_set(flatname)) then begin
               plottitle = 'MJD='+strtrim(string(mjd),2) $
                + ' PLATE='+platestr + ': '

               spreduce, flatname, arcname, objname, $
                pixflatname=pixflatname, plugfile=plugfile, lampfile=lampfile, $
                indir=inputDir, plugdir=plugDir, outdir=extractDir, $
                summarystruct=summarystruct, plottitle=plottitle
            endif

            splog, 'Time to reduce camera ', camnames[icam], ' = ', $
             systime(1)-stime2, ' seconds', format='(a,a,a,f6.0,a)'

            heap_gc   ; garbage collection
         endif

         splog, prelog=''
      endfor ; End loop for camera number

      splog, 'Time to reduce all cameras = ', $
       systime(1)-stime1, ' seconds', format='(a,f6.0,a)'

      ;----------
      ; Combine all red+blue exposures for a given sequence

      if (NOT keyword_set(nocombine)) then begin

         spawn, 'mkdir -p '+combineDir

         stime3 = systime(1)

         for side=1, 2 do begin

            j = where(allseq.seqid EQ seqid[iseq])
            icam = where(specnums EQ side)
            rawfiles = allseq[j].name[icam] ; Raw file names
            j = where(rawfiles NE 'UNKNOWN', nfile)
            if (nfile GT 0) then rawfiles = rawfiles[j]

            ; Now search for any extracted spectra on disk, overwriting
            ; FILES[i] with '' if a file does not exist.
            ; Add an 's' because we now output '.fits'

            files = 'spSpec2d-' + strmid(rawfiles, 4)+'s'

            for i=0, nfile-1 do $
             files[i] = findfile(filepath(files[i], root_dir=extractDir))
            j = where(files NE '', nfile)
            if (nfile GT 0) then files = files[j] $
             else files = 0

            splog, 'Combining ' + strtrim(string(nfile),2) $
             + ' files for side ' + strtrim(string(side),2)

            if (nfile GT 0) then begin
               outputroot = string('spMerge2d-',mjd,'-',platestr, $
                format='(a,i5.5,a1,a4)')
               individualroot  = string('spInd2d-',mjd,'-',platestr, $
                format='(a,i5.5,a1,a4)')

               for i=0, nfile-1 do splog, 'Combine file ', files[i]
               combine2dout, files, filepath(outputroot, root_dir=combineDir), $
                side, wavemin=alog10(3750.0d), window=100, $
                individual=filepath(individualroot, root_dir=combineDir)
            endif

         endfor

         splog, 'Time to combine sequence ', seqid[iseq], ' = ', $
          systime(1)-stime3, ' seconds', format='(a,i5,a,f6.0,a)'
      endif

   endfor ; End loop for sequence number

   splog, 'Total time for SPALLREDUCE = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPALLREDUCE at ', systime()

   ;----------
   ; Close log files

   if (keyword_set(logfile)) then splog, /close
   if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then begin
      device, /close
      set_plot, 'x'
   endif

   return
end
;------------------------------------------------------------------------------
