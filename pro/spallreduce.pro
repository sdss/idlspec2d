;+
; NAME:
;   spallreduce
;
; PURPOSE:
;   Calling script for SPREDUCE and COMBINE2DOUT that reduces a night
;   of data according to a plan file.
;
; CALLING SEQUENCE:
;   spallreduce, planfile=planfile, [ combineonly=combineonly, docams=docams ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name of output plan file; default to 'spPlan2d.par'
;   docams     - Cameras to reduce; default to ['b1', 'r2', 'b2', 'r1']
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Rather than use FINDFILE, we should know what all the idlout* files are.
;   Need to implement lampfile.
;   Pass QA flags to SPREDUCE
;   Try to avoid using SPAWN.
;   Should pair up flats+arcs like Scott had done previously before passing
;     to SPREDUCE.
;
; PROCEDURES CALLED:
;   combine2dout
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

pro spallreduce, planfile=planfile, combineonly=combineonly, docams=docams, $
        display=display

   if (NOT keyword_set(planfile)) then planfile = 'spPlan2d.par'

   docomb = 0
   if (NOT keyword_set(docams)) then begin
      docams = ['b1', 'r2', 'b2', 'r1'] ; do all cameras
      docomb = 1
   endif
   ndo = N_elements(docams)

   yanny_read, planfile, pdata, hdr=hdr

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

   if (N_elements(allseq) EQ 0) then $
    message, 'No ONEEXP structures in plan file ' + planfile

   ; Find keywords from the header
   plugDir = yanny_par(hdr, 'plugDir')
   extractDir = yanny_par(hdr, 'extractDir')
   combDir = yanny_par(hdr, 'combDir')
   if (NOT keyword_set(combDir)) then combDir=extractDir

   inputDir = yanny_par(hdr, 'inputDir')
   flatDir = yanny_par(hdr, 'flatDir')
   mjd = yanny_par(hdr, 'MJD')
   run = yanny_par(hdr, 'run')
   run = strtrim(string(run), 2)

   logfile = yanny_par(hdr, 'logfile')
   if (keyword_set(logfile)) then $
    logfile = filepath(logfile, root_dir=extractDir)
   plotfile = yanny_par(hdr, 'plotfile')
   if (keyword_set(plotfile)) then $
    plotfile = filepath(plotfile, root_dir=extractDir)

   spawn, 'mkdir -p '+extractDir

   verybegin = systime(1)
   ; Open log files for output
   if (keyword_set(logfile)) then begin
      splog, filename=logfile, /append
      splog, 'Log file ', logfile, ' opened ', systime()
   endif
   if (keyword_set(plotfile)) then begin
      set_plot, 'ps'
      device, filename=plotfile, /color
      splog, 'Plot file ', plotfile, ' opened ', systime()
   endif

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()

   camnames = ['b1', 'r2', 'b2', 'r1']
   camnums = ['01', '02', '03', '04']
   ncam = N_elements(camnames)

   ; Find all the sequence IDs
   seqid = allseq[ sort(allseq.seqid) ].seqid
   seqid = allseq[ uniq(allseq.seqid) ].seqid

   for iseq=0, N_elements(seqid)-1 do begin

      ; Get the plate ID number from any (e.g., the first) exposure with
      ; this sequence ID number
      j = where(allseq.seqid EQ seqid[iseq])
      plateid = allseq[j[0]].plateid
      splog, 'Begin plate ', strtrim(string(plateid),2), ' at ', systime()

      plateDir=filepath(strtrim(string(plateid),2)+'/2d_'+ $
       run,root_dir=extractDir)
      combineDir=filepath(strtrim(string(plateid),2)+'/comb_'+ $
       run,root_dir=combDir)

      if (NOT keyword_set(combineonly)) then begin

      ; Find the corresponding plug map file
      j = where(allplug.seqid EQ seqid[iseq] $
            AND allplug.plateid EQ plateid )
      if (j[0] NE -1) then plugfile = allplug[j[0]].name $
       else message, 'No plug map file for SEQID= ' $
        + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
        + strtrim(string(plateid),2)

      for ido=0, ndo-1 do begin

         icam = (where(camnames EQ docams[ido], camct))[0]
         splog, camname=camnames[icam]
         if (camct NE 1) then message, 'Non-unique camera ID '

         ; Find the corresponding pixel flat
         pixflatname = pixflats.name[icam]
         if (pixflatname EQ 'UNKNOWN') then $
          message, 'No pixel flat for CAMERA= ' + camnums[icam]

         j = where(allseq.seqid EQ seqid[iseq] $
               AND allseq.flavor EQ 'science' $
               AND allseq.name[icam] NE 'UNKNOWN' )

         if (j[0] NE -1) then begin

            ; String array with all science exposures at this sequence + camera
            objname = allseq[j].name[icam]

            ;------
            ; Select **all** flat exposures at this sequence + camera

            j = where(allseq.seqid EQ seqid[iseq] $
                  AND allseq.flavor EQ 'flat' $
                  AND allseq.name[icam] NE 'UNKNOWN', nflat )
            if (nflat GT 0) then flatname = allseq[j].name[icam] $
             else message, 'No flat for SEQID= ' $
              + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
              + strtrim(string(plateid),2) + ', CAMERA= ' + camnums[icam]

            ;------
            ; Select **all** arc exposures at this sequence + camera

            j = where(allseq.seqid EQ seqid[iseq] $
                  AND allseq.flavor EQ 'arc' $
                  AND allseq.name[icam] NE 'UNKNOWN', narc )
            if (narc GT 0) then arcname = allseq[j].name[icam] $
             else message, 'No arc for SEQID= ' $
              + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
              + strtrim(string(plateid),2) + ', CAMERA= ' + camnums[icam]

            ;-----
            ; Read the time of each flat + arc

            arctime = dblarr(narc)
            for iarc=0, narc-1 do begin
               sdssproc, arcname[iarc], indir=inputDir, hdr=hdr
               arctime[iarc] = sxpar(hdr, 'TAI')
            endfor

            flattime = dblarr(nflat)
            for iflat=0, nflat-1 do begin
               sdssproc, flatname[iflat], indir=inputDir, hdr=hdr
               flattime[iflat] = sxpar(hdr, 'TAI')
            endfor

            ;-----
            ; Enforce a one-to-one match between arcs + flats.
            ; Do this by pairing each arc with the flat w/ the nearest SEQID.
            ; This will create a list of flats that is the same length
            ; as the list of arcs

            flatsortname = strarr(narc)
            for iarc=0, narc-1 do begin
               junk = min( abs(flattime - arctime[iarc]), iclose)
               flatsortname[iarc] = flatname[iclose]
               splog, 'Pair arc ', arcname[iarc], ' with flat ', flatname[iclose]
            endfor

            ; Get full name of pixel flat
            pixflatname = filepath(pixflatname, root_dir=flatDir)

            spawn, 'mkdir -p '+plateDir

            spreduce, flatsortname, arcname, objname, $
             pixflatname=pixflatname, plugfile=plugfile, lampfile=lampfile, $
             indir=inputDir, plugdir=plugDir, outdir=plateDir

            heap_gc   ; garbage collection
         endif

         splog, camname=''
      endfor ; End loop for camera number

   endif

   ; Combine all red+blue exposures for a given sequence

   if (docomb) then begin

      spawn, 'mkdir -p '+combineDir
      startcombtime = systime(1)

      for side=1, 2 do begin

            outputroot = 'idlout-'+string(format='(i1,a,i4.4,a)',side, $
             '-',plateid,'-')

            expres = string(format='(a,i1,a)', 's-', side, '*.fit')
            files = findfile(filepath(expres, root_dir=plateDir), count=nfile)
            splog, 'Combining ' + strtrim(string(nfile),2) $
             + ' files for side ' + strtrim(string(side),2)

            if (nfile GT 0) then begin
               for i=0, nfile-1 do $
                splog, 'Combine file ', files[i]
               combine2dout, files, filepath(outputroot, root_dir=combineDir), $
                wavemin=alog10(3750.0), display=display, window=100
            endif

      endfor

      splog, 'Finished combining sequence', seqid[iseq], ' in', $
       systime(1)-startcombtime, ' seconds'
   endif

   endfor ; End loop for sequence number

   veryend = systime(1)
   splog, 'Total time for elapsed ', veryend-verybegin,' seconds'

   ; Close log files
   if (keyword_set(logfile)) then splog, /close
   if (keyword_set(plotfile)) then begin
      device, /close
      set_plot, 'x'
   endif

   return
end
;------------------------------------------------------------------------------
