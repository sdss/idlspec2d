;+
; NAME:
;   spreduce2d
;
; PURPOSE:
;   Calling script for SPREDUCE that reduces a night of data according
;   to a plan file.
;
; CALLING SEQUENCE:
;   spreduce2d, [ planfile, docams=, /do_telluric, /xdisplay ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name(s) of output plan file; default to reducing all
;                plan files matching 'spPlan2d*.par'
;   docams     - Cameras to reduce; default to ['b1', 'b2', 'r1', 'r2']
;   do_telluric- Passed to EXTRACT_OBJECT
;   xdisplay   - Send plots to X display rather than to plot file
;
; OUTPUT:
;
; COMMENTS:
;   The following environment variables must be set:
;      RAWDATA_DIR
;      SPECLOG_DIR
;      SPECFLAT_DIR
;   Look for raw FITS data files in RAWDATA_DIR/MJD.
;   Look for plug map files in SPECLOG_DIR/MJD.
;   Look for spectroscopic flat files in SPECFLAT_DIR.
;
; EXAMPLES:
;
; BUGS:
;   This routine spawns the Unix command 'mkdir'.
;
; PROCEDURES CALLED:
;   cpbackup
;   idlspec2d_version()
;   idlutils_version()
;   splog
;   spreduce
;   yanny_free
;   yanny_par()
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------

pro spreduce2d, planfile, docams=docams, do_telluric=do_telluric, $
 xdisplay=xdisplay

   if (NOT keyword_set(planfile)) then planfile = findfile('spPlan2d*.par')

   ;----------
   ; If multiple plan files exist, then call this script recursively
   ; for each such plan file.

   if (N_elements(planfile) GT 1) then begin
      for i=0, N_elements(planfile)-1 do $
       spreduce2d, planfile[i], docams=docams, do_telluric=do_telluric, $
        xdisplay=xdisplay
      return
   endif

   if (NOT keyword_set(docams)) then docams = ['b1', 'b2', 'r1', 'r2']

   thismem = memory()
   maxmem = 0

   ;----------
   ; Read environment variables for RAWDATA_DIR, SPECLOG_DIR, SPECFLAT_DIR

   rawdata_dir = getenv('RAWDATA_DIR')
   if (NOT keyword_set(rawdata_dir)) then $
    message, 'Must set environment variable RAWDATA_DIR'

   speclog_dir = getenv('SPECLOG_DIR')
   if (NOT keyword_set(speclog_dir)) then $
    message, 'Must set environment variable SPECLOG_DIR'

   specflat_dir = getenv('SPECFLAT_DIR')
   if (NOT keyword_set(specflat_dir)) then $
    message, 'Must set environment variable SPECFLAT_DIR'

   ;----------
   ; Strip path from plan file name, and change to that directory

   thisplan = fileandpath(planfile[0], path=thispath)
   cd, thispath, current=origdir
   if (NOT keyword_set(thispath)) then cd, origdir

   ;----------
   ; Find the SPEXP structure

   yanny_read, thisplan, pdata, hdr=hdr
   for i=0, N_elements(pdata)-1 do begin
      if (tag_names(*pdata[i], /structure_name) EQ 'SPEXP') then $
       allseq = *pdata[i]
   endfor
   yanny_free, pdata

   if (N_elements(allseq) EQ 0) then begin
      splog, 'ABORT: No SPEXP structures in plan file ' + thisplan
      cd, origdir
      return
   endif

   ;----------
   ; Find keywords from the header

   mjd = long(yanny_par(hdr, 'MJD'))
   mjdstr = string(mjd, format='(i05.5)')

   inputdir = concat_dir(rawdata_dir, mjdstr)
   plugdir = concat_dir(speclog_dir, mjdstr)
   extractdir = yanny_par(hdr, 'extractdir')
   logfile = yanny_par(hdr, 'logfile')
   plotfile = yanny_par(hdr, 'plotfile')

   if (keyword_set(extractdir)) then $
    spawn, 'mkdir -p ' + extractdir

   stime0 = systime(1)

   ;----------
   ; Open log files for output

   if (keyword_set(logfile)) then begin
      cpbackup, logfile
      splog, filename=logfile
      splog, 'Log file ' + logfile + ' opened ' + systime()
   endif
   if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then begin
      cpbackup, plotfile
      set_plot, 'ps'
      device, filename=plotfile, /color
      splog, 'Plot file ' + plotfile
   endif
   splog, 'IDL version: ' + string(!version,format='(99(a," "))')
   spawn, 'uname -a', uname
   splog, 'UNAME: ' + uname[0]
   splog, 'DISPLAY=' + getenv('DISPLAY')

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()
   spawn, 'specflat_version', flatvers
   splog, 'specflat version ' + flatvers[0]
   spawn, 'speclog_version', slogvers
   splog, 'speclog version ' + slogvers[0]

   splog, 'Plan file ' + thisplan
   splog, 'DOCAMS = ', docams

   camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   ;----------
   ; Find all the unique plate plugging names

   allnames = allseq[ sort(allseq.mapname) ].mapname
   allnames = allnames[ uniq(allnames) ]

   for imap=0, N_elements(allnames)-1 do begin

      ;----------
      ; Get the plate ID number from any (e.g., the first) exposure with
      ; this sequence ID number

      thismap = allnames[imap]
      j = where(allseq.mapname EQ thismap)
      plateid = allseq[j[0]].plateid
      platestr = string(plateid, format='(i4.4)')

      stime1 = systime(1)
      splog, 'Begin plate ' + platestr + ' at ' + systime()

      ;----------
      ; Find the corresponding plug map file

      plugfile = 'plPlugMapM-' + thismap + '.par'
      splog, 'Plug map file = ', plugfile

      for ido=0, n_elements(docams)-1 do begin

         icam = (where(camnames EQ docams[ido], camct))[0]
         splog, prelog=camnames[icam]
         if (camct NE 1) then message, 'Non-unique camera ID: ' + docams[ido]

         ;----------
         ; Find the corresponding pixel flat

         j = where(allseq.mapname EQ thismap $
               AND (allseq.flavor EQ 'science' OR allseq.flavor EQ 'smear') $
               AND allseq.name[icam] NE 'UNKNOWN' )

         if (j[0] NE -1) then begin

            ; String array with all science exposures at this sequence + camera
            objname = allseq[j].name[icam]

            ;-----------
            ; Select **all** flat exposures at this sequence + camera

            j = where(allseq.mapname EQ thismap $
                  AND allseq.flavor EQ 'flat' $
                  AND allseq.name[icam] NE 'UNKNOWN', nflat )
            if (nflat GT 0) then begin
               flatname = allseq[j].name[icam]
            endif else begin
               flatname = ''
               splog, 'ABORT: No flat for MAPNAME= ' + thismap $
                + ', PLATEID= ' + platestr + ', CAMERA= ' + camnames[icam]
            endelse

            ;-----------
            ; Select **all** arc exposures at this sequence + camera

            j = where(allseq.mapname EQ thismap $
                  AND allseq.flavor EQ 'arc' $
                  AND allseq.name[icam] NE 'UNKNOWN', narc )
            if (narc GT 0) then begin
               arcname = allseq[j].name[icam]
            endif else begin
               arcname = ''
               splog, 'ABORT: No arc for MAPNAME= ' + thismap $
                + ', PLATEID= ' + platestr + ', CAMERA= ' + camnames[icam]
            endelse

            ;----------
            ; Get full name of pixel flat

            stime2 = systime(1)

            ;----------
            ; Reduce this set of frames (all objects w/same plate + camera)

            if (keyword_set(arcname) AND keyword_set(flatname)) then begin
               plottitle = 'PLATE='+platestr $
                + ' MJD='+strtrim(string(mjd),2)+' '

               spreduce, flatname, arcname, objname, $
                plugfile=plugfile, lampfile=lampfile, $
                indir=inputdir, plugdir=plugdir, outdir=extractdir, $
                summarystruct=summarystruct, plottitle=plottitle, $
                do_telluric=do_telluric
            endif

            splog, 'Time to reduce camera ', camnames[icam], ' = ', $
             systime(1)-stime2, ' seconds', format='(a,a,a,f6.0,a)'

            heap_gc   ; garbage collection
         endif

         splog, prelog=''
      endfor ; End loop for camera number

      splog, 'Time to reduce all cameras = ', $
       systime(1)-stime1, ' seconds', format='(a,f6.0,a)'

   endfor ; End loop for plugging name

   ; Track memory usage
   thismem = memory()
   maxmem = maxmem > thismem[3]
   splog, 'Max memory usage = ', string(maxmem/1e6,format='(f7.1)'), ' MB'

   splog, 'Total time for SPREDUCE2D = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPREDUCE2D at ' + systime()

   ;----------
   ; Close log files

   if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then begin
      device, /close
      set_plot, 'x'
   endif

   if (keyword_set(logfile)) then splog, /close

   ; Change back to original directory
   cd, origdir
   return
end
;------------------------------------------------------------------------------
