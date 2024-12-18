;+
; NAME:
;   spreduce2d
;
; PURPOSE:
;   Calling script for SPREDUCE that reduces a night of data according
;   to a plan file.
;
; CALLING SEQUENCE:
;   spreduce2d, [ planfile, docams=, /do_telluric, /xdisplay, $
;    /writeflatmodel, /writearcmodel, /bbspec ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name(s) of output plan file; default to reducing all
;                plan files matching 'spPlan2d*.par'
;   docams     - Cameras to reduce; default to ['b1', 'b2', 'r1', 'r2']
;   do_telluric- Passed to EXTRACT_OBJECT
;   xdisplay   - Send plots to X display rather than to plot file
;   writeflatmodel - passed to SPCALIB via SPREDUCE to trigger writing
;                    out of flat model info to file.
;   writearcmodel  - passed to SPCALIB via SPREDUCE to trigger writing
;                    out of arc model info to file.
;   bbspec         - use bbspec extraction code
;   noreject       - Override Bad calibration rejection (use with caution)
;
; Optional Keywords:
;   MWM_fluxer  - Utilize MWM optional settings (ie gaia reddening and different S/N cuts)
;
;
; OUTPUT:
;
; COMMENTS:
;   The following environment variables must be set:
;      BOSS_SPECTRO_DATA
;      SDSSCORE
;      SPECFLAT_DIR
;   Look for raw FITS data files in BOSS_SPECTRO_DATA/MJD.
;   Look for obsSummary files in SDSSCORE/MJD.
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
;      Apr-2010  Added "write[flat,arc]model" pass-through (A. Bolton, Utah)
;   15-Aug-2011  Added pass-through for spatial split of sky model (A. Bolton, Utah)
;   15-Nov-2018: Modified for use only one spectrograph for the BHM (HJIM)
;-
;------------------------------------------------------------------------------

pro spreduce2d, planfile, docams=docams, do_telluric=do_telluric, saveraw=saveraw,$
 xdisplay=xdisplay, writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, $
 bbspec=bbspec, nitersky=nitersky, lco=lco, plates=plates, legacy=legacy, $
 gaiaext=gaiaext, MWM_fluxer=MWM_fluxer, map3d = map3d, force_arc2trace=force_arc2trace,$
 corrline=corrline,  no_db=no_db, debug=debug, noreject=noreject
 
 RESOLVE_ALL, /QUIET, /SKIP_EXISTING, /CONTINUE_ON_ERROR
 CPU, TPOOL_NTHREADS = 1

   if (NOT keyword_set(planfile)) then planfile = findfile('spPlan2d*.par')
   if (NOT keyword_set(nitersky)) then nitersky = 2 

   ;----------
   ; If multiple plan files exist, then call this script recursively
   ; for each such plan file.

   if (N_elements(planfile) GT 1) then begin
      for i=0, N_elements(planfile)-1 do $
       spreduce2d, planfile[i], docams=docams, do_telluric=do_telluric, $
        xdisplay=xdisplay, writeflatmodel=writeflatmodel, noreject=noreject, $
        writearcmodel=writearcmodel, bbspec=bbspec, nitersky=nitersky, $
        plates=plates, legacy=legacy, corr_line=corr_line,MWM_fluxer=MWM_fluxer, $
        no_db=no_db, saveraw=saveraw, debug=debug, force_arc2trace=force_arc2trace
      return
   endif

   thismem = memory()
   maxmem = 0

   ;----------
   ; Strip path from plan file name, and change to that directory

   thisplan = fileandpath(planfile[0], path=thispath)
   cd, thispath, current=origdir
   if (NOT keyword_set(thispath)) then cd, origdir

   ;----------
   ; Find the SPEXP structure

   allseq = yanny_readone(thisplan, 'SPEXP', hdr=hdr, /anon)
   if (N_elements(allseq) EQ 0) then begin
      splog, 'ABORT: No SPEXP structures in plan file ' + thisplan
      cd, origdir
      return
   endif

   ;----------
   ; Find keywords from the header
    
   fieldid = yanny_par(hdr,'fieldname')
   run2d = strtrim(string(yanny_par(hdr,'RUN2D')),2)
   mjd = long(yanny_par(hdr, 'MJD'))
   obs = strtrim(yanny_par(hdr,'OBS'),2)
   mjdstr = string(mjd, format='(i05.5)')

   get_field_type, fieldid=fieldid, mjd=mjd, legacy=legacy, plates=plates, fps=fps

   ;----------
   ; Read environment variables for BOSS_SPECTRO_DATA, SDSSCORE, SPECFLAT_DIR
   rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   splog, getenv('IDLSPEC2D_DIR')
   if (NOT keyword_set(rawdata_dir)) then $
        message, 'Must set environment variable BOSS_SPECTRO_DATA'
   if keyword_set(legacy) or keyword_set(plates) then begin
        speclog_dir = getenv('SPECLOG_DIR')
        if (NOT keyword_set(speclog_dir)) then $
            message, 'Must set environment variable SPECLOG_DIR'
        splog, 'Setting SPECLOG_DIR=', speclog_dir
        plugdir = concat_dir(speclog_dir, mjdstr)
   endif else begin
        sdsscore_dir = getenv('SDSSCORE_DIR')
        if (NOT keyword_set(sdsscore_dir)) then $
            message, 'Must set environment variable SDSSCORE_DIR'
        sdsscore_dir  = concat_dir(sdsscore_dir, strlowcase(obs))
        sdsscore_dir  = concat_dir(sdsscore_dir, 'summary_files')
        plugdir = sdsscore_dir
   endelse
   specflat_dir = getenv('SPECFLAT_DIR')
   if (NOT keyword_set(specflat_dir)) then $
        message, 'Must set environment variable SPECFLAT_DIR'

   ;; HJIM -- Change the default number of cameras
   if keyword_set(legacy) then begin
      if (NOT keyword_set(docams)) then docams = ['b1', 'r1', 'b2', 'r2']
      camnames = ['b1', 'r1', 'b2', 'r2']
   endif else begin
      if strmatch(obs, 'apo',/fold_case) eq 1 then begin
        if (NOT keyword_set(docams)) then docams = ['b1', 'r1']
        camnames = ['b1', 'r1']
      endif else begin
        if (NOT keyword_set(docams)) then docams = ['b2', 'r2']
        camnames = ['b2', 'r2']
      endelse
   endelse
   ncam = N_elements(camnames)


   inputdir = concat_dir(rawdata_dir, mjdstr)


   fieldmjd = field_to_string(yanny_par(hdr,'fieldname')) + '-' + mjdstr
   logfile = 'spDiag2d-' + fieldmjd + '.log'
   plotfile = 'spDiag2d-' + fieldmjd + '.ps'
   
   foo = fileandpath(planfile, path=outdir)
   if (keyword_set(outdir)) then $
    spawn, 'mkdir -p ' + outdir

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
      device, filename=plotfile, /color, /portrait
      splog, 'Plot file ' + plotfile
   endif
   splog, 'IDL version: ' + string(!version,format='(99(a," "))')
   spawn, 'uname -a', uname
   splog, 'UNAME: ' + uname[0]
   splog, 'DISPLAY=' + getenv('DISPLAY')

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()
   spawn, 'specflat_version', flatvers, /noshell
   splog, 'specflat version ' + flatvers[0]
   spawn, 'speclog_version', slogvers, /noshell
   splog, 'speclog version ' + slogvers[0]
   splog, 'Plan file ' + thisplan
   splog, 'DOCAMS = ', docams


   if keyword_set(legacy) or keyword_set(plates) then begin
     ;----------
     ; Find all the unique plate plugging names
     allnames = allseq[ sort(allseq.mapname) ].mapname
     allnames = allnames[ uniq(allnames) ]
     for imap=0, N_elements(allnames)-1 do begin
       ;----------
       ; Get the field ID number from any (e.g., the first) exposure with
       ; this sequence ID number
       thismap = allnames[imap]
       j = where(allseq.mapname EQ thismap)
       fieldid = allseq[j[0]].fieldid
       fieldstr = field_to_string(fieldid) ;- JEB plate number problem
       stime1 = systime(1)
       splog, 'Begin field ' + fieldstr + ' at ' + systime()
       ;----------
       ; Find the corresponding plug map file
       plugfile = 'plPlugMapM-' + thismap + '.par'
       splog, 'Plug map file = ', plugfile
       for ido=0, n_elements(docams)-1 do begin
         icam = (where(camnames EQ docams[ido], camct))[0]
         splog, prelog=camnames[icam]
         if (camct NE 1) then message, 'Non-unique camera ID: ' + docams[ido]
         ;----------
         ; Set the flag for splitting the sky model between spatial CCD halves:
         ; (Re: ticket #1388: strange r2 amplifier-boundary break)
         ; JEB : Strange behaviour of r1 after changing it on summer 2014
         if ((camnames[icam] eq 'r2') and (mjd ge 55300)) or $
           ((camnames[icam] eq 'r1') and (mjd ge 56858)) then $
           splitsky = 1B else splitsky = 0B
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
               + ', FIELDID= ' + fieldstr + ', CAMERA= ' + camnames[icam]
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
               + ', FIELDID= ' + fieldstr + ', CAMERA= ' + camnames[icam]
           endelse
           ;----------
           ; Get full name of pixel flat
           stime2 = systime(1)
           ;----------
           ; Reduce this set of frames (all objects w/same field + camera)
           if (keyword_set(arcname) AND keyword_set(flatname)) then begin
             plottitle = 'FIELDID='+fieldstr $
               + ' MJD='+strtrim(string(mjd),2)+' '
             spreduce, flatname, arcname, objname, run2d=run2d, $
                    plugfile=plugfile, lampfile=lampfile, noreject=noreject, $
                    indir=inputdir, plugdir=plugdir, outdir=outdir, $
                    plottitle=plottitle, do_telluric=do_telluric, $
                    writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, $
                    bbspec=bbspec, splitsky=splitsky, nitersky=nitersky, $
                    plates=plates, legacy=legacy,saveraw=saveraw,$
                    gaiaext=gaiaext,map3d = map3d, MWM_fluxer=MWM_fluxer, $
                    corrline=corrline,no_db=no_db, debug=debug
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
   endif else begin
     ;----------
     ; Find all the unique fields names
     allfields = allseq[ sort(allseq.fieldid) ].fieldid
     allfields = allfields[ uniq(allfields) ]
     for imap=0, N_elements(allfields)-1 do begin
        ;----------
        ; Get the fieldstr ID number from any (e.g., the first) exposure with
        ; this sequence ID number
        thisfield = allfields[imap]
        j = where(allseq.fieldid EQ thisfield)
        fieldstr = thisfield ; HJIM: fieldid is a string
        stime1 = systime(1)
        splog, 'Begin BHM field ' + fieldstr + ' at ' + systime()
        for ido=0, n_elements(docams)-1 do begin
           icam = (where(camnames EQ docams[ido], camct))[0]
           splog, prelog=camnames[icam]
           if (camct NE 1) then message, 'Non-unique camera ID: ' + docams[ido]
           j = where(allseq.fieldid EQ thisfield $;HJIM change mapname for field id
                 AND (allseq.flavor EQ 'science' OR allseq.flavor EQ 'smear') $
                 AND allseq.name[icam] NE 'UNKNOWN' )
           if (j[0] NE -1) then begin
              ; String array with all science exposures at this sequence + camera
              objname = allseq[j].name[icam]
              objmap = allseq[j].mapname
              for obmap = 0, n_elements(objname) -1 do begin
                if obmap EQ 0 then begin
                  confile = (findfile(filepath('confSummaryF-' + objmap[obmap] + '.par',$
                                                root_dir=plugdir, subdir='*/*'), count=ct))[0]
                  print,filepath('confSummaryF-' + objmap[obmap] + '.par',root_dir=plugdir, subdir='*/*')
                  if ct ne 0 then objobssfile = 'confSummaryF-' + objmap[obmap] + '.par' $
                             else objobssfile = 'confSummary-' + objmap[obmap] + '.par'
                  splog, 'confSummary file = ', objobssfile
                endif
                if obmap GT 0 then begin
                  confile = (findfile(filepath('confSummaryF-' + objmap[obmap] + '.par',$
                                                root_dir=plugdir, subdir='*/*'), count=ct))[0]
                  if ct ne 0 then objobssfile1 = 'confSummaryF-' + objmap[obmap] + '.par' $
                             else objobssfile1 = 'confSummary-' + objmap[obmap] + '.par'
                  splog, 'confSummary file = ', objobssfile1
                  objobssfile=[objobssfile,  objobssfile1]
                endif
              endfor
              ;objobssfile = objobssfile[ uniq(objobssfile) ]

              ;-----------
              ; Select **all** flat exposures at this sequence + camera
              j = where(allseq.fieldid EQ thisfield $;HJIM change mapname for field id
                    AND allseq.flavor EQ 'flat' $
                    AND allseq.name[icam] NE 'UNKNOWN', nflat )
              if (nflat GT 0) then begin
                 flatname = allseq[j].name[icam]
                 flatmap = allseq[j].mapname;[icam]
                 ;----------
                 ; Find the corresponding plug map file
                 for fmap = 0, nflat -1 do begin
                    if fmap EQ 0 then begin
                       confile = (findfile(filepath('confSummaryF-' + flatmap[fmap] + '.par',$
                                                    root_dir=plugdir, subdir='*/*'), count=ct))[0]
                       if ct ne 0 then calobjobssfile = 'confSummaryF-' + flatmap[fmap] + '.par' $
                                  else calobjobssfile = 'confSummary-' + flatmap[fmap] + '.par'
                       splog, 'confSummary file = ', calobjobssfile
                    endif
                    if fmap GT 0 then begin
                       confile = (findfile(filepath('confSummaryF-' + flatmap[fmap] + '.par',$
                                                    root_dir=plugdir, subdir='*/*'), count=ct))[0]
                       if ct ne 0 then calobjobssfile1 = 'confSummaryF-' + flatmap[fmap] + '.par' $
                                  else calobjobssfile1 = 'confSummary-' + flatmap[fmap] + '.par'                                  
                       splog, 'confSummary file = ', calobjobssfile1
                       calobjobssfile=[calobjobssfile,  calobjobssfile1]
                    endif
                 endfor
              endif else begin
                   flatname=''
                   splog, ' No flat for FIELDID= ' + thisfield $
                    + ', CAMERA= ' + camnames[icam];
                   splog, ' Using BOSSFLATSLIB'
                 
              endelse
              ;calobjobssfile = calobjobssfile[ uniq(calobjobssfile) ]

              ;-----------
              ; Select **all** arc exposures at this sequence + camera
              j = where(allseq.fieldid EQ thisfield $
                    AND allseq.flavor EQ 'arc' $
                    AND allseq.name[icam] NE 'UNKNOWN', narc )
              if (narc GT 0) then begin
                 arcname = allseq[j].name[icam]
              endif else begin
                 splog, ' No arc for FIELDID= ' + thisfield $
                  + ', CAMERA= ' + camnames[icam]; + ' using ' + arcname
              endelse
              ;----------
              ; Get full name of pixel flat
              stime2 = systime(1)
              ;----------
              ; Reduce this set of frames (all objects w/same plate + camera)
              if (keyword_set(arcname) AND keyword_set(flatname)) then begin
                 plottitle = ' FIELDID='+fieldstr+' ' + ' MJD='+strtrim(string(mjd),2)+' '
                 spreduce, flatname, arcname, objname, run2d=run2d, $
                        plugfile=objobssfile,calobjobssfile=calobjobssfile, $
                        lampfile=lampfile,saveraw=saveraw, noreject=noreject,$
                        indir=inputdir, plugdir=plugdir, outdir=outdir, $
                        plottitle=plottitle, do_telluric=do_telluric, $
                        writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, $
                        bbspec=bbspec, splitsky=splitsky, nitersky=nitersky,$
                        gaiaext=gaiaext, map3d = map3d, MWM_fluxer=MWM_fluxer,$
                        corrline=corrline, no_db=no_db, debug=debug,$
                        force_arc2trace=force_arc2trace
              endif
              splog, 'Time to reduce camera ', camnames[icam], ' = ', $
               systime(1)-stime2, ' seconds', format='(a,a,a,f6.0,a)'
              heap_gc   ; garbage collection
           endif
           splog, prelog=''
        endfor ; End loop for camera number
        splog, 'Time to reduce all cameras = ', $
         systime(1)-stime1, ' seconds', format='(a,f6.0,a)'
     endfor ; End loop for fps config name
   endelse

   if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then begin
      device, /close
      ps2pdf, plotfile
      set_plot, 'x'
   endif

   ; Track memory usage
   thismem = memory()
   maxmem = maxmem > thismem[3]
   splog, 'Max memory usage = ', string(maxmem/1e6,format='(f7.1)'), ' MB'

   splog, 'Total time for SPREDUCE2D = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPREDUCE2D at ' + systime()

   ;----------
   ; Close log files
   if (keyword_set(logfile)) then splog, /close

   ; Change back to original directory
   cd, origdir
   return
end
;------------------------------------------------------------------------------
