;+
; NAME:
;   spcombine
;
; PURPOSE:
;   Calling script for SPCOADD_FRAMES.
;
; CALLING SEQUENCE:
;   spcombine, [ planfile, docams=, /xdisplay ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name(s) of output plan file; default to reducing all
;                plan files matching 'spPlancomb*.par'
;   docams     - Cameras to combine; default to ['b1', 'b2', 'r1', 'r2']
;   adderr     - Additional error to add to the formal errors, as a
;                fraction of the flux; default to 0.03 (3 per cent).
;   xdisplay   - Send plots to X display rather than to plot file
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   This routine spawns the Unix command 'mkdir'.
;
; PROCEDURES CALLED:
;   cpbackup
;   djs_filepath()
;   idlspec2d_version()
;   idlutils_version()
;   spcoadd_frames
;   splog
;   yanny_free
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   06-Jul-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------

pro spcombine, planfile, docams=docams, adderr=adderr, xdisplay=xdisplay, $
      noflux=noflux

   if (NOT keyword_set(planfile)) then planfile = findfile('spPlancomb*.par')
;   if (NOT keyword_set(adderr)) then adderr = 0.03

   ;----------
   ; If multiple plan files exist, then call this script recursively
   ; for each such plan file.

   if (N_elements(planfile) GT 1) then begin
      for i=0, N_elements(planfile)-1 do $
       spcombine, planfile[i], docams=docams, xdisplay=xdisplay
      return
   endif

   if (NOT keyword_set(docams)) then docams = ['b1', 'b2', 'r1', 'r2']

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

   extractdir = yanny_par(hdr, 'extractdir')
   combinedir = yanny_par(hdr, 'combinedir')
   logfile = yanny_par(hdr, 'logfile')
   plotfile = yanny_par(hdr, 'plotfile')
   fcalibprefix = yanny_par(hdr, 'fcalibprefix')
   combinefile = yanny_par(hdr, 'combinefile')

   if (keyword_set(combinedir)) then $
    spawn, 'mkdir -p ' + combinedir

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
   splog, 'Plan file ', thisplan
   splog, 'DOCAMS = ', docams

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()

   camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   ;----------
   ; Select frames that match the cameras specified by DOCAM, then trim
   ; to files that aren't names UNKNOWN, and that actually exist on disk.

   for ido=0, n_elements(docams)-1 do begin
      ii = (where(camnames EQ docams[ido], camct))[0]
      if (camct NE 1) then message, 'Non-unique camera ID: ' + docams[ido]
      if (ido EQ 0) then icams = ii $
       else icams = [icams,ii]
   endfor

   objname = allseq.name[icams]
   j = where(objname NE 'UNKNOWN')
   if (j[0] EQ -1) then begin
      splog, 'ABORT: All file names are UNKNOWN in plan file ' + thisplan
      cd, origdir
      return
   endif
   objname = objname[j]

   for ifile=0, n_elements(objname)-1 do $
    if ((findfile(objname[ifile]))[0] EQ '') then objname[ifile] = ''
   j = where(objname)
   if (j[0] EQ -1) then begin
      splog, 'ABORT: No files on disk for plan file ' + thisplan
      cd, origdir
      return
   endif
   objname = objname[j]

   ;----------
   ; Compute the spectro-photometry

   if NOT keyword_set(noflux) then $
     spflux, objname, fcalibprefix, adderr=adderr

   ;----------
   ; Co-add the fluxed exposures

   spcoadd_frames, djs_filepath(objname, root_dir=extractdir), $
    djs_filepath(combinefile, root_dir=combinedir), $
    fcalibprefix=fcalibprefix, adderr=adderr

   platesn, djs_filepath(combinefile, root_dir=combinedir), 

   heap_gc   ; garbage collection

   splog, 'Total time for SPCOMBINE = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPCOMBINE at ', systime()

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
