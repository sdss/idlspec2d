;+
; NAME:
;   newspcombine
;
; PURPOSE:
;   Flux calibrates and combines exposures with a common plugmap.
;   (Calling script for SPCOADD_FLUXED_FRAMES.)
;
; CALLING SEQUENCE:
;   newspcombine, [ planfile, docams=, adderr=, /xdisplay, minsn2=, $
;     best_exposure= ]
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
;   minsn2     - Minimum S/N^2 to include science frame in coadd (default 0.2)
;   best_exposure - string containing the exposure ID of the "best" exposure.
;                This is the exposure which all exposures will be scaled to
;                match before the images are coadded.  If this is not set, the
;                exposure with the best S/N and spectrophotometry is chosen.
;
; OUTPUT:
;   spPlate-PPPP-MMMMM.fits file
;   Diganostic plots in spDiagcomb-PPPP-MMMMM.ps, spSN2d-PPPP-MMMMM.ps
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   This routine spawns the Unix command 'mkdir'.
;   For a given exposure to be combined data from all 4 cameras must be
;   present
;
; PROCEDURES CALLED:
;   djs_filepath()
;   cpbackup
;   idlspec2d_version()
;   idlutils_version()
;   lookforgzip()
;   spcoadd_fluxed_frames
;   splog
;   yanny_free
;   yanny_par()
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   12-Aug-2003  Adapted from spcombine by C. Tremonti, Steward Observatory
;-
;------------------------------------------------------------------------------
pro newspcombine, planfile, docams=docams, adderr=adderr, xdisplay=xdisplay, $
   minsn2=minsn2, best_exposure = best_exposure

   if (NOT keyword_set(planfile)) then planfile = findfile('spPlancomb*.par')
   if (n_elements(adderr) EQ 0) then adderr = 0.03
   if (n_elements(minsn2) EQ 0) then minsn2 = 0.2

   ;----------
   ; If multiple plan files exist, then call this script recursively
   ; for each such plan file.

   if (N_elements(planfile) GT 1) then begin
      for i=0, N_elements(planfile)-1 do $
       newspcombine, planfile[i], docams=docams, xdisplay=xdisplay
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
   plotsnfile = yanny_par(hdr, 'plotsnfile')
   fcalibprefix = yanny_par(hdr, 'fcalibprefix')
   combinefile = yanny_par(hdr, 'combinefile')
   thismjd = long(yanny_par(hdr, 'MJD'))
   if (NOT keyword_set(thismjd)) then $
    thismjd = max(allseq.mjd)

   if (keyword_set(combinedir)) then $
    spawn, 'mkdir -p ' + combinedir

   stime0 = systime(1)

   ;----------
   ; Open log files for output

   if (keyword_set(logfile)) then begin
      cpbackup, djs_filepath(logfile, root_dir=combinedir)
      splog, filename=djs_filepath(logfile, root_dir=combinedir)
      splog, 'Log file ' + logfile + ' opened ' + systime()
      splog, 'IDL version: ' + string(!version,format='(99(a," "))')
      spawn, 'uname -a', uname
      splog, 'UNAME: ' + uname[0]
   endif
   if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then begin
      cpbackup, djs_filepath(plotfile, root_dir=combinedir)
      if keyword_set(plotsnfile) then $
        cpbackup, djs_filepath(plotsnfile, root_dir = combinedir)
      set_plot, 'ps'
      device, filename=djs_filepath(plotfile, root_dir=combinedir), /color, $
        /portrait, xsize=8.0, ysize=9.5, xoffset=0.25, yoffset=0.5, /inch
      splog, 'Plot file ' + plotfile
   endif
   splog, 'Plan file ', thisplan
   splog, 'DOCAMS = ', docams

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()

   ; This is the assumed ordering of the cameras in the par file (should check)
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

   ; If all file names are UNKNOWN, then abort.
   j = where(objname NE 'UNKNOWN')
   if (j[0] EQ -1) then begin
      splog, 'ABORT: All file names are UNKNOWN in plan file ' + thisplan
      cd, origdir
      return
   endif

   ; Replace all UNKNOWN file names with nulls.
   j = where(objname EQ 'UNKNOWN')
   if (j[0] NE -1) then objname[j] = ''

   ; Replace all file names that do not exist on disk with nulls.
   for ifile=0, n_elements(objname)-1 do $
    if (keyword_set(objname[ifile])) then $
     objname[ifile] = (lookforgzip(djs_filepath(objname[ifile], $
      root_dir=extractdir)))[0]

   ; Log missing files
   planname = allseq.name[icams]
   missing = where(objname eq '')

   if missing[0] ne -1 then $
      splog, 'Files not found on disk: ' + planname[missing]
 
   ; Now all the file names in ALLSEQ.NAME should exist or be set to null.
   allseq.name[icams] = objname

   ;-------------------------------
   ;  For simplicity require data for all 4 cameras -- if this is not done
   ; then some changes need to be made to "frame_flux_tweak"  
   fourcam = (allseq.name[0] ne '') + (allseq.name[2] ne '') + $
             (allseq.name[1] ne '') + (allseq.name[3] ne '')
   missing1 = where(fourcam gt 0 and fourcam lt 4)
   if missing1[0] ne -1 then allseq[missing1].name[*] = ''

   j = where(allseq.name[icams])
   if (j[0] EQ -1) then begin
      splog, 'ABORT: No files on disk for plan file ' + thisplan
      cd, origdir
      return
   endif
   objname = (allseq.name[icams])[j]

   ;----------
   ; Separate science & smear frames

   isci = where(allseq.flavor EQ 'science', nsci)
   ismear = where(allseq.flavor EQ 'smear')

   if (isci[0] EQ -1) then begin 
      splog, 'ABORT: No science frames in this plan ' + thisplan
      cd, origdir
      return
   endif

   ;----------
   ; Check for Minimum S/N in science frame  -- NEW: require min S/N for 
   ; all four frames of an exposure
   
   if keyword_set(minsn2) then begin
      lowsn2 = fltarr(nsci)

      for i=0,nsci-1 do begin
         checkhdr1 = headfits(allseq[isci[i]].name[0])
         checkhdr2 = headfits(allseq[isci[i]].name[1])
         checkhdr3 = headfits(allseq[isci[i]].name[2])
         checkhdr4 = headfits(allseq[isci[i]].name[3])
         if size(checkhdr1,/tname) NE 'INT' then begin
           lowsn2[i] = (sxpar(checkhdr1,'FRAMESN2') LT minsn2) + $
                       (sxpar(checkhdr2, 'FRAMESN2') LT minsn2) + $
                       (sxpar(checkhdr3, 'FRAMESN2') LT minsn2) + $
                       (sxpar(checkhdr4, 'FRAMESN2') LT minsn2) 
         endif
      endfor

      goodexp = where(lowsn2 eq 0, ngood)
      if ngood lt nsci then begin 
        splog, 'Excluded ', total(lowsn2), $
             ' frames with SN^2 less than ', minsn2, format='(a,i4,a,f7.2)'
      endif
      if ngood eq 0 then begin
        splog, 'ABORT: All Frames rejected due to minimum S/N limit'
        cd, origdir
        return
      endif
      isci = isci[goodexp] 
   endif

   sciname = allseq[isci].name[icams]
   j = where(sciname)
   if j[0] EQ -1 then begin 
      splog, 'ABORT: No science frames in this plan ' + thisplan
      cd, origdir
      return
   endif
   sciname = sciname[j]

  ;----------------------------
  ; Check for existance of smear frames

   if (ismear[0] EQ -1) then begin 
      splog, 'No smear frames in this plan ' + thisplan
      smearname = ''
   endif else begin
     smearname = allseq[ismear].name[icams]
     j = where(smearname) 
     if j[0] EQ -1 then begin 
       splog, 'No smear frames in this plan ' + thisplan
       smearname = ''
     endif else smearname = smearname[j] 
   endelse

   ;----------
   ; Check that all files to be combined have the same plugmap

   nsci = n_elements(sciname)
   mapname = strarr(nsci)
   for i = 0, nsci - 1 do begin
     checkhdr = headfits(sciname[i])
     mapname[i] = strtrim(sxpar(checkhdr, 'NAME'), 2) 
   endfor
   diff = where(mapname ne mapname[0], ndiff)
   if ndiff gt 0 then begin
     splog, 'ABORT: Files have different plugmaps!!'
      cd, origdir
      return
   endif
  
   ;----------
   ; Check for the presence of a new tsObj file

   tsobjname = djs_filepath('tsObj-' + mapname[0] + '.fit', $
               root_dir = extractdir)
   if file_test(tsobjname) then begin
      splog, 'Using fibermags from tsObj-' + mapname[0] 
   endif else begin
      splog, 'No new tsObj found ... using fibermags from plugmap'
      tsobjname = ''
   endelse 

   ;----------
   ; Co-add the fluxed exposures

   spcoadd_fluxed_frames, sciname, combinefile, mjd=thismjd, $
    combinedir=combinedir, fcalibprefix=fcalibprefix, adderr=adderr, $
    docams=docams, plotsnfile=plotsnfile, tsobjname = tsobjname, $
    smearname = smearname, best_exposure = best_exposure

   ;----------
   ; Close plot file - S/N plots are then put in the PLOTSNFILE file.

   if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then begin
      set_plot, 'x'
   endif

   heap_gc   ; garbage collection

   splog, 'Total time for SPCOMBINE = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPCOMBINE at ', systime()


   ;----------
   ; Close log files and change to original directory

   if (keyword_set(logfile)) then splog, /close
   cd, origdir

   return
end
;------------------------------------------------------------------------------
