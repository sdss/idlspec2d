;+
; NAME:
;   spcombine_v5
;
; PURPOSE:
;   Calling script for SPCOADD_FRAMES.
;
; CALLING SEQUENCE:
;   spcombine_v5, [ planfile, docams=, adderr=, /xdisplay, minsn2= ]
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
;   dfpsclose
;   dfpsplot
;   idlspec2d_version()
;   idlutils_version()
;   spcoadd_frames
;   spflux
;   splog
;   yanny_free
;   yanny_par()
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   06-Jul-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro spcombine_v5, planfile, docams=docams, adderr=adderr, xdisplay=xdisplay, $
 minsn2=minsn2

   if (NOT keyword_set(planfile)) then planfile = findfile('spPlancomb*.par')
   if (n_elements(adderr) EQ 0) then adderr = 0.03
   if (n_elements(minsn2) EQ 0) then minsn2 = 0.2

   ;----------
   ; If multiple plan files exist, then call this script recursively
   ; for each such plan file.

   if (N_elements(planfile) GT 1) then begin
      for i=0, N_elements(planfile)-1 do $
       spcombine, planfile[i], docams=docams, adderr=adderr, $
        xdisplay=xdisplay, minsn=minsn
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
      set_plot, 'ps'
      dfpsplot, djs_filepath(plotfile, root_dir=combinedir), /color
      splog, 'Plot file ' + plotfile
   endif
   splog, 'Plan file ', thisplan
   splog, 'DOCAMS = ', docams

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()

   camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   ;----------
   ; Select frames that match the cameras specified by DOCAM.

   for ido=0, n_elements(docams)-1 do begin
      ii = (where(camnames EQ docams[ido], camct))[0]
      if (camct NE 1) then message, 'Non-unique camera ID: ' + docams[ido]
      if (ido EQ 0) then icams = ii $
       else icams = [icams,ii]
   endfor

   ;----------
   ; Compute a score for each frame and each exposure.
   ; Replace all UNKNOWN file names with nulls.
   ; The score will be zero if the file name is set to "NULL" or does not exist.

   dims = size(allseq)
   nexp = n_elements(allseq)
   ndocam = n_elements(icams)
   score = fltarr(ndocam, nexp)
   camspecid = lonarr(ndocam, nexp)
   expnum = lonarr(ndocam, nexp)
   for i=0L, nexp-1 do begin
      for j=0L, ndocam-1 do begin
         if (allseq[i].name[icams[j]] EQ 'UNKNOWN') then begin
            allseq[i].name[icams[j]] = ''
         endif else begin
            thisfile = (lookforgzip(djs_filepath(allseq[i].name[icams[j]], $
             root_dir=extractdir)))[0]
            if (keyword_set(thisfile)) then begin
               hdr = headfits(allseq[i].name[j])
               score[j,i] = sxpar(hdr, 'FRAMESN2')
               cameras = strtrim(sxpar(hdr, 'CAMERAS'),2)
               camspecid[j,i] = strmid(cameras, 1, 1)
               expnum[j,i] = sxpar(hdr, 'EXPOSURE')
            endif else begin
               allseq[i].name[icams[j]] = ''
            endelse
         endelse
      endfor
   endfor

   ; Discard the smear exposures by setting their scores equal to zero
   for iexp=0L, nexp-1 do $
    score[*,iexp] = score[*,iexp] * (allseq[iexp].flavor NE 'smear')

; WHAT ABOUT A MINIMUM S/N CUT, OR MINIMUM RELATIVE TO THE MAX !!!???

   ;----------
   ; Select the "best" exposure based upon the minimum score in all cameras

   expscore = fltarr(nexp)
   for iexp=0L, nexp-1 do $
    expscore[iexp] = min([score[*,iexp]])
   bestscore = max(expscore, ibest)
   splog, 'Best exposure = ', expnum[0,ibest], ' score = ', bestscore

   ;----------
   ; Compute the spectro-photometry

   i1 = where(camspecid EQ 1 AND score GT 0, ct1)
   i2 = where(camspecid EQ 2 AND score GT 0, ct2)
   objname = allseq.name[icams]

   if (ct1 GT 0) then $
    spflux_v5, objname[i1], adderr=adderr, combinedir=combinedir
   if (ct2 GT 0) then $
    spflux_v5, objname[i2], adderr=adderr, combinedir=combinedir

   ;----------
   ; Compute the flux-correction vectors

   if (ct1 GT 0) then $
    spfluxcorr_v5, objname[i1], adderr=adderr, combinedir=combinedir, $
     bestexpnum=expnum[0,ibest]
   if (ct2 GT 0) then $
    spfluxcorr_v5, objname[i2], adderr=adderr, combinedir=combinedir, $
     bestexpnum=expnum[0,ibest]

   ;----------
   ; Close plot file - S/N plots are then put in the PLOTSNFILE file.

   if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then dfpsclose

   ;----------
   ; Co-add the fluxed exposures

   ii = where(score GT 0)
   spcoadd_v5, objname[ii], combinefile, mjd=thismjd, combinedir=combinedir, $
    adderr=adderr, docams=docams, plotsnfile=plotsnfile

   heap_gc   ; garbage collection

   splog, 'Total time for SPCOMBINE = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPCOMBINE at ' + systime()

   ;----------
   ; Close log files and change to original directory

   if (keyword_set(logfile)) then splog, /close
   cd, origdir

   return
end
;------------------------------------------------------------------------------
