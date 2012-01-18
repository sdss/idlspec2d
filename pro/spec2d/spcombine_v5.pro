;+
; NAME:
;   spcombine_v5
;
; PURPOSE:
;   Calling script for SPCOADD_V5.
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
;   minsn2     - Minimum S/N^2 to include science frame in coadd; default
;                to 0 to only include those with S/N > 0.
;                Note that all exposures with a score less than 0.2 times
;                the score of the best exposure are discarded; for those
;                purposes, the score used is the worst of all 4 cameras.
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   We currently hard-wire the rejection of all smears and
;   any with (S/N)^2 less than 20% of the best exposure.
;
; PROCEDURES CALLED:
;   cpbackup
;   dfpsclose
;   dfpsplot
;   headfits
;   idlspec2d_version()
;   idlutils_version()
;   spcoadd_v5
;   spflux_v5
;   spfluxcorr_v5
;   splog
;   sxpar()
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
  if (n_elements(minsn2) EQ 0) then minsn2 = 0.
 
  thismem = memory()
  maxmem = 0
 
  ;----------
  ; If multiple plan files exist, then call this script recursively
  ; for each such plan file.

  if (N_elements(planfile) GT 1) then begin
    for i=0, N_elements(planfile)-1 do $
      spcombine_v5, planfile[i], docams=docams, adderr=adderr, $
      xdisplay=xdisplay, minsn=minsn
    return
  endif

  if (NOT keyword_set(docams)) then docams = ['b1', 'b2', 'r1', 'r2']

  ;----------
  ; Strip path from plan file name, and change to that directory

  thisplan = fileandpath(planfile[0], path=outdir)
  cd, outdir, current=origdir
  if (NOT keyword_set(outdir)) then cd, origdir

  ;----------
  ; Find the SPEXP structure

  allseq = yanny_readone(thisplan, 'SPEXP', hdr=hdr, /anon)
  if (N_elements(allseq) EQ 0) then begin
    splog, 'ABORT: No SPEXP structures in plan file ' + thisplan
    cd, origdir
    return
  endif

  ;----------
  ; Find keywords from the header and construct output file names

  thismjd = long(yanny_par(hdr, 'MJD'))
  if (NOT keyword_set(thismjd)) then $
   thismjd = max(allseq.mjd)
  platemjd = string(yanny_par(hdr,'plateid'),format='(i4.4)') $
   + '-' + string(thismjd,format='(i5.5)')
  logfile = 'spDiagcomb-' + platemjd + '.log'
  plotfile = 'spDiagcomb-' + platemjd + '.ps'
  fcalibprefix = 'spFluxcalib-' + platemjd
  combinefile = 'spPlate-' + platemjd + '.fits'
  plotsnfile = 'spSN2d-' + platemjd + '.ps'

  stime0 = systime(1)

  ;----------
  ; Open log files for output

  if (keyword_set(logfile)) then begin
    cpbackup, djs_filepath(logfile, root_dir=outdir)
    splog, filename=djs_filepath(logfile, root_dir=outdir)
    splog, 'Log file ' + logfile + ' opened ' + systime()
    splog, 'IDL version: ' + string(!version,format='(99(a," "))')
    spawn, 'uname -a', uname
    splog, 'UNAME: ' + uname[0]
  endif
  if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then begin
    cpbackup, djs_filepath(plotfile, root_dir=outdir)
    set_plot, 'ps'
    dfpsplot, djs_filepath(plotfile, root_dir=outdir), /color
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
  ; The score will be MINSN2 if the file name is set to "NULL"
  ; or does not exist.

  dims = size(allseq)
  nexp = n_elements(allseq)
  ndocam = n_elements(icams)
  score = fltarr(ndocam, nexp) - (minsn2<0)
  camspecid = lonarr(ndocam, nexp)
  expnum = lonarr(ndocam, nexp)
  camerasarr= strarr(ndocam, nexp)
  for i=0L, nexp-1 do begin
    for j=0L, ndocam-1 do begin
      if (allseq[i].name[icams[j]] EQ 'UNKNOWN') then begin
        allseq[i].name[icams[j]] = ''
      endif else begin
        thisfile = (lookforgzip(djs_filepath(allseq[i].name[icams[j]], $
          root_dir=outdir)))[0]
        if (keyword_set(thisfile)) then begin
          hdr = headfits(thisfile)
          score[j,i] = sxpar(hdr, 'FRAMESN2')
          cameras = strtrim(sxpar(hdr, 'CAMERAS'),2)
          camspecid[j,i] = strmid(cameras, 1, 1)
          camerasarr[j,i] = cameras
          expnum[j,i] = sxpar(hdr, 'EXPOSURE')
        endif else begin
          expnum[j,i] = long(strmid(allseq[i].name[icams[j]],11,8))
          allseq[i].name[icams[j]] = ''
        endelse
      endelse
    endfor
  endfor

  ;----------
  ; If all data is missing from one of the cameras, then discard
  ; both cameras from that spectrograph.

  ; Case where we discard spectrograph #1
  ; If there are no exposures with both b1+b2 cameras
  ; then discard all spectrograph#1 data
  if (total( total(camerasarr EQ 'b1',1) GT 0 AND $
   total(camerasarr EQ 'r1',1) GT 0 ) EQ 0) then begin
     ii = where(docams EQ 'b2' OR docams EQ 'r2', ct)
     if (ct GT 0) then begin
        splog, 'Discarding spectro-1 data'
        camerasarr = camerasarr[ii,*]
        camspecid = camspecid[ii,*]
        docams = docams[ii]
        expnum = expnum[ii,*]
        icams = icams[ii]
        score = score[ii,*]
     endif
  endif

  ; Case where we discard spectrograph #2
  ; If there are no exposures with both b1+b2 cameras
  ; then discard all spectrograph#1 data
  if (total( total(camerasarr EQ 'b2',1) GT 0 AND $
   total(camerasarr EQ 'r2',1) GT 0 ) EQ 0) then begin
     ii = where(docams EQ 'b1' OR docams EQ 'r1', ct)
     if (ct GT 0) then begin
        splog, 'Discarding spectro-2 data'
        camerasarr = camerasarr[ii,*]
        camspecid = camspecid[ii,*]
        docams = docams[ii]
        expnum = expnum[ii,*]
        icams = icams[ii]
        score = score[ii,*]
     endif
  endif

  ; Discard the smear exposures by setting their scores equal to (MINSN2<0)
  qsmear = allseq.flavor EQ 'smear'
  for iexp=0L, nexp-1 do $
    score[*,iexp] = score[*,iexp] * (qsmear[iexp] EQ 0) $
    + (minsn2<0) * qsmear[iexp]

  ;----------
  ; Select the "best" exposure based upon the minimum score in all cameras

  expscore = fltarr(nexp)
  for iexp=0L, nexp-1 do $
    expscore[iexp] = min([score[*,iexp]])
  bestscore = max(expscore, ibest)
  splog, 'Best exposure = ', expnum[0,ibest], ' score = ', bestscore

  ;----------
  ; Discard exposures whose score is less than some fraction of the
  ; best exposure, or whose score is less than some absolute value.
  ; These numbers are hard-wired!!!???

  ibad = where(expscore LE minsn2 OR expscore LT 0.20*bestscore, nbad)
  if (nbad GT 0) then begin
    for j=0, nbad-1 do splog, 'WARNING: Discarding ' $
      + allseq[ibad[j]].flavor + ' exposure #', $
      expnum[0,ibad[j]], ' with score=', expscore[ibad[j]]
    score[*,ibad] = (minsn2<0)
  endif

  ;----------
  ; Compute the spectro-photometry

  i1 = where(camspecid EQ 1 AND score GT minsn2, ct1)
  i2 = where(camspecid EQ 2 AND score GT minsn2, ct2)
  objname = allseq.name[icams]

  configuration=obj_new("configuration",thismjd)
  splog, prename='sp1'
  if (ct1 GT 0) then begin
    spflux_v5, objname[i1], adderr=adderr, combinedir=outdir, $
     minfracthresh=configuration->spflux_v5_minfracthresh()
  endif
  splog, prename='sp2'
  if (ct2 GT 0) then begin
    spflux_v5, objname[i2], adderr=adderr, combinedir=outdir, $
     minfracthresh=configuration->spflux_v5_minfracthresh()
  endif
  splog, prename=''

  ; Track memory usage
  thismem = memory()
  maxmem = maxmem > thismem[3]
  splog, 'Max memory usage = ', string(maxmem/1e6,format='(f7.1)'), ' MB'

  ;----------
  ; Compute the flux-correction vectors

  if (ct1 GT 0) then $
   spfluxcorr_v5, objname[i1], adderr=adderr, combinedir=outdir, $
    bestexpnum=expnum[0,ibest]
  if (ct2 GT 0) then $
   spfluxcorr_v5, objname[i2], adderr=adderr, combinedir=outdir, $
    bestexpnum=expnum[0,ibest]

  ; Track memory usage
  thismem = memory()
  maxmem = maxmem > thismem[3]
  splog, 'Max memory usage = ', string(maxmem/1e6,format='(f7.1)'), ' MB'

  ;----------
  ; Close plot file - S/N plots are then put in the PLOTSNFILE file.

  if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then dfpsclose

  ;----------
  ; Co-add the fluxed exposures

  ii = where(score GT minsn2, ct)
  if (ct GT 0) then begin
     spcoadd_v5, objname[ii], combinefile, mjd=thismjd, combinedir=outdir, $
      adderr=adderr, docams=docams, plotsnfile=plotsnfile, $
      bestexpnum=expnum[0,ibest]
  endif else $
     splog, 'ABORT: No exposures with SCORE > ' + strtrim(string(minsn2),2)
  obj_destroy,configuration    
  heap_gc   ; garbage collection

  ; Track memory usage
  thismem = memory()
  maxmem = maxmem > thismem[3]
  splog, 'Max memory usage = ', string(maxmem/1e6,format='(f7.1)'), ' MB'

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
