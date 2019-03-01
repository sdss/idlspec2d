;+
; NAME:
;   rm_spcombine_v5
;
; PURPOSE:
;   Calling script for SPCOADD_V5.
;
; CALLING SEQUENCE:
;   rm_spcombine_v5, [ planfile, docams=, adderr=, /xdisplay, minsn2= ]
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
;   skipfluxing- Skip the step to generate spFluxcalib* files
;   nofcorr    - Skip the step to generate and use the spFluxcorr* files
;   nodist     - Skip the step to generate and use the spFluxdistort* files
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
;   18-Jan-2019  Modified by Hector Ibarra
;   16-Jan-2014  Modified by Yue Shen to test improved RM spectrophotometry
;   06-Jul-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro rm_spcombine_v5, planfile, docams=docams, adderr=adderr, xdisplay=xdisplay, $
 minsn2=minsn2, topdir=topdir,finaldir=finaldir,nprox=nprox, oneexp=oneexp, $
 skipfluxing=skipfluxing, nofcorr=nofcorr,nodist=nodist,useairmass=useairmass, $
 xyfit=xyfit, skipfcorr=skipfcorr, loaddesi=loaddesi, lco=lco

  if (NOT keyword_set(planfile)) then planfile = findfile('spPlancomb*.par')
  if (n_elements(adderr) EQ 0) then adderr = 0.03
  if (n_elements(minsn2) EQ 0) then minsn2 = 0.0
  if keyword_set(lco) then begin
    obsdir='LCO'
  endif else begin
    obsdir='APO'
  endelse
  ;----------
 
  

  ; Default is to do the spectrophotometry on all exposures each
  if n_elements(oneexp) eq 0 then oneexp = 1L

  ; Default is to use airmass in determining the local fluxing std array
  if n_elements(useairmass) eq 0 then useairmass = 1L

  thismem = memory()
  maxmem = 0
 
  ;----------
  ; If multiple plan files exist, then call this script recursively
  ; for each such plan file.

  if (N_elements(planfile) GT 1) then begin
    for i=0, N_elements(planfile)-1 do $
      rm_spcombine_v5, planfile[i], docams=docams, adderr=adderr, $
      xdisplay=xdisplay, minsn=minsn, topdir=topdir, nprox=nprox, $
      oneexp=oneexp, finaldir=finaldir, skipfluxing=skipfluxing,$
      nofcorr=nofcorr,nodist=nodist,useairmass=useairmass,xyfit=xyfit, $
      skipfcorr=skipfcorr,loaddesi=loaddesi
    return
  endif
  obsdir='';coment this line for the final version HJIM
  if not keyword_set(topdir) then $
    topdir=getenv('BOSS_SPECTRO_REDUX') + '/' + obsdir + '/' + getenv('RUN2D') $
    + '/' + strmid(planfile,11,4) + '/'


  if (NOT keyword_set(docams)) then docams = ['b1', 'r1']

  ;----------
  ; Strip path from plan file name, and change to that directory

  thisplan = fileandpath(topdir+planfile[0], path=outdir)
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
  fieldmjd = string(yanny_par(hdr,'fieldid'),format='(i4.4)') $
   + '-' + string(thismjd,format='(i5.5)')
  for i=0, n_elements(allseq.mjd)-1 do begin
    if i EQ 0 then begin
      ;logfile = 'spDiagcomb-'+fieldmjd+'-'+string(i,format='(i2.2)')+'.log'
      ;plotfile = 'spDiagcomb-'+fieldmjd+'-'+string(i,format='(i2.2)')+'.ps'
      fcalibprefix = 'spFluxcalib-'+fieldmjd+'-'+string(i,format='(i2.2)')
      ;plotsnfile = 'spSN2d-'+fieldmjd+'-'+string(i,format='(i2.2)')+'.ps'
    endif else begin
      ;logfile1 = 'spDiagcomb-'+fieldmjd+'-'+string(i,format='(i2.2)')+'.log'
      ;plotfile1 = 'spDiagcomb-'+fieldmjd+'-'+string(i,format='(i2.2)')+'.ps'
      fcalibprefix1 = 'spFluxcalib-'+fieldmjd+'-'+string(i,format='(i2.2)')
      ;plotsnfile1 = 'spSN2d-'+fieldmjd+'-'+string(i,format='(i2.2)')+'.ps'
      ;logfile=[logfile,logfile1]
      ;plotfile=[plotfile,plotfile1]
      fcalibprefix=[fcalibprefix,fcalibprefix1]
      ;plotsnfile=[plotsnfile,plotsnfile1]
    endelse
  endfor
  combinefile = 'spField-' + fieldmjd+'.fits'
  plotfile = 'spDiagcomb-'+fieldmjd+'.ps'
  ;logfile = 'spDiagcomb-'+fieldmjd+'-'+string(i,format='(i2.2)')+'.log'
  logfile = 'spDiagcomb-'+fieldmjd+'.log'
  plotsnfile = 'spSN2d-'+fieldmjd+'-X.ps'
  if keyword_set(skipfluxing) then begin
    ;for i=0, n_elements(allseq.mjd)-1 do begin
    ;  if i EQ 0 then begin
         logfile = 'spDiagcomb-'+fieldmjd+'-quick.log'
    ;  endif else begin
    ;     logfile1 = 'spDiagcomb-'+fieldmjd+'-'+string(i,format='(i2.2)')+'-quick.log'
    ;     logfile=[logfile,logfile1]
    ;  endelse
    ;endfor
  endif

  ; modify the RM recalibration combinedir
  if not keyword_set(finaldir) then outdir = outdir $;'recalib/' $
    else outdir = outdir + finaldir

  stime0 = systime(1)

  ;----------
  ; Open log files for output

  if (keyword_set(logfile)) then begin
    for i=0, n_elements(logfile)-1 do begin
      cpbackup, djs_filepath(logfile[i], root_dir=outdir)
      splog, filename=djs_filepath(logfile[i], root_dir=outdir)
      splog, 'Log file ' + logfile[i] + ' opened ' + systime()
      splog, 'IDL version: ' + string(!version,format='(99(a," "))')
      spawn, 'uname -a', uname
      splog, 'UNAME: ' + uname[0]
    endfor
  endif

  if not keyword_set(skipfluxing) then begin
    if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then begin
      for i=0, n_elements(plotfile)-1 do begin
        cpbackup, djs_filepath(plotfile[i], root_dir=outdir)
        set_plot, 'ps'
        dfpsplot, djs_filepath(plotfile[i], root_dir=outdir), /color
        splog, 'Plot file ' + plotfile[i]
      endfor
    endif
  endif

  splog, 'Plan file ', thisplan
  splog, 'DOCAMS = ', docams

  splog, 'idlspec2d version ' + idlspec2d_version()
  splog, 'idlutils version ' + idlutils_version()

  camnames = ['b1', 'r1']
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
          root_dir=topdir)))[0]
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

  ibad = where(expscore LE minsn2 OR expscore LT 0.20*bestscore, nbad, complemen =igood)
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

  ; stop here and run diagnosis
  ; message, 'Stop here and run diagonis' 


  if not keyword_set(skipfluxing) then begin
     splog, prename='sp1'
     if not keyword_set(oneexp) then begin
        splog, 'Do Fluxcalib vectors using all exposures'
        if (ct1 GT 0) then begin
           rm_spflux_v5, objname[i1], adderr=adderr, combinedir=outdir, $
            minfracthresh=configuration->spflux_v5_minfracthresh(),nprox=nprox, $
            useairmass=useairmass,bestexpnum=bestexpnum_sp1,xyfit=xyfit, $
            loaddesi=loaddesi
           bestexp_b1 = expnum[0,igood[bestexpnum_sp1[0]]]
           bestexp_r1 = expnum[1,igood[bestexpnum_sp1[1]]]
           splog, 'Best exposure for spectrophotometry (blue1): ', bestexp_b1
           splog, 'Best exposure for spectrophotometry (red1): ', bestexp_r1

        endif
     ; instead of spline all the exps, do the spectro-photometry on one exp at a time
     endif else begin
        splog, 'Do Fluxcalib vectors using individual exposures'
        if (ct1 GT 0) then begin
           for ido=0L, ct1/2 - 1L do begin
              rm_spflux_v5, objname[i1[ido*2:ido*2+1]], adderr=adderr, $
               combinedir=outdir, nprox=nprox, $
               minfracthresh=configuration->spflux_v5_minfracthresh(), $
               useairmass=useairmass,xyfit=xyfit,loaddesi=loaddesi
            endfor
        endif
     endelse


     ;splog, prename='sp2'
     ;if not keyword_set(oneexp) then begin
     ;   splog, 'Do Fluxcalib vectors using all exposures'
     ;   if (ct2 GT 0) then begin
     ;      rm_spflux_v5, objname[i2], adderr=adderr, combinedir=outdir, $
     ;       minfracthresh=configuration->spflux_v5_minfracthresh(),nprox=nprox, $
     ;       useairmass=useairmass,bestexpnum=bestexpnum_sp2,xyfit=xyfit, $
     ;       loaddesi=loaddesi
     ;      bestexp_b2 = expnum[2,igood[bestexpnum_sp2[0]]]
     ;      bestexp_r2 = expnum[3,igood[bestexpnum_sp2[1]]]
     ;      splog, 'Best exposure for spectrophotometry (blue2): ', bestexp_b2
     ;      splog, 'Best exposure for spectrophotometry (red2): ', bestexp_b2


     ;   endif
     ;endif else begin
     ;   splog, 'Do Fluxcalib vectors using individual exposures'
     ;   if (ct2 GT 0) then begin
     ;      for ido=0L, ct2/2 - 1L do $
     ;         rm_spflux_v5, objname[i2[ido*2:ido*2+1]], adderr=adderr, $
     ;          combinedir=outdir, nprox = nprox, $
     ;          minfracthresh=configuration->spflux_v5_minfracthresh(), $
     ;          useairmass=useairmass,xyfit=xyfit,loaddesi=loaddesi
     ;   endif
     ;endelse
     splog, prename=''

     ; Track memory usage
     thismem = memory()
     maxmem = maxmem > thismem[3]
     splog, 'Max memory usage = ', string(maxmem/1e6,format='(f7.1)'), ' MB'

  endif else splog, 'Skip fluxcalibration (assume already done earlier)'

  ;----------
  ; Compute the flux-correction vectors

  if not keyword_set(nofcorr) and not keyword_set(skipfcorr) then begin
     ; for the flux-correction vectors, still use all the exps at once
 
     if not keyword_set(oneexp) then begin
       splog, 'Do Flux-correction vectors for all exposures'
       if (ct1 GT 0) then $
          rm_spfluxcorr_v5, objname[i1], adderr=adderr, combinedir=outdir, $
            bestexpnum=[bestexp_b1,bestexp_r1]  ;expnum[0,ibest]
     endif else begin
       splog, 'Do Flux-correction vectors for individual exposures'
          if (ct1 GT 0) then begin
            for ido=0L, ct1/2 - 1L do begin
              rm_spfluxcorr_v5, objname[i1[ido*2:ido*2+1]], adderr=adderr, $ 
                combinedir=outdir;, indf=(ido+1);, bestexpnum=[bestexp_b1,bestexp_r1]
            endfor
          endif
     endelse
     ;if (ct2 GT 0) then $
     ; rm_spfluxcorr_v5, objname[i2], adderr=adderr, combinedir=outdir, $
     ;  bestexpnum=[bestexp_b2,bestexp_r2]  ;expnum[0,ibest]

     ; Track memory usage
     thismem = memory()
     maxmem = maxmem > thismem[3]
     splog, 'Max memory usage = ', string(maxmem/1e6,format='(f7.1)'), ' MB'

     ;----------
     ; Close plot file - S/N plots are then put in the PLOTSNFILE file.

     if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then dfpsclose
   endif else splog, 'Skip flux-corretion'

  ;exit
  ;----------
  ; Co-add the fluxed exposures

  ii = where(score GT minsn2, ct)
  if (ct GT 0) then begin
     ; set values for wavemin and wavemax so that all coadded plates/spectra are
     ; on the same wavelength grid; wavemin and wavemax are Log-10 wavelength of 
     ; the first and last pixel
     wavemin=3.5523 & wavemax=4.0171
     rm_spcoadd_v5, objname[ii], combinefile, mjd=thismjd, combinedir=outdir, $
      adderr=adderr, docams=docams, plotsnfile=plotsnfile, $
      bestexpnum=expnum[0,ibest],nofcorr=nofcorr,nodist=nodist, $
      wavemin=wavemin, wavemax=wavemax
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
