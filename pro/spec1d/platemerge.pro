;+
; NAME:
;   platemerge
;
; PURPOSE:
;   Merge all Spectro-1D outputs with photoPosPlate,spInspect files
;
; CALLING SEQUENCE:
;   platemerge, [ plate=, mjd=, except_tags=, indir=, outroot=, $
;    run2d=, /include_bad, /calc_noqso, /skip_line ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   plate       - Plates to include; default to all files
;                 specified by the PLATELIST routine.
;   mjd         - Optional MJDs corresponding to the specified PLATEs;
;                 if specified, then PLATE and MJD should have the same
;                 number of elements.
;   except_tags - Tag names to exclude; default to '*COVAR'.
;   indir       - Input directory with platelist.fits file; passed to
;                 platelist topir option which defaults to $BOSS_SPECTRO_REDUX
;   outroot     - Root name for output files; default to
;                 $BOSS_SPECTRO_REDUX/$RUN2D/spAll; the files are then
;                 spAll-$RUN2D.fits, spAll-$RUN2D.dat, spAllLine-$RUN2D.dat.
;   run2d       - List of RUN2D subdirectories to merge, one set of output
;                 files per name in $RUN2D; default to all values of RUN2D
;                 returned by PLATELIST.
;   include_bad  - If set, then include bad plates
;   calc_noqso  - If set, then also include redshift info for best non-QSO
;                 redshift fits.  Defaults to being set.
;   skip_line    - If set, skip the generation of spAllLine.fits
;
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Depends upon the platelist.fits file written by PLATELIST.
;   Trims to only 'good' plates, or those in a public data release.
;
;   The SPECPRIMARY output element is used to select a unique set of
;   objects in the case of duplicate observations.  Any objects observed
;   multiple times will have SPECPRIMARY=1 for one instance only, and =0
;   for all other instances.  The criteria (in order of importance) are
;   as follows:
;     1) Prefer observations with positive SN_MEDIAN in r-band
;     2) Prefer PLATEQUALITY='good' over any other plate quality
;     3) Prefer observations with ZWARNING=0
;     4) Prefer objects with larger SN_MEDIAN in r-band
;
;   Temporary files are created first, such as 'spAll.fits.tmp', which
;   are renamed at the end of the routine to 'spAll.fits', etc.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   copy_struct
;   copy_struct_inx
;   djs_filepath()
;   headfits()
;   hogg_mrdfits()
;   mrdfits()
;   mwrfits_chunks
;   platelist
;   readspec
;   repstr
;   spheregroup
;   splog
;   struct_print
;   sxaddpar
;   sxpar()
;
; REVISION HISTORY:
;   30-Oct-2000  Written by D. Schlegel, Princeton
;   29-Jul-2010  Added EXCLUDE_CLASS and SKIP_LINE, A. Bolton, Utah
;   17-Mar-2011  Changed EXCLUDE_CLASS behavior, A. Bolton, Utah
;   30-Jun-2011  Changed EXCLUDE_CLASS to more specific and correct
;                CALC_NOQSO, including proper rchi2diff, A. Bolton, Utah
;   04-Oct-2012  Added SPECBOSS tag (ASB, Utah).
;------------------------------------------------------------------------------

pro platemerge1, plate=plate, mjd=mjd, except_tags=except_tags1, $
 indir=indir, outroot=outroot1, run2d=run2d, include_bad=include_bad, $
 calc_noqso=calc_noqso, skip_line=skip_line

   dtheta = 2.0 / 3600.

   if (n_elements(except_tags1) GT 0) then except_tags = except_tags1 $
    else except_tags = '*COVAR'
   if (keyword_set(outroot1)) then begin
      outroot = [outroot1, outroot1+'Line']
   endif else begin
      outroot = ['spAll','spAllLine']
      if (keyword_set(run2d)) then outroot = outroot + '-' + run2d
      outroot = djs_filepath(outroot, root_dir=getenv('BOSS_SPECTRO_REDUX'), $
       subdir=run2d)
   endelse
   if (n_elements(calc_noqso) eq 0) then calc_noqso = 1B

   t1 = systime(1)
   thismem = memory()

   ;----------
   ; Read platelist
      
   platelist, plist=plist, topdir=indir, run2d=run2d
   if (NOT keyword_set(plist)) then return
      
   ;----------
   ; Find out if this plist includes dereddened SN2 values
   
   plist_tags = tag_names(plist)
   ii = where(plist_tags EQ 'DEREDSN2', n)
   if n GT 0 then has_deredsn2 = 1 else has_deredsn2 = 0

   ;----------
   ; Trim to good (or non-bad public) plates

   if (keyword_set(plate)) then begin
      nplate = n_elements(plist)
      if (keyword_set(mjd) AND n_elements(mjd) NE nplate) then $
       message, 'Number of elements in PLATE and MJD must agree'

      qkeep = bytarr(nplate)
      if (keyword_set(mjd)) then begin
         for i=0L, n_elements(plate)-1 do begin
            for j=0L, n_elements(mjd)-1 do begin
               qkeep = qkeep OR (plist.plate EQ plate[i] AND plist.mjd EQ mjd[j])
            endfor
         endfor
      endif else begin         
         for i=0L, n_elements(plate)-1 do $
          qkeep = qkeep OR plist.plate EQ plate[i]
      endelse
      ikeep = where(qkeep, nkeep)
      if (nkeep EQ 0) then return
      plist = plist[ikeep]
   endif

   if (keyword_set(run2d)) then begin
      qkeep = strmatch(strtrim(plist.run2d),run2d)
      ikeep = where(qkeep, nkeep)
      if (nkeep EQ 0) then return
      plist = plist[ikeep]
   endif

   qdone = strtrim(plist.status1d,2) EQ 'Done'
   if (NOT keyword_set(include_bad)) then begin
      qdone = qdone AND $
       (strtrim(plist.platequality,2) EQ 'good' $
       OR strtrim(plist.platequality,2) EQ 'marginal' $
       OR (strtrim(plist.public,2) NE '' AND $
           strtrim(plist.platequality,2) NE 'bad') )
   endif
   indx = where(qdone, ct)
   if (ct EQ 0) then return
   plist = plist[indx]

   nfile = n_elements(plist)
   nout = total(plist.n_total)
   splog, 'Total number of objects = ', nout

   ;----------
   ; Find the first tsObj file that exists for use in constructing the
   ; output structure.

   ifile = 0
   while (NOT keyword_set(tsobj0)) do begin
      readspec, plist[ifile].plate, mjd=plist[ifile].mjd, $
       run2d=strtrim(plist[ifile].run2d), tsobj=tsobj0, /silent
      if (keyword_set(tsobj0)) then begin
         tsobj0 = tsobj0[0]
      endif else begin
         ifile = ifile + 1
         if (ifile EQ nfile) then $
          message, 'No photoPosPlate files found!'
      endelse
   endwhile

   ;----------
   ; Create the additional tags to add to the output structure

   pstuff = create_struct( $
    'programname' , ' ', $
    'chunk'       , ' ', $
    'platequality', ' ', $
    'platesn2'    , 0.0, $
    'deredsn2'    , 0.0, $
    'primtarget'  ,  0L, $
    'sectarget'   ,  0L, $
    'lambda_eff'  , 0.0, $
    'bluefiber'   ,  0L, $
    'zoffset'     , 0.0, $
    'xfocal'      , 0.0, $
    'yfocal'      , 0.0, $
    'boss_target1',  0LL, $
    'boss_target2',  0LL, $
    'ancillary_target1',  0LL, $
    'ancillary_target2',  0LL, $
    'specprimary' ,  0B, $
    'specboss' ,  0B, $
    'boss_specobj_id'  ,  0L, $
    'nspecobs'    ,   0, $
;;- SB Oct 2012: remove QSO VAC inputs for DR10
;;    'z_person'    , 0.0, $
;;    'class_person',  0L, $
;;    'z_conf_person', 0L, $
;;    'comments_person', '', $
    'calibflux'   , fltarr(5), $
    'calibflux_ivar', fltarr(5) )

   ;----------
   ; Loop through each file

   splog, 'Reading ZANS files'
   for ifile=0L, nfile-1 do begin
      print, 'Reading ZANS file ',ifile+1, ' of ', nfile

      readspec, plist[ifile].plate, mjd=plist[ifile].mjd, $
       run2d=strtrim(plist[ifile].run2d), run1d=strtrim(plist[ifile].run1d), $
       zans=zans, $  ;; zmanual=zmanual, 
       plugmap=plugmap, /silent
      
      zans = struct_selecttags(zans, except_tags='OBJID')

; ASB 2011 Mar: append info on best non-galaxy and non-qso redshifts/classes:
; ASB 2011 Jun: changed to do the "no-qso" case exclusively, and more correctly.
; ASB 2011 Jul: moved this into SPREDUCE1D.  Retained for now, but
;               test to see if the tags are already present.  To be
;               removed following verification.
;;      if keyword_set(calc_noqso) then begin
      if (keyword_set(calc_noqso) and (tag_exist(zans, 'z_noqso') eq 0)) then begin
         print, '  Finding non-QSO redshift info.'
         pstring = string(plist[ifile].plate, format='(i4.4)')
         mstring = string(plist[ifile].mjd, format='(i5.5)')
         zallfile = getenv('BOSS_SPECTRO_REDUX') + '/' + $
                    strtrim(plist[ifile].run2d, 2) + '/' + $
                    pstring + '/' + strtrim(plist[ifile].run1d, 2) + $
                    '/spZall-' + pstring + '-' + mstring + '.fits'
         zall = mrdfits(zallfile,1)
         nfib = max(zall.fiberid) - min(zall.fiberid) + 1L
         nzall = n_elements(zall) / nfib
         zall = reform(zall, nzall, nfib)
         class_all = strtrim(zall.class,2)
         id_noqso = replicate(-1L, nfib)
         rchi2diff_noqso = replicate(0., nfib)
         for ii = 0L, nfib-1 do begin
            wh_noqso = where(class_all[*,ii] ne 'QSO')
            wh_noqso = (wh_noqso[sort(wh_noqso)])[0:1]
            id_noqso[ii] = wh_noqso[0]
            rchi2diff_noqso[ii] = total(zall[wh_noqso[0]:wh_noqso[1]-1,ii].rchi2diff)
         endfor
         zans_noqso = zall[id_noqso,lindgen(nfib)]
         noqso_struc = replicate( $
                       {z_noqso: 0., z_err_noqso: 0., zwarning_noqso: 0L, $
                        class_noqso: ' ', subclass_noqso: ' ', $
                        rchi2diff_noqso: 0.}, nfib)
         noqso_struc.z_noqso = zans_noqso.z
         noqso_struc.z_err_noqso = zans_noqso.z_err
         noqso_struc.class_noqso = zans_noqso.class
         noqso_struc.subclass_noqso = zans_noqso.subclass
         noqso_struc.rchi2diff_noqso = rchi2diff_noqso
; Re-set the small-delta-chi2 bit:
         minrchi2diff = 0.01
         small_rchi2diff = rchi2diff_noqso lt minrchi2diff
         zw_new = zans_noqso.zwarning
         zflagval = sdss_flagval('ZWARNING', 'SMALL_DELTA_CHI2')
         zw_new = zw_new - (zw_new and zflagval)
         zw_new = zw_new or (zflagval * small_rchi2diff)
         noqso_struc.zwarning_noqso = zw_new
         zans = struct_addtags(zans, noqso_struc)
         noqso_struc = 0
         zall = 0
         zans_noqso = 0
      endif

      if (ifile EQ 0) then begin
         outdat1 = create_struct(pstuff, zans[0])
         struct_assign, {junk:0}, outdat1 ; Zero-out all elements
         outdat = replicate(outdat1, nout)
      endif

      indx = lindgen(plist[ifile].n_total)
      if (ifile GT 0) then indx += total(plist[0:ifile-1].n_total)

      ; The following is very slow, so we do this differently...
;      copy_struct_inx, zans, outdat, index_to=indx
      tmpdat = outdat[indx]
      copy_struct, zans, tmpdat
      outdat[indx] = tmpdat

      ; Fill in the first columns of this output structure
      outdat[indx].programname = plist[ifile].programname
      outdat[indx].chunk = plist[ifile].chunk
      outdat[indx].platequality = plist[ifile].platequality
      outdat[indx].platesn2 = plist[ifile].platesn2
      if has_deredsn2 then $
       outdat[indx].deredsn2 = plist[ifile].deredsn2

      ; Read the following from the manual inspection
      ;- SB Oct 2012: removed for DR10
      ;; if (keyword_set(zmanual[0])) then begin
      ;;    outdat[indx].z_person = zmanual.z_person
      ;;    outdat[indx].class_person = zmanual.class_person
      ;;    outdat[indx].z_conf_person = zmanual.z_conf_person
      ;;    outdat[indx].comments_person = zmanual.comments
      ;; endif

      ; Get PRIMTARGET+SECTARGET with those values from
      ; the plug-map structure in spPlate file.
      outdat[indx].primtarget = plugmap.primtarget
      outdat[indx].sectarget = plugmap.sectarget
      outdat[indx].lambda_eff = plugmap.lambda_eff
      if (tag_exist(plugmap,'zoffset')) then $
       outdat[indx].zoffset = plugmap.zoffset
      if (tag_exist(plugmap,'xfocal')) then $
       outdat[indx].xfocal = plugmap.xfocal
      if (tag_exist(plugmap,'yfocal')) then $
       outdat[indx].yfocal = plugmap.yfocal
      if (tag_exist(plugmap,'bluefiber')) then $
       outdat[indx].bluefiber = plugmap.bluefiber
      if (tag_exist(plugmap,'boss_target1')) then $
       outdat[indx].boss_target1 = plugmap.boss_target1
      if (tag_exist(plugmap,'boss_target2')) then $
       outdat[indx].boss_target2 = plugmap.boss_target2
      if (tag_exist(plugmap,'ancillary_target1')) then $
       outdat[indx].ancillary_target1 = plugmap.ancillary_target1
      if (tag_exist(plugmap,'ancillary_target2')) then $
       outdat[indx].ancillary_target2 = plugmap.ancillary_target2

      ; Read the following from the plug-map if those tags exist
      if (tag_exist(plugmap,'CALIBFLUX')) then $
       outdat[indx].calibflux = plugmap.calibflux
      if (tag_exist(plugmap,'CALIBFLUX_IVAR')) then $
       outdat[indx].calibflux_ivar = plugmap.calibflux_ivar
   endfor

   splog, 'Time to read data = ', systime(1)-t1, ' sec'

   ;----------
   ; Set the SPECPRIMARY flag to 0 or 1

   t2 = systime(1)

   ; Determine the score for each object
   ; 1) Prefer observations with positive SN_MEDIAN in r-band
   ; 2) Prefer PLATEQUALITY='good' over any other plate quality
   ; 3) Prefer observations with ZWARNING=0
   ; 4) Prefer objects with larger SN_MEDIAN in r-band
; ASBjuly2011: test against ZWARNING_NOQSO for GALAXY targets:
   zw_primtest = outdat.zwarning
   if tag_exist(outdat, 'ZWARNING_NOQSO') then begin
      wh_galtarget = where(strmatch(outdat.objtype, 'GALAXY*'), ngaltarget)
      if (ngaltarget gt 0) then zw_primtest[wh_galtarget] = outdat[wh_galtarget].zwarning_noqso
   endif
   if (n_elements(outdat[0].sn_median) EQ 1) then jfilt = 0 $
    else jfilt = 2
   score = 4 * (outdat.sn_median[jfilt] GT 0) $
    + 2 * (strmatch(outdat.platequality,'good*') EQ 1) $
;;;    + 1 * (outdat.zwarning EQ 0) $ ; replaced with line below ASBjuly2011
    + 1 * (zw_primtest EQ 0) $
    + (outdat.sn_median[jfilt]>0) / max(outdat.sn_median[jfilt]+1.)

   ingroup = spheregroup(outdat.plug_ra, outdat.plug_dec, dtheta, $
   multgroup=multgroup, firstgroup=firstgroup, nextgroup=nextgroup)

   ; Set the unique object IDs
   outdat.boss_specobj_id = ingroup + 1L

   for j=0L, n_elements(firstgroup)-1L do begin
      if (firstgroup[j] NE -1) then begin
         if (multgroup[j] EQ 1) then begin
            outdat[firstgroup[j]].specprimary = 1
            outdat[firstgroup[j]].nspecobs = 1
         endif else begin
            indx = lonarr(multgroup[j])
            indx[0] = firstgroup[j]
            for k=0L, multgroup[j]-2L do indx[k+1] = nextgroup[indx[k]]
            foo = max(score[indx], ibest)
            outdat[indx[ibest]].specprimary = 1
            outdat[indx].nspecobs = multgroup[j]
         endelse
      endif
   endfor

   ; ASB: Copy specprimary into specboss
   ; (Thinking is that specprimary can be superseded downstream.)
   outdat.specboss = outdat.specprimary

   splog, 'Time to assign primaries = ', systime(1)-t2, ' sec'

   ;----------
   ; Pre-condition to FITS structure to have same-length strings
   ; (for any given tag name) by concatenating spaces.

   ntag = n_tags(outdat)
   tags = tag_names(outdat)
   for itag=0L, ntag-1L do begin
      if (size(outdat[0].(itag), /tname) EQ 'STRING') then begin
         if (NOT keyword_set(silent)) then $
          print, 'Padding whitespace for string array ' + tags[itag]
         taglen = strlen(strtrim(outdat.(itag)))
         maxlen = max(taglen)
         padspace = string('', format='(a'+string(maxlen)+')')
         outdat.(itag) = strmid(outdat.(itag) + padspace, 0, maxlen)
      endif
   endfor

   ;----------
   ; Write the output FITS file, writing one plate at a time

   ; Don't allow duplicate tags between the tsObj structure and what
   ; is already in the output structure.  For ex, MJD is in both.
   tsobj0 = struct_selecttags(tsobj0, except_tags=tag_names(outdat))
   platedat1 = create_struct(outdat[0], tsobj0)
   if (keyword_set(except_tags)) then $
    platedat1 = struct_selecttags(platedat1, except_tags=except_tags)
   struct_assign, {junk:0}, platedat1 ; Zero-out all elements

   splog, 'Writing FITS file ' + outroot[0]+'.fits'
   for ifile=0L, nfile-1 do begin
      print, 'Writing plate ', ifile+1, ' of ', nfile

      platedat = replicate(platedat1, plist[ifile].n_total)
      indx = lindgen(plist[ifile].n_total)
      if (ifile GT 0) then indx += total(plist[0:ifile-1].n_total)
      readspec, plist[ifile].plate, mjd=plist[ifile].mjd, $
       run2d=strtrim(plist[ifile].run2d), tsobj=tsobj, /silent
      if (keyword_set(tsobj)) then $
       copy_struct, tsobj, platedat $
      else $
       splog, 'WARNING: No tsObj file found for plate ', outdat[indx[0]].plate
      copy_struct, outdat[indx], platedat

      ; All strings must be the same length, or appending to the FITS file
      ; will result in corruption.  The only string in the tsobj structure
      ; is for the RERUN.
      platedat.rerun = string(platedat.rerun+'   ',format='(a3)')

      mwrfits_chunks, platedat, outroot[0]+'.fits.tmp', $
       create=(ifile EQ 0), append=(ifile GT 0)
   endfor

   outdat = 0 ; Clear memory

   ;----------
   ; Create the structure for ASCII output

   adat1 = create_struct( $
    'plate'      ,  0L, $
    'mjd'        ,  0L, $
    'fiberid'    ,  0L, $
    'class'      ,  '', $
    'subclass'   ,  '', $
    'z'          , 0.0, $
    'z_err'      , 0.0, $
    'zwarning'   ,  0L, $
    'plug_ra'    , 0.0d, $
    'plug_dec'   , 0.0d, $
    'specprimary',  0L, $
    'chunk'      ,  '', $
    'platesn2'   ,  0.0, $
    'deredsn2'   ,  0.0, $
    'objtype'    ,  '', $
    'boss_target1', 0LL, $
    'ancillary_target1', 0LL, $
    'tileid'     ,  0L, $
    'objc_type'  ,  '', $
    'modelflux'  ,  fltarr(5) )
    ;; 'z_person'   , 0.0, $
    ;; 'class_person', 0L, $
    ;; 'z_conf_person', 0L )

   tag_alias = [['SPECPRIMARY','PRIMARY'], $
    ['FIBERID','FIBER'], $
    ['BOSS_TARGET1','BOSS1'], $
    ['ANCILLARY_TARGET1','ANCILLARY1']]

   ; Read the tags that we need from the FITS file
   outdat = hogg_mrdfits(outroot[0]+'.fits.tmp', 1, nrowchunk=10000L, $
    columns=tag_names(adat1))
   adat = replicate(adat1, n_elements(outdat))
   copy_struct, outdat, adat

   ; Replace any blank strings for CLASS with "".
   ii = where(strtrim(adat.class,2) EQ '')
   if (ii[0] NE -1) then adat[ii].class = '""'

   ; Replace any blank strings for SUBCLASS with "".
   ; If SUBCLASS contains several words, then use a plus sign between
   ; the words rather than a space.
   adat.subclass = strtrim(adat.subclass,2)
   ii = where(adat.subclass EQ '')
   if (ii[0] NE -1) then adat[ii].subclass = '""'
   adat.subclass = repstr(adat.subclass, ' ', '+')

   objtypes = ['UNKNOWN', 'CR', 'DEFECT', 'GALAXY', 'GHOST', 'KNOWNOBJ', $
    'STAR', 'TRAIL', 'SKY']
   adat.objc_type = objtypes[outdat.objc_type]

   outdat = 0 ; Clear memory

   splog, 'Writing ASCII file ' + outroot[0]+'.dat'
   struct_print, adat, filename=outroot[0]+'.dat.tmp', alias=tag_alias

   adat = 0 ; Clear memory

   ;----------
   ; Create the merged line data

   if (not keyword_set(skip_line)) then begin
       splog, 'Writing FITS zline file ' + outroot[1]+'.fits'
       for ifile=0L, nfile-1 do begin
           splog, 'Writing zline ', ifile+1, ' of ', nfile
           readspec, plist[ifile].plate, mjd=plist[ifile].mjd, $
             run2d=strtrim(plist[ifile].run2d), run1d=strtrim(plist[ifile].run1d), $
             zline=linedat, /silent

           if (ifile EQ 0) then begin
               nobj = total(plist.n_total)
               nper = n_elements(linedat) / plist[0].n_total
               sxaddpar, linehdr, 'DIMS0', nper, ' Number of emission lines'
               sxaddpar, linehdr, 'DIMS1', nobj, ' Number of objects'
               linedat1 = linedat[0]
               struct_assign, {junk:0}, linedat1
           endif

      ; Demand that the structure has the same format as the first
      ; one written.
           linedat_out = replicate(linedat1, n_elements(linedat))
           struct_assign, linedat, linedat_out

           mwrfits_chunks, linedat_out, outroot[1]+'.fits.tmp', linehdr, $
            create=(ifile EQ 0), append=(ifile GT 0)
       endfor
   endif

   ;----------
   ; Rename temporary files

   spawn, ['mv', outroot[0]+'.fits.tmp', outroot[0]+'.fits'], /noshell
   spawn, ['gzip', outroot[0]+'.dat.tmp'], /noshell
   spawn, ['mv', outroot[0]+'.dat.tmp.gz', outroot[0]+'.dat.gz'], /noshell
   if (not keyword_set(skip_line)) then $
    spawn, ['mv', outroot[1]+'.fits.tmp', outroot[1]+'.fits'], /noshell

   thismem = memory()
   maxmem = thismem[3]
   splog, 'Maximum memory usage = ', maxmem/1.d6, ' MB'
   splog, 'Total time = ', systime(1)-t1, ' sec'

   return
end

;------------------------------------------------------------------------------
pro platemerge, run2d=run2d, indir=indir, _EXTRA=Extra

   platelist, plist=plist, topdir=indir, run2d=run2d
   
   if (NOT keyword_set(plist)) then return

   alldir = strtrim(plist.run2d)
   alldir = alldir[uniq(alldir, sort(alldir))]
   if (keyword_set(run2d)) then begin
      nmatch = lonarr(n_elements(alldir))
      for i=0, n_elements(run2d)-1 do $
       nmatch += strmatch(alldir,run2d[i])
      indx = where(nmatch GT 0, ct)
      if (ct EQ 0) then return
      alldir = alldir[indx]
   endif

   for i=0, n_elements(alldir)-1 do $
    platemerge1, run2d=alldir[i], indir=indir, _EXTRA=Extra

   return
end
;------------------------------------------------------------------------------
