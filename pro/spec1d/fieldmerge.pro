;+
; NAME:
;   fieldmerge
;
; PURPOSE:
;   Merge all Spectro-1D outputs with photoPosPlate,spInspect files
;
; CALLING SEQUENCE:
;   fieldmerge, [ field=, mjd=, except_tags=, indir=, outroot=, $
;    run2d=, /include_bad, /calc_noqso, /skip_line ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   field       - fields to include; default to all files
;                 specified by the fieldLIST routine.
;   mjd         - Optional MJDs corresponding to the specified fields;
;                 if specified, then field and MJD should have the same
;                 number of elements.
;   except_tags - Tag names to exclude; default to '*COVAR'.
;   indir       - Input directory with fieldlist.fits file; passed to
;                 fieldlist topir option which defaults to $BOSS_SPECTRO_REDUX
;   outroot     - Root name for output files; default to
;                 $BOSS_SPECTRO_REDUX/$RUN2D/spAll; the files are then
;                 spAll-$RUN2D.fits, spAll-$RUN2D.dat, spAllLine-$RUN2D.dat.
;   run2d       - List of RUN2D subdirectories to merge, one set of output
;                 files per name in $RUN2D; default to all values of RUN2D
;                 returned by fieldLIST.
;   include_bad - If set, then include bad fields
;   calc_noqso  - If set, then also include redshift info for best non-QSO
;                 redshift fits.  Defaults to being set.
;   skip_line   - If set, skip the generation of spAllLine.fits
;   mergerun2d  - If set, generate a single $BOSS_SPECTRO_REDUX/spAll.fits
;                 file combining all RUN2D versions in
;                 $BOSS_SPECTRO_REDUX/fieldlist.fits.
;                 Ignores run2d, $RUN2D, and does *not* write a separate
;                 file for each RUN2D.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Depends upon the fieldlist.fits file written by fieldLIST.
;   Trims to only 'good' fields, or those in a public data release.
;
;   The SPECPRIMARY output element is used to select a unique set of
;   objects in the case of duplicate observations.  Any objects observed
;   multiple times will have SPECPRIMARY=1 for one instance only, and =0
;   for all other instances.  The criteria (in order of importance) are
;   as follows:
;     1) Prefer observations with positive SN_MEDIAN in r-band
;     2) Prefer fieldQUALITY='good' over any other field quality
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
;   fieldlist
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



pro fieldmerge1, field=field, mjd=mjd, except_tags1=except_tags1, $
 indir=indir, outroot1=outroot1, run2d=run2d, include_bad=include_bad, $
 calc_noqso=calc_noqso, skip_line=skip_line, plist=plist, legacy=legacy, $
 plates=plates, photo_file=photo_file, XCSAO=XCSAO, $
 lite=lite, skip_specprimary=skip_specprimary

   dtheta = 2.0 / 3600.

   if (n_elements(except_tags1) GT 0) then except_tags = except_tags1 $
    else except_tags = '*COVAR'
   if (keyword_set(outroot1)) then begin
      outroot = [outroot1, outroot1+'Line']
   endif else begin
      outroot = ['spAll','spAllLine']
      if (keyword_set(field) and keyword_set(mjd)) then begin
        outroot='spectra/full/'+field_to_string(field)+'/'+strtrim(string(mjd),2)+'/'+outroot+'-'+field_to_string(field)+'-'+strtrim(string(mjd),2)
      endif else begin
        if (keyword_set(run2d)) then outroot = outroot + '-' + repstr(run2d,'/','-')
      endelse
      outroot = djs_filepath(outroot, root_dir=getenv('BOSS_SPECTRO_REDUX'), $
       subdir=run2d)
   endelse
   if (n_elements(calc_noqso) eq 0) then calc_noqso = 1B

   t1 = systime(1)
   thismem = memory()
   splog, outroot

   ;----------
   ; Read fieldlist if needed
     
   if (NOT keyword_set(plist)) then $ 
    conflist, plist=plist, topdir=indir, run2d=run2d
   if (NOT keyword_set(plist)) then return
   
   ;----------
   ; Find out if this plist includes dereddened SN2 values
   
   plist_tags = tag_names(plist)
   ii = where(plist_tags EQ 'DEREDSN2', n)
   if n GT 0 then has_deredsn2 = 1 else has_deredsn2 = 0

   ;----------
   ; Trim to good (or non-bad public) plates

   if (keyword_set(field) and keyword_set(mjd)) then begin
      nplate = n_elements(plist)

      field_plate=plist.field
      qkeep = bytarr(nplate)
      if (keyword_set(mjd)) then begin
        qkeep = qkeep OR (field_plate EQ field AND plist.mjd EQ mjd)
      endif
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
   ;HJIM decomet this part for the final version
   if (NOT keyword_set(include_bad)) then begin
      qdone = qdone AND $
       (strtrim(plist.fieldquality,2) EQ 'good' $
       OR strtrim(plist.fieldquality,2) EQ 'marginal' $
       OR (strtrim(plist.public,2) NE '' AND $
           strtrim(plist.fieldquality,2) NE 'bad') )
   endif

   indx = where(qdone eq 1, ct)
   ;print, indx
   if (ct EQ 0) then return
   plist = plist[indx]

   nfile = n_elements(plist)
   nout = total(plist.n_total)
   splog, 'Total number of objects = ', nout

   ;----------
   ; Find the first tsObj file that exists for use in constructing the
   ; output structure.
   ifile = 0
   ;;HJIM coment the next line for the final version
   ;tsobj0=1
   brake_t=0
   no_photo_file=1
   ; If the photoPlate files are chosen, add the tsObj structure
   if keyword_set(photo_file) then begin
    while (NOT keyword_set(tsobj0)) and brake_t eq 0 do begin
      ;HJIM Coment the upper line and decoment the lower line
      readspec, plist[ifile].field, mjd=plist[ifile].mjd, $
        run2d=strtrim(plist[ifile].run2d), tsobj=tsobj0, $
        legacy=legacy, plates=plates, /silent
       tsobj0 = tsobj0[0]
       brake_t=1
       no_photo_file=0
      if (NOT keyword_set(tsobj0)) then begin
         brake_t=0
         no_photo_file=1
         ifile = ifile + 1
         if (ifile EQ nfile) then $
          brake_t=1
      endif
    endwhile
   endif

   ;----------
   ; Create the additional tags to add to the output structure

   pstuff = create_struct( $
    'programname' , ' ', $
    'survey'      , ' ', $
    'fieldquality', ' ', $
    'fieldsn2'    , 0.0, $
    'exp_disp_med', 0.D, $
    'fiberid_list', ' ', $
    'lambda_eff'  , 0.0, $
    'bluefiber'   ,  0L, $
    'zoffset'     , 0.0, $
    'xfocal'      , ' ', $
    'yfocal'      , ' ', $
    'calibflux'   , fltarr(5), $
    'calibflux_ivar', fltarr(5),$
    'gaia_bp', 0.0, $
    'gaia_rp', 0.0, $
    'gaia_g', 0.0, $
    'cadence', ' ', $
    'firstcarton', ' ', $
    'carton_to_target_pk', ' ', $
    'RACAT', 0.D,$
    'DECCAT', 0.D,$
    'COORD_EPOCH', 0.0,$
    'PMRA', 0.0,$
    'PMDEC', 0.0,$
    'PARALLAX', 0.0,$
    'catalogid'      , long64(0), $
    'catalogid_v0'   , long64(0), $
    'catalogid_v0p5' , long64(0), $
    ;'catalogid_v1'   , long64(0), $
    'gaia_id_dr2', long64(0), $
    'FIBER2MAG', fltarr(5), $
    'PSFMAG', fltarr(5), $
    'obs', '',$
    'field', 0L, $
    'plate', 0L, $
    'designs', '', $
    'configs', '', $
    'nexp', 0, $
    'exptime', 0, $
    'airmass', 0.0, $
    'Seeing20', 0.0, $
    'Seeing50', 0.0, $
    'Seeing80', 0.0, $
    'Assigned', ' ', $
    'on_target', ' ', $
    'valid', ' ', $
    'healpix', 0L, $
    'healpixgrp', 0, $
    'healpix_path', ' ', $
    'mjd_final', 0.D, $
    'mjd_list', ' ', $
    'tai_list', ' ', $
    'fieldsnr2g_list', ' ', $
    'fieldsnr2r_list', ' ', $
    'fieldsnr2i_list', ' ', $
    'RA_LIST', ' ', $
    'DEC_LIST', ' ', $
    'DELTA_RA_LIST', ' ',$
    'DELTA_DEC_LIST', ' ',$
    'moon_dist', ' ', $
    'moon_phase', ' ', $
    'ebv', 0.0, $
    'ebv_type', '', $
    'wise_mag', fltarr(4), $
    'twomass_mag', fltarr(3), $ 
    'guvcat_mag', fltarr(2), $
    ;'gaia_parallax', 0.0, $
    ;'gaia_pmra', 0.0, $
    ;'gaia_pmdec', 0.0,$
    'fiber_offset', 0.0,$
    'SPEC_FILE', ' ')

    if not keyword_set(skip_specprimary) then $
        pstuff = struct_addtags(pstuff, {specprimary:0B, specboss:0B, $
                                         boss_specobj_id:0L, nspecobs:0})

    if keyword_set(legacy) then $
        pstuff = struct_addtags(pstuff, {chunk:'',DEREDSN2:'',primtarget:0L,sectarget:0L})
    
    if keyword_set(XCSAO) then $
        pstuff = struct_addtags(pstuff, {XCSAO_rv: !values.f_nan, XCSAO_erv: !values.f_nan,$
                                        XCSAO_Rxc: !values.f_nan, $
                                        XCSAO_Teff: !values.f_nan, XCSAO_eteff: !values.f_nan,$
                                        XCSAO_Logg: !values.f_nan, XCSAO_elogg: !values.f_nan,$
                                        XCSAO_Feh: !values.f_nan, XCSAO_efeh: !values.f_nan})
   ;----------
   ; Loop through each file

   splog, 'Reading ZANS files'
   for ifile=0L, nfile-1 do begin
      splog, 'Reading ZANS file ',ifile+1, ' of ', nfile
      field_plate=field_to_string(plist[ifile].field)
     
      readspec, field_plate, mjd=plist[ifile].mjd, $
       run2d=strtrim(plist[ifile].run2d), run1d=strtrim(plist[ifile].run1d), $
       zans=zans, objhdr=objhdr, $  ;; zmanual=zmanual, 
       plugmap=plugmap, legacy=legacy, plates=plates, /silent, unsigned=(ifile EQ 0)

      zans = struct_selecttags(zans, except_tags=['OBJID','TILE'])
      if tag_exist(zans,'SPEC2_G') then zans = struct_selecttags(zans, except_tags='SPEC2_G')
      if tag_exist(zans,'SPEC2_R') then zans = struct_selecttags(zans, except_tags='SPEC2_R')
      if tag_exist(zans,'SPEC2_I') then zans = struct_selecttags(zans, except_tags='SPEC2_I')



; ASB 2011 Mar: append info on best non-galaxy and non-qso redshifts/classes:
; ASB 2011 Jun: changed to do the "no-qso" case exclusively, and more correctly.
; ASB 2011 Jul: moved this into SPREDUCE1D.  Retained for now, but
;               test to see if the tags are already present.  To be
;               removed following verification.
      if (keyword_set(calc_noqso) and (tag_exist(zans, 'z_noqso') eq 0)) then begin
         splog, '  Finding non-QSO redshift info.'

		;;- JB : Change field string format PLATEPROBLEM
         pstring = field_to_string(field_plate)
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

      if keyword_set(xcsao) then begin
        pstring = field_to_string(field_plate)
        mstring = string(plist[ifile].mjd, format='(i5.5)')
        XCSAOfile = getenv('BOSS_SPECTRO_REDUX') + '/' + $
                strtrim(plist[ifile].run2d, 2) + '/' + $
                pstring + '/' + strtrim(plist[ifile].run1d, 2) + $
                '/spXCSAO-' + pstring + '-' + mstring + '.fits'
        if  FILE_TEST(XCSAOfile) then begin
            splog, 'Reading XCSAO file: spXCSAO-' + pstring + '-' + mstring + '.fits'
            XCSAO = mrdfits(XCSAOfile,1)
        endif else splog, 'No XCSAO file for spXCSAO-' + pstring + '-' + mstring + '.fits'
      endif
      if (ifile EQ 0) then begin
         htags = ['SDSSV_BOSS_TARGET0']
         pstuff = struct_addtags(pstuff, $
          struct_selecttags(plugmap[0], select_tags=htags))
         ;pstuff = create_struct(pstuff, plutt[0])
         outdat1 = create_struct(pstuff, struct_selecttags(zans[0], except_tags=['field','fiberid_list']))
         struct_assign, {junk:0}, outdat1 ; Zero-out all elements
         if keyword_set(xcsao) then begin
            outdat1.xcsao_rv=!values.f_nan
            outdat1.xcsao_erv=!values.f_nan
            outdat1.XCSAO_Rxc=!values.f_nan
            outdat1.XCSAO_Teff=!values.f_nan
            outdat1.XCSAO_eteff=!values.f_nan
            outdat1.XCSAO_Logg=!values.f_nan
            outdat1.XCSAO_elogg=!values.f_nan
            outdat1.XCSAO_Feh=!values.f_nan
            outdat1.XCSAO_efeh=!values.f_nan
         endif
         outdat = replicate(outdat1, nout)
      endif

      indx = lindgen(plist[ifile].n_total)
      if (ifile GT 0) then indx += total(plist[0:ifile-1].n_total)

      tmpdat = outdat[indx]
      copy_struct, zans, tmpdat
      outdat[indx] = tmpdat

      ; Fill in the first columns of this output structure
;      if tag_exist(plist,'programname') then begin
;         outdat[indx].programname = plist[ifile].programname
;      endif else outdat[indx].programname=plugmap.program
 
      outdat = strct_to_struct(plugmap,'*','PROGRAM',outdat,indx, outtag='programname')
      outdat = strct_to_struct(plist,ifile,'CHUNK',outdat,indx)
      outdat = strct_to_struct(plist,ifile,'PLATEQUALITY',outdat,indx, altTag='fieldquality')
      outdat = strct_to_struct(plist,ifile,'PLATESN2',outdat,indx,outtag='fieldsn2', altTag='fieldsn2')
      if has_deredsn2 then outdat = strct_to_struct(plist,ifile,'deredsn2',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','EXP_DISP_MED', outdat,indx)
;      outdat = strct_to_struct(plist,ifile,'exptime',outdat,indx)
      outdat = strct_to_struct(plist,ifile,'DESIGNID',outdat,indx)
;      outdat = strct_to_struct(plist,ifile,'DESIGNS',outdat,indx)
;      outdat = strct_to_struct(plist,ifile,'CONFIGS',outdat,indx)
;      outdat = strct_to_struct(plist,ifile,'SURVEY',outdat,indx)
      outdat = strct_to_struct(plist,ifile,'OBS',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','SURVEY',outdat,indx)

      outdat = strct_to_struct(plugmap,'*','AIRMASS',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','SEEING20',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','SEEING50',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','SEEING80',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','MJDLIST',outdat,indx, outTag='mjd_list')
      outdat = strct_to_struct(plugmap,'*','TAILIST',outdat,indx, outTag='tai_list')
      outdat = strct_to_struct(plugmap,'*','DESIGNS',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','CONFIGS',outdat,indx)



      ; Get PRIMTARGET+SECTARGET with those values from
      ; the plug-map structure in spfield file.
      ; HJIM decoment the next three lines for the final version
      outdat = strct_to_struct(plugmap,'*','primtarget',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','sectarget',outdat,indx)

      outdat = strct_to_struct(plugmap,'*','lambda_eff',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','zoffset',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','XFOCAL_LIST',outdat,indx,outTag='xfocal')
      outdat = strct_to_struct(plugmap,'*','YFOCAL_LIST',outdat,indx,outTag='yfocal')
      outdat = strct_to_struct(plugmap,'*','bluefiber',outdat,indx)

      outdat = strct_to_struct(plugmap,'*','CALIBFLUX',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','CALIBFLUX_IVAR',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','CALIBFLUX_IVAR',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','GAIA_BP',outdat,indx,altTag='BP_MAG',altOutTag='GAIA_BP')
      outdat = strct_to_struct(plugmap,'*','GAIA_RP',outdat,indx,altTag='RP_MAG',altOutTag='GAIA_RP')
      outdat = strct_to_struct(plugmap,'*','GAIA_G',outdat,indx,altTag='GAIA_G_MAG',altOutTag='GAIA_G')
      outdat = strct_to_struct(plugmap,'*','SDSSV_BOSS_TARGET0',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','CADENCE',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','FIRSTCARTON',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','CARTON_TO_TARGET_PK',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','ASSIGNED_LIST',outdat,indx, outTag='ASSIGNED')
      outdat = strct_to_struct(plugmap,'*','ON_TARGET_LIST',outdat,indx, outTag='ON_TARGET')
      outdat = strct_to_struct(plugmap,'*','VALID_LIST',outdat,indx, outTag='VALID')
      outdat = strct_to_struct(plugmap,'*','ICATALOGID',outdat,indx, outTag='CATALOGID')
      
      catversion = plugmap.catversion
      v0 = where(strmatch(catversion, '0.0*') eq 1, ct)
      if ct gt 0 then outdat = strct_to_struct(plugmap, v0, 'ICATALOGID', outdat, indx[v0], outTag='CATALOGID_V0')
      
      v0p5 = where(strmatch(catversion, '0.5*') eq 1, ct)
      if ct gt 0 then outdat = strct_to_struct(plugmap, v0p5, 'ICATALOGID', outdat, indx[v0p5], outTag='CATALOGID_V0p5')
      
;      v1 = where(strmatch(catversion, '1.*') eq 1, ct)
;      if ct gt 0 then outdat = strct_to_struct(plugmap, v1, 'ICATALOGID', outdat, indx[v1], outTag='CATALOGID_V1')
      
      outdat = strct_to_struct(plugmap,'*','gaia_id_dr2',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','FIBER2MAG',outdat,indx,altTag='mag')
      outdat = strct_to_struct(plugmap,'*','PSFMAG',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','NEXP',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','EXPTIME',outdat,indx)
;      outdat = strct_to_struct(plugmap,'*','field',outdat,indx, altTag='plate', altOutTag='field')
      outdat = strct_to_struct(plugmap,'*','mjd_final',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','RACAT',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','RACAT',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','DECCAT',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','COORD_EPOCH',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','PMRA',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','PMDEC',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','PARALLAX',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','TAI_LIST',outdat,indx)

      outdat = strct_to_struct(plugmap,'*','fieldsnr2g_list',outdat,indx,altTag='platesnr2g_list', altOutTag='fieldsnr2i_list')
      outdat = strct_to_struct(plugmap,'*','fieldsnr2r_list',outdat,indx,altTag='platesnr2r_list', altOutTag='fieldsnr2i_list')
      outdat = strct_to_struct(plugmap,'*','fieldsnr2i_list',outdat,indx,altTag='platesnr2i_list', altOutTag='fieldsnr2i_list')

      outdat = strct_to_struct(plugmap,'*','RA_LIST',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','DEC_LIST',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','DELTA_RA_LIST',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','DELTA_DEC_LIST',outdat,indx)

      outdat = strct_to_struct(plugmap,'*','MOON_DIST',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','MOON_PHASE',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','EBV',outdat,indx,altTag='SFD_EBV', altOutTag='EBV')
      outdat = strct_to_struct(plugmap,'*','EBV_TYPE',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','WISE_MAG',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','TWOMASS_MAG',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','GUVCAT_MAG',outdat,indx)
      ;outdat = strct_to_struct(plugmap,'*','GAIA_PARALLAX',outdat,indx)
      ;outdat = strct_to_struct(plugmap,'*','GAIA_PMRA',outdat,indx)
      ;outdat = strct_to_struct(plugmap,'*','GAIA_PMDEC',outdat,indx)
      outdat = strct_to_struct(plugmap,'*','fiber_offset',outdat,indx)

    
      if keyword_set(xcsao) then begin
        if FILE_TEST(XCSAOfile) then begin
            outdat = strct_to_struct(XCSAO,'*','rv',outdat,indx,outTag='XCSAO_rv')
            outdat = strct_to_struct(XCSAO,'*','erv',outdat,indx,outTag='XCSAO_erv')
            outdat = strct_to_struct(XCSAO,'*','R',outdat,indx,outTag='XCSAO_Rxc')
            outdat = strct_to_struct(XCSAO,'*','Teff',outdat,indx,outTag='XCSAO_Teff')
            outdat = strct_to_struct(XCSAO,'*','Teff',outdat,indx,outTag='XCSAO_Teff')
            outdat = strct_to_struct(XCSAO,'*','eteff',outdat,indx,outTag='XCSAO_eteff')
            outdat = strct_to_struct(XCSAO,'*','Logg',outdat,indx,outTag='XCSAO_Logg')
            outdat = strct_to_struct(XCSAO,'*','elogg',outdat,indx,outTag='XCSAO_elogg')
            outdat = strct_to_struct(XCSAO,'*','Feh',outdat,indx,outTag='XCSAO_Feh')
            outdat = strct_to_struct(XCSAO,'*','efeh',outdat,indx,outTag='XCSAO_efeh')
        endif
      endif
      
      spec_file = outdat[indx].spec_file
      for fid = 0L, n_elements(spec_file)-1 do begin
                spec_file[fid] = 'spec-'+ field_to_string((outdat[indx].field)[fid]) + '-' + $
                                 strtrim(string(plist[ifile].mjd),2) + '-' +$
                                 strtrim(plugmap[fid].catalogid,2)+'.fits'
      endfor
      outdat[indx].spec_file = spec_file
      
      ; Add the cas-styled specobjid to output
      if not tag_exist(outdat, 'specobjid') then $
        outdat = struct_addtags(outdat, replicate({specobjid:0ULL},n_elements(outdat)))
      words= STREGEX(STRTRIM(outdat[indx].run2d,2),'^v([0-9]+)_([0-9]+)_([0-9]+)', /SUB, /EXTRACT)
      ; did it parse as vXX_YY_ZZ?
      if words[0] ne '' then begin
        rerun= (long(words[1,*])-5L)*10000L+ (long(words[2,*])*100L)+ (long(words[3,*]))
      endif else begin
        splog, "WARNING: Unable to parse RERUN from ", (outdat[indx])[uniq(outdat[indx].run2d)].run2d, " for CAS-style SPECOBJID; Using 0 instead"
        rerun= intarr(n_elements(outdat[indx].field))
      endelse
      outdat[indx].specobjid = sdss_specobjid_17(outdat[indx].field,outdat[indx].target_index,$
                                                    outdat[indx].mjd,rerun)
      
      
      healpix_now=1
      if keyword_set(healpix_now) then begin
        mwm_root='$MWM_HEALPIX'
        healpix_path_t=outdat[indx].healpix_path
        plt_t=field_to_string(outdat[indx].field)
        healp=coords_to_healpix(zans.fiber_ra,zans.fiber_dec)
        outdat[indx].healpix=healp.healpix
        outdat[indx].healpixgrp=healp.healpixgrp

        for fid = 0L, n_elements(zans)-1 do begin
            if plugmap[fid].icatalogid  eq 0 then healpix_path_t[fid] = '' $
            else begin
                healpix_path_t[fid]=mwm_root + '/'+ $
                  strtrim(string(healp[fid].healpixgrp),2) + '/' + $
                  strtrim(string(healp[fid].healpix),2) + '/boss/' + $
                  strtrim(plist[ifile].run2d, 2)+ '/' + $
                  'spec-' + strtrim(string(plt_t[0]),2) + '-' + $
                  strtrim(string(plist[ifile].mjd),2) + $
                  '-' + strtrim(plugmap[fid].catalogid,2)+'.fits'
            endelse
        endfor
        outdat[indx].healpix_path=healpix_path_t
      endif else begin
      
        outdat = strct_to_struct(plugmap,'*','HEALPIX',outdat,indx)
        outdat = strct_to_struct(plugmap,'*','HEALPIXGRP',outdat,indx)
        if (tag_exist(plugmap,'HEALPIX_PATH') or tag_exist(plugmap,'HEALPIX_DIR')) then begin
          if tag_exist(plugmap,'HEALPIX_DIR') then healpix_path_t=plugmap.healpix_dir $
          else healpix_path_t=plugmap.healpix_path
          for fid = 0L, n_elements(zans)-1 do begin
            healpix_path_t[fid] = repstr(healpix_path_t[fid],'XXXX',strtrim(string(plist[ifile].mjd),2))
          endfor
          outdat[indx].healpix_path=healpix_path_t
        endif
      endelse
   endfor
      
   splog, 'Time to read data = ', systime(1)-t1, ' sec'

   ;----------
   ; Set the SPECPRIMARY flag to 0 or 1
   if not keyword_set(skip_specprimary) then begin
       t2 = systime(1)

       ; Determine the score for each object
       ; 1) Prefer observations with positive SN_MEDIAN in r-band
       ; 2) Prefer fieldQUALITY='good' over any other field quality
       ; 3) Prefer observations with ZWARNING=0
       ; 4) Prefer objects with larger SN_MEDIAN in r-band
     ; ASBjuly2011: test against ZWARNING_NOQSO for GALAXY targets:
       zw_primtest = outdat.zwarning
       if tag_exist(outdat, 'ZWARNING_NOQSO') then begin
          wh_galtarget = where(strmatch(outdat.objtype, 'GALAXY*'), ngaltarget)
          if (ngaltarget gt 0) then zw_primtest[wh_galtarget] = outdat[wh_galtarget].zwarning_noqso
       endif
       if (n_elements(outdat[0].sn_median) EQ 1) then jfilt = 0 else jfilt = 2
       score = 4 * (outdat.sn_median[jfilt] GT 0) $
             + 2 * (strmatch(outdat.fieldquality,'good*') EQ 1) $
             + 1 * (zw_primtest EQ 0) $
             + (outdat.sn_median[jfilt]>0) / max(outdat.sn_median[jfilt]+1.)

       ingroup = spheregroup(outdat.fiber_ra, outdat.fiber_dec, dtheta, $
                             multgroup=multgroup, firstgroup=firstgroup, $
                             nextgroup=nextgroup)

       ; Set the unique object IDs
       if tag_exist(outdat, 'boss_specobj_id') then outdat.boss_specobj_id = ingroup + 1L

       if tag_exist(outdat, 'specprimary') then sp = 1 else sp = 0
       if tag_exist(outdat, 'nspecobs') then ns =1 else ns = 0
       for j=0L, n_elements(firstgroup)-1L do begin
          if (firstgroup[j] NE -1) then begin
             if (multgroup[j] EQ 1) then begin
                if sp then outdat[firstgroup[j]].specprimary = 1
                if ns then outdat[firstgroup[j]].nspecobs = 1
             endif else begin
                indx = lonarr(multgroup[j])
                indx[0] = firstgroup[j]
                for k=0L, multgroup[j]-2L do indx[k+1] = nextgroup[indx[k]]
                foo = max(score[indx], ibest)
                if sp then outdat[indx[ibest]].specprimary = 1
                if ns then outdat[indx].nspecobs = multgroup[j]
             endelse
          endif
       endfor

   
       ; ASB: Copy specprimary into specboss
       ; (Thinking is that specprimary can be superseded downstream.)
       if tag_exist(outdat, 'specboss') then outdat.specboss = outdat.specprimary

       splog, 'Time to assign primaries = ', systime(1)-t2, ' sec'
   endif 

   ;----------
   ; Pre-condition to FITS structure to have same-length strings
   ; (for any given tag name) by concatenating spaces.

   ntag = n_tags(outdat)
   tags = tag_names(outdat)
   for itag=0L, ntag-1L do begin
      if (size(outdat[0].(itag), /tname) EQ 'STRING') then begin
         if (NOT keyword_set(silent)) then $
          splog, 'Padding whitespace for string array ' + tags[itag]
         taglen = strlen(strtrim(outdat.(itag)))
         maxlen = max(taglen)
         padspace = string('', format='(a'+string(maxlen)+')')
         outdat.(itag) = strmid(outdat.(itag) + padspace, 0, maxlen)
      endif
   endfor

   ;----------
   ; Write the output FITS file, writing one field at a time

   ; Don't allow duplicate tags between the tsObj structure and what
   ; is already in the output structure.  For ex, MJD is in both.
   if no_photo_file eq 0 then begin
     tsobj0 = struct_selecttags(tsobj0, except_tags=tag_names(outdat));;HJIM Decoment this part
     platedat1 = create_struct(outdat[0], tsobj0)
   endif else begin
     platedat1 = outdat[0]
   endelse
   if (keyword_set(except_tags)) then $
    platedat1 = struct_selecttags(platedat1, except_tags=except_tags)
   struct_assign, {junk:0}, platedat1 ; Zero-out all elements

   splog, 'Writing FITS file ' + outroot[0]+'.fits'
   for ifile=0L, nfile-1 do begin
      field_plate=plist[ifile].field
      splog, 'Writing field ', ifile+1, ' of ', nfile
      str_p='field'

      platedat = replicate(platedat1, plist[ifile].n_total)
      indx = lindgen(plist[ifile].n_total)
      if (ifile GT 0) then indx += total(plist[0:ifile-1].n_total)
      if keyword_set(photo_file) then begin
        readspec, field_plate, mjd=plist[ifile].mjd, $
         run2d=strtrim(plist[ifile].run2d), tsobj=tsobj, $
         legacy=legacy, plates=plates, /silent
        if (keyword_set(tsobj)) then $
         copy_struct, tsobj, platedat $
        else $
         splog, 'WARNING: No tsObj file found for '+str_p+' ', outdat[indx[0]].field
      endif else begin
        readspec, field_plate, mjd=plist[ifile].mjd, $
         run2d=strtrim(plist[ifile].run2d), $
         legacy=legacy, plates=plates, /silent
      endelse
      copy_struct, outdat[indx], platedat

      ; All strings must be the same length, or appending to the FITS file
      ; will result in corruption.  The only string in the tsobj structure
      ; is for the RERUN.
      ; HJIM Decoment the next line
      if no_photo_file eq 0 then $
        platedat.rerun = string(platedat.rerun+'   ',format='(a3)')

      mwrfits_chunks, platedat, outroot[0]+'.fits.tmp', $
       create=(ifile EQ 0), append=(ifile GT 0)
   endfor

   outdat = 0 ; Clear memory

   ;----------
   ; Create the structure for ASCII output field

   adat1 = create_struct( $
    str_p        ,  0L, $
    'mjd'        ,  0L, $
    'fiberid'    ,  0L, $
    'class'      ,  '', $
    'subclass'   ,  '', $
    'z'          , 0.0, $
    'z_err'      , 0.0, $
    'zwarning'   ,  0L, $
    'fiber_ra'    , 0.0d, $
    'fiber_dec'   , 0.0d, $
    'specprimary',  0L, $
    'fieldsn2'   ,  0.0, $
    'objtype'    ,  '', $
    'tileid'     ,  0L, $
    'objc_type'  ,  '', $
    'modelflux'  ,  fltarr(5) )


    if keyword_set(legacy) then $
        adat1 = struct_addtags(adat1, {chunk:'',DEREDSN2:'',$
                                       specprimary:0B,specboss:0B,boss_specobj_id:0L,$
                                       nspecobs:0})

   tag_alias = [['SPECPRIMARY','PRIMARY'], $
    ['FIBERID','FIBER']];, $

   ; Read the tags that we need from the FITS file
   outdat = hogg_mrdfits(outroot[0]+'.fits.tmp', 1, nrowchunk=10000L, $
    columns=tag_names(adat1), /unsigned)
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
   ;HJIM decomet the next line
   if no_photo_file eq 0 then $
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
           field_plate=plist[ifile].field
           splog, 'Writing zline ', ifile+1, ' of ', nfile
           readspec, field_plate, mjd=plist[ifile].mjd, $
             run2d=strtrim(plist[ifile].run2d), run1d=strtrim(plist[ifile].run1d), $
             zline=linedat, legacy=legacy, plates=plates, /silent             
           if not (tag_exist(linedat,'field')) then begin
              adatag = create_struct('field',field_plate)
           endif  
           if (ifile EQ 0) then begin
               nobj = total(plist.n_total)
               nper = n_elements(linedat) / plist[0].n_total
               sxaddpar, linehdr, 'DIMS0', nper, ' Number of emission lines'
               sxaddpar, linehdr, 'DIMS1', nobj, ' Number of objects'
               linedat1 = linedat[0]
               struct_assign, {junk:0}, linedat1
               if keyword_set(adatag) then begin
                 linedat1=create_struct(adatag,linedat1)
               endif
           endif
        ; Demand that the structure has the same format as the first
        ; one written.
           linedat_out = replicate(linedat1, n_elements(linedat))
           ;adatag_in = replicate(adatag, n_elements(linedat))
           struct_assign, linedat, linedat_out
           ;struct_assign, adatag_in, linedat_out
           if keyword_set(adatag) then begin
              linedat_out.field=field_plate
           endif
           mwrfits_chunks, linedat_out, outroot[1]+'.fits.tmp', linehdr, $
            create=(ifile EQ 0), append=(ifile GT 0)
       endfor
   endif

   if keyword_set(lite) then begin
        exclude_spall_tags=['FIELDID_LIST','XFOCAL','YFOCAL','DESIGNS','CONFIGS',$
                            'MJD_LIST','TAI_LIST','FIELDSN2G_LIST','FIELDSN2R_LIST',$
                            'FIELDSN2I_LIST','RA_LIST', 'DEC_LIST', 'BLUEFIBER',$
                            'DELTA_RA_LIST', 'DELTA_DEC_LIST',$
                            'ZOFFSET','SDSSV_BOSS_TARGET0','ASSIGNED',$
                            'ON_TARGET','VALID','MOON_DIST','MOON_PHASE','FIBERID_LIST',$
                            'FIELDSNR2G_LIST','FIELDSNR2R_LIST','FIELDSNR2I_LIST',$
                            'TCOLUMN','NPOLY','THETA','FRACNSIGMA','FRACNSIGHI',$
                            'FRACNSIGLO','HEALPIX_PATH','carton_to_target_pk']
        
        spall = mrdfits(outroot[0]+'.fits.tmp',1)
        spall_lite = struct_selecttags(spall,except_tags=exclude_spall_tags)
        spall_lite = struct_addtags(spall_lite, replicate(create_struct('ASSIGNED',0S,$
                                                          'ON_TARGET',0S,$
                                                          'VALID',0S,$
                                                          'MOON_DIST',0.0,$
                                                          'MOON_PHASE',0.0,$
                                                          'carton_to_target_pk',0LL),n_elements(spall)))
        for i=0, n_elements(spall)-1 do begin
            spall_lite[i].ASSIGNED   = min(FIX(strsplit(spall[i].ASSIGNED,/extract)))
            spall_lite[i].ON_TARGET  = min(FIX(strsplit(spall[i].ON_TARGET,/extract)))
            spall_lite[i].VALID      = min(FIX(strsplit(spall[i].VALID,/extract)))
            spall_lite[i].MOON_DIST  = avg(float(strsplit(spall[i].MOON_DIST,/extract)))
            spall_lite[i].MOON_PHASE = avg(float(strsplit(spall[i].MOON_PHASE,/extract)))
            spall_lite[i].carton_to_target_pk = LONG64(spall[i].carton_to_target_pk)
        endfor
        spall_name=REPSTR(outroot[0],'spAll','spAll-lite')
        MWRFITS, spall_lite, spall_name+'.fits.tmp', fits_hdr, Status=Status
        spawn, ['gzip', spall_name+'.fits.tmp'], /noshell
   endif

   ;----------
   ; Rename temporary files
   ;print, outroot[0]
   spawn, ['mv', outroot[0]+'.fits.tmp', outroot[0]+'.fits'], /noshell
   spawn, ['gzip', outroot[0]+'.dat.tmp'], /noshell
   spawn, ['mv', outroot[0]+'.dat.tmp.gz', outroot[0]+'.dat.gz'], /noshell
   if (not keyword_set(skip_line)) then $
    spawn, ['mv', outroot[1]+'.fits.tmp', outroot[1]+'.fits'], /noshell
   if keyword_set(lite) then $
    spawn, ['mv',   spall_name+'.fits.tmp.gz', spall_name+'.fits.gz'], /noshell
   thismem = memory()
   maxmem = thismem[3]
   splog, 'Maximum memory usage = ', maxmem/1.d6, ' MB'
   splog, 'Total time = ', systime(1)-t1, ' sec'

   return
end

;------------------------------------------------------------------------------
pro fieldmerge, run2d=run2d, indir=indir, mergerun2d=mergerun2d, programs=programs, $
  legacy=legacy, plates=plates, skip_specprimary=skip_specprimary, lite=lite, _EXTRA=Extra

RESOLVE_ALL, /QUIET, /SKIP_EXISTING, /CONTINUE_ON_ERROR
CPU, TPOOL_NTHREADS = 1  
  
  undefine, logfile
  if TAG_EXIST(Extra,'field',/QUIET) then begin
     logfile='fieldmerge'
     logfile=logfile+'-'+field_to_string(Extra.field)
     if TAG_EXIST(Extra,'mjd',/QUIET) then logfile=logfile+'-'+strtrim(Extra.mjd,2)
     logfile = djs_filepath(logfile, root_dir='')

     logfile=logfile+'.log'
     cpbackup, logfile
     splog, filename=logfile
     splog, 'Log file ' + logfile + ' opened ' + systime()
  endif 
  
   if keyword_set(mergerun2d) then begin
       conflist, outdir=getenv('BOSS_SPECTRO_REDUX'), plist=plist
       fieldmerge1, plist=plist, legacy=legacy, skip_specprimary=skip_specprimary, $
         plates=plates, lite=lite, _EXTRA=Extra

   endif else begin

       conflist, plist=plist, topdir=indir, run2d=run2d
       ;print, plist;indir
       if (NOT keyword_set(plist)) then return

       if keyword_set(programs) then begin
          splog, 'Selecting only plates with programname:'
          splog, programs
          nmatch = lonarr(n_elements(plist))         
          for i=0, n_elements(programs)-1 do $
             nmatch+= strmatch(strtrim(plist.PROGRAMNAME,2), strtrim(programs[i],2) ) 
          indx = where( nmatch GT 0, ct)
          if (ct EQ 0) then begin
             splog, 'No plates found with programnames : '
             splog, programs
             return 
          endif
          plist=plist[indx]
          fieldmerge1, plist=plist, run2d=run2d, legacy=legacy, skip_specprimary=skip_specprimary, $
            plates=plates, lite=lite, _EXTRA=Extra
          return
       endif



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
       splog, alldir
       for i=0, n_elements(alldir)-1 do $
        fieldmerge1, run2d=alldir[i], indir=indir, legacy=legacy, skip_specprimary=skip_specprimary, $
         plates=plates, lite=lite, _EXTRA=Extra

   endelse
   
   splog, 'Successful completion of FIELDMERGE at ' + systime()
   ;----------
   ; Close log files
   splog, /close
end
;------------------------------------------------------------------------------
