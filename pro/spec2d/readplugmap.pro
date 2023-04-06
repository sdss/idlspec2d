;+
; NAME:
;   readplugmap
;
; PURPOSE:
;   Read plugmap file and append tags as requested
;
; CALLING SEQUENCE:
;   plugmap = readplugmap( plugfile, [ spectrographid, plugdir=, $
;    /apotags, /deredden, /calibobj, exptime=, hdr=, fibermask=, _EXTRA= ] )
;
; INPUTS:
;   plugfile  - Name of Yanny-parameter plugmap file
;
; OPTIONAL INPUTS:
;   spectrographid  - The spectrograph number, either 1 or 2;
;                     if not set (or 0), then return all object fibers
;   plugdir   - Directory for PLUGFILE
;   apotags   - If set, then add a number of tags to the output structure
;               constructed from the Yanny header.  These tags are:
;               CARTID, PLATEID, TILEID, RAPLATE, DECPLATE, REDDEN_MED.
;               Also add the tags FIBERSN[3], SYTHMAG[3] which are used
;               by the on-the-mountain reductions.
;   deredden  - If set, then deredden the MAG fields using the median
;               reddening value for the entire plate as found in the
;               Yanny header of the plugmap file; this is done for the
;               on-the-mountain reductions.
;   exptime   - Default exposure time SCI_EXPTIME to add to the output;
;               if there are multiple pointings and EXPTIME is set, then
;               the exposure time in each pointing is scaled such that
;               their sum is EXPTIME.
;   calibobj  - If set, then add a CALIBFLUX,CALIBFLUX_IVAR entries based upon
;               the calibObj (deprecated) or photoPlate files.
;               For stellar objects, this contains the
;               PSF fluxes in nMgy.  For galaxies, it contains the fiber fluxes
;               multiplied by the median (PSF/fiber) flux ratio for stars.
;               The MAG fields are left unchanged.
;               For objects with no calibObj entry, simply set these fields as:
;                 CALIBFLUX = 22.5 - 2.5*alog10(MAG), CALIBFLUX_IVAR = 0.
;               We apply the putative AB corrections to these fluxes
;               (but not to the MAG values).
;               Also add the SFD reddening values as SFD_EBV.
;   _EXTRA    - Keywords for PLUG2TSOBJ(), such as MJD,INDIR
;
; OUTPUTS:
;   plugmap   - Plugmap structure
;
; OPTIONAL OUTPUTS:
;   hdr       - Header from Yanny-formatted plugmap file
;   fibermask - Byte array with bits set for unknown fibers
;
; COMMENTS:
;   Do not use the calibObj structure if more than 10% of the non-sky
;   objects do not have fluxes.
;
;   Reads $IDLSPEC2D_DIR/opfiles/washers.par for ZOFFSET status overrides.
;   The original plugmap files reflect what we wanted to do; the overrides
;   and the return of this function reflect what we actually did.
;
; EXAMPLES:
;
; BUGS:
;   The AB corrections are hard-wired to be the same as in the photoop
;   product as of 18 Feb 2004.
;
; PROCEDURES CALLED:
;   djs_filepath()
;   euler
;   splog
;   yanny_par_fc()
;   yanny_read
;   FILE_TEST
;   prerun_readplugmap
;
; INTERNAL FUNCTIONS CALLED:
;   fits_to_yanny_hdr
;   get_parent_plugfile
;   clean_fibermap
;
; REVISION HISTORY:
;   29-Jan-2001  Written by S. Burles, FNAL
;   07-Aug-2012  Added ZOFFSET overrides; S. Bailey, LBL
;   21-Jun-2021  Editted by S Morrison to add correction check for catalog
;   15-Oct-2021  Modified by S Morrison to prepare for FPS:
;                       -adding readFPSobsSummary, ApogeeToBOSSRobomap, BossToApogeeRobomap
;                       -breaking readPlateplugMap in to a new function to preseve for FPS
;   15-Nov-2021  Modified by S Morrison to read and fill plugmap once and save it in fits file
;-
;------------------------------------------------------------------------------

function fits_to_yanny_hdr, hdr
    @plugmapkeys.idl

    nhead=n_elements(hdr)
    keys=strmid(hdr,0,8)
    yanny_hdr=strarr(nhead)
    for i=0, nhead-1 do begin
        key = strtrim(keys[i],2)
        if key eq 'EXTNAME' then continue
        matched_id=where(key_match_dict.Values() EQ key, ct)
        if ct eq 0 then begin
            matched_id2=where(key_match_dict2.Values() EQ key, ct2)
            if (ct2 eq 0) then continue
            matched_key = (key_match_dict2.Keys())[matched_id2]
            yanny_hdr[i]= matched_key+' '+string(sxpar(hdr, key))
        endif else begin
            matched_key = (key_match_dict.Keys())[matched_id]
            yanny_hdr[i]= matched_key+' '+string(sxpar(hdr, key))
        endelse
    endfor
    yanny_hdr=yanny_hdr[where(strlen(yanny_hdr) gt 0)]
    return, yanny_hdr
end

;------------------------------------------------------------------------------

function merge_fmap_datamodel, fibermap, plates=plates, legacy=legacy
    cols=Dictionary("assigned", 1, "on_target", 1, "valid", 1, "decollided", 0, $
                    "carton_to_target_pk", -999L, "cadence", "plates", $
                    "RACAT", !Values.D_NAN, "DECCAT", !Values.D_NAN, $
                    "COORD_EPOCH", 2000.0)

    foreach key, cols.Keys() do begin
        if not tag_exist(fibermap, key) then $
            fibermap = struct_addtags(fibermap,replicate(create_struct(key, cols[key]), $
                                        n_elements(fibermap)))
    endforeach
    
    if keyword_set(plates) or keyword_set(legacy) then begin
        masked = where(fibermap.fibermask eq 1, ctmaksed)
        if ctmaksed gt 0 then begin
            fibermap[masked].assigned = 0
            fibermap[masked].on_target = 0
            fibermap[masked].valid = 0
            fibermap[masked].decollided = 1
        endif
;        fibermap.RACAT = fibermap.RA
;        fibermap.DECCAT = fibermap.DEC

    endif
       
    return, fibermap
end

;------------------------------------------------------------------------------

function get_parent_plugfile, pf, cloned_from, plugdir=plugdir
   ;plugdir = file_dirname(file_dirname(pf),/MARK_DIRECTORY)
   if not keyword_set(plugdir) then begin
      thisplug = FILE_DIRNAME(FILE_DIRNAME(pf))
      confile = (findfile(filepath('confSummaryF-'+strtrim(cloned_from,2)+'.par',$
                                    root_dir=thisplug, subdir='*'), count=ct))[0]
      if (ct ne 0) then ppf = 'confSummaryF-'+strtrim(cloned_from,2)+'.par' $
                   else ppf = 'confSummary-'+strtrim(cloned_from,2)+'.par'
      pfp = ppf
   endif else begin
      confile = (findfile(filepath('confSummaryF-'+strtrim(cloned_from,2)+'.par',$
                                    root_dir=plugdir, subdir='*'), count=ct))[0]
      if (ct ne 0) then ppf = 'confSummaryF-'+strtrim(cloned_from,2)+'.par' $
                   else ppf = 'confSummary-'+strtrim(cloned_from,2)+'.par'
      ppf = (findfile(djs_filepath(ppf, root_dir=plugdir, subdir='*'), count=ct))[0]
   endelse
   return, ppf
end
;------------------------------------------------------------------------------

function clean_fibermap, fibermap, plates=plates
   
    float_cols=['RACAT','DECCAT','PMRA','PMDEC','PARALLAX','LAMBDA_EFF',$
                    'COORD_EPOCH','BP_MAG','GAIA_G_MAG','RP_MAG','H_MAG',$
                    'GAIA_PARALLEX', 'GAIA_PMRA', 'GAIA_PMDEC']
    float_lists=['MAG','CATDB_MAG','FIBER2MAG','PSFMAG','WISE_MAG','TWOMASS_MAG','GUVCAT_MAG']
    int_cols=['CARTON_TO_TARGET_PK']
    
    colnames = TAG_NAMES(fibermap)
        
    splog, 'Cleaning Fibermap'
    
    foreach col, float_cols do begin
        if tag_exist(fibermap, col) eq 0 then CONTINUE
        temp = fibermap.(where(colnames eq col))
        inan = where(temp eq -999., ct)
        if ct eq 0 then continue
        temp[inan] = !Values.F_NAN
        fibermap.(where(colnames eq col)) = temp
    endforeach
    
    foreach col, float_lists do begin
        if tag_exist(fibermap, col) eq 0 then CONTINUE
        temp = fibermap.(where(colnames eq col))
        nfilt = n_elements(temp[*,0])
        for ifilt = 0,nfilt-1 do begin
            tempf = temp[ifilt,*]
            inan = where(tempf eq -999., ct)
            if ct eq 0 then continue
            tempf[inan] = !Values.F_NAN
            temp[ifilt, *] = tempf
        endfor
        fibermap.(where(colnames eq col)) = temp
    endforeach

    return, fibermap
end

;------------------------------------------------------------------------------


function readplugmap, plugfile, spectrographid, plugdir=plugdir, savdir=savdir, $
    apotags=apotags, deredden=deredden, exptime=exptime, calibobj=calibobj, $
    hdr=hdr, fibermask=fibermask, plates=plates, legacy=legacy, gaiaext=gaiaext, $
    MWM_fluxer=MWM_fluxer, nfiles=nfiles, ccd=ccd, clobber=clobber, cartid=cartid,$
    no_db=no_db, _EXTRA=KeywordsForPhoto
    
    if keyword_set(plates) or keyword_set(legacy) then begin
        yanny_read, (findfile(djs_filepath(plugfile[0], root_dir=plugdir), count=ct))[0], junk, hdr=filehdr, /anonymous
        fieldid = field_to_string((yanny_par_fc(filehdr, 'plateId'))[0])
        if keyword_set(KeywordsForPhoto) then begin
            junk = where(strmatch(tag_names(KeywordsForPhoto), 'mjd',/fold_case),ct)
            if ct ne 0 then mjd = strtrim(KeywordsForPhoto.MJD,2)
        endif else mjd = (yanny_par_fc(filehdr, 'fscanMJD'))[0]
        if not keyword_set(cartid) then cartid=(yanny_par_fc(filehdr, 'cartridgeId'))[0]
        confid = string((yanny_par_fc(filehdr, 'fscanId'))[0],format='(i2.2)')
        if strmatch(plugfile, 'plPlugMapM-*-*-*[az].par', /FOLD_CASE) then  $
                    confid = confid+(yanny_par_fc(filehdr, 'pointing'))[0]
        if keyword_set(apotags) AND keyword_set(ccd) then begin
            mapfits_name = 'spfibermap-'+fieldid+'-'+mjd+'-'+ccd+'.fits'
        endif else begin
            mapfits_name = 'spfibermap-'+fieldid+'-'+mjd+'.fits'
        endelse
        if keyword_set(savdir) then mapfits_name=djs_filepath(mapfits_name, root_dir=savdir)
        
    endif else begin
        if keyword_set(plugdir) then $
            yanny_read, (findfile(djs_filepath(plugfile[0], root_dir=plugdir, subdir='*'), count=ct))[0], junk, hdr=filehdr, /anonymous $
        else yanny_read, plugfile[0], junk, hdr=filehdr, /anonymous
        fieldid = field_to_string(long(yanny_par_fc(filehdr, 'field_id')))
        mjd = (yanny_par_fc(filehdr, 'MJD'))[0]
        confid = (yanny_par_fc(filehdr, 'configuration_id'))[0]
        if keyword_set(apotags) AND keyword_set(ccd) then begin
            mapfits_name = 'spfibermap-'+fieldid+'-'+mjd+'-'+ccd+'.fits'
        endif else begin
            mapfits_name = 'spfibermap-'+fieldid+'-'+mjd+'.fits'
        endelse
        if keyword_set(savdir) then mapfits_name=djs_filepath(mapfits_name, root_dir=savdir)
     endelse

    if keyword_set(clobber) then file_delete, mapfits_name, /ALLOW_NONEXISTENT, /QUIET

    fits_fibermap = FILE_TEST(mapfits_name)
    if n_elements(plugfile) gt 1 then nfiles=1 else nfiles=0

    hdr = ['cut']
    plugmap = []
    fibermask = [-100]
    if keyword_set(plugdir) then plugdir_raw = plugdir
    foreach pf, plugfile do begin
        if keyword_set(plugdir_raw) then plugdir = plugdir_raw
        if keyword_set(plugdir) then begin
            pff = findfile(djs_filepath(pf, root_dir=plugdir, subdir ='*'), count = ct)
            if ct eq 0 then pff = findfile(djs_filepath(pf, root_dir=plugdir), count = ct)
            yanny_read, pff[0], junk, hdr=filehdr, /anonymous
        endif else yanny_read, pf, junk, hdr=filehdr, /anonymous
        if keyword_set(plates) or keyword_set(legacy) then begin
            ppf = pf
        endif else begin
            cloned_from_conf = yanny_par_fc(filehdr, 'cloned_from', count = cnt)
            cloned_from_conf = strtrim(string(cloned_from_conf,format='(i)'),2)
            cloned_to_conf = yanny_par_fc(filehdr, 'configuration_id')
            cloned_to_design = yanny_par_fc(filehdr, 'design_id')
            if not (keyword_set(apotags) and  (not FILE_TEST(mapfits_name))) then begin
                if (( cloned_from_conf NE '-999') AND (cnt NE 0)) then begin
                    ppf = get_parent_plugfile(pf,cloned_from_conf,plugdir=plugdir)
                    if keyword_set(plugdir) then ppf = file_basename(ppf)
                    splog, file_basename(pf), ' is cloned from ',file_basename(ppf)
                endif else ppf = pf
             endif else ppf = pf
        endelse

        fits_fibermap = FILE_TEST(mapfits_name)
        create=0
        exist_ext=0
        if (not fits_fibermap) then begin
            splog, 'No fits fiber map exists: '+mapfits_name
            create=1
            exist_ext=0
        endif else begin
            hdr_struct=MRDFITS(mapfits_name, 'Summary', sumhdr,/silent, status=st)
            if st ne 0 then begin
                ct = 0
                file_delete, mapfits_name, /ALLOW_NONEXISTENT
            endif else begin
                map_ext = where(strmatch(hdr_struct.EXTNAME, file_basename(pf)+"*", /fold_case),ct)
            endelse
            map_ext = where(strmatch(hdr_struct.EXTNAME, file_basename(pf)+"*", /fold_case),ct)
            if ct NE 0 then begin
                splog, 'Reading fits fiber map extension for '+file_basename(pf)+' from '+mapfits_name
                fibermap = mrdfits(mapfits_name, file_basename(pf), fhdr1,/silent)
                plugmap = [plugmap, fibermap]
                hdr1 = struct_to_yannyhdr(file_basename(pf),hdr_struct=hdr_struct)
                hdr = [hdr, hdr1,'cut']
                fibermask = [fibermask, fibermap.fibermask, -100]
                create=0
                exist_ext=file_basename(pf)
            endif else begin
                create = 1
                exist_ext=0
            endelse
            if ((pf ne ppf) AND (create eq 1)) then begin
                map_ext = where(strmatch(hdr_struct.EXTNAME, file_basename(ppf)+"*", /fold_case),ct)
                if ct NE 0 then begin
                    splog, 'Reading fits fiber map extension '+file_basename(ppf)+' cloned for '+file_basename(pf)+' from '+mapfits_name

                    fibermap = mrdfits(mapfits_name, file_basename(ppf), hdr1,/silent)
                    undefine, fits_hdr
                    sxaddpar, fits_hdr, 'EXTNAME', FILE_BASENAME(pf), ' Complete Plugmap/confSummary'
                    MWRFITS, fibermap, mapfits_name, fits_hdr, Status=Status, /silent

                    if keyword_set(plugdir) then begin
                        tpf = findfile(djs_filepath(pf, root_dir=plugdir, subdir='*'), count = ct)
                        if ct eq 0 then tpf = findfile(djs_filepath(pf, root_dir=plugdir), count = ct)
                        pf = tpf
                    endif
                    
                    yanny_read, pf, junk, hdr=hdr1, /anonymous
                    hdr_struct=yannyhdr_to_struct(hdr1, pf, plugdir=plugdir, hdr_struct=hdr_struct,cartid=cartid)

                    update_fitstable, mapfits_name, 'Summary', ' Table of plugmap/confSummary header parameters', hdr_struct
                    ;modfits, mapfits_name, hdr_struct, EXTNAME='Summary'

                    plugmap = [plugmap, fibermap]
                    hdr = [hdr, hdr1,'cut']
                    fibermask = [fibermask, fibermap.fibermask, -100]
                    create = 0
                    exist_ext=file_basename(pf)
                endif
            endif
        endelse
        if keyword_set(create) then begin
            splog, 'Creating New fits fiber map extension from '+ppf

            fibermap = prerun_readplugmap(ppf, mapfits_name, plugdir=plugdir, apotags=apotags, $
                                          exptime=exptime, hdr=hdr1, fibermask=fibermask1, $
                                          cartid=cartid, nfiles=nfiles, no_db=no_db, $
                                          plates=plates, legacy=legacy, _EXTRA=KeywordsForPhoto)

            plugmap = [plugmap, fibermap]
            ;hdr = [hdr, fits_to_yanny_hdr(hdr1),'cut']
            hdr = [hdr, hdr1,'cut']
            fibermask = [fibermask, fibermap.fibermask, -100]
            Undefine, fibermask1
        endif
        Undefine, hdr1
    endforeach
    hdr = [hdr,'cut', '  ']
    fibermask=[fibermask,-100]
            
    fieldid = (yanny_par_fc(hdr, 'field_id'))[0]
    ra_field=float(yanny_par_fc(hdr, 'raCen'))
    dec_field=float(yanny_par_fc(hdr, 'decCen'))

    if keyword_set(plates) then programname = yanny_par_fc(hdr, 'programname')

    ;if keyword_set(plates) and keyword_set(programname) then $
    ;          stdtype = 'SPECTROPHOTO_STD' else stdtype = 'standard_boss'
    stdtype = 'SPECTROPHOTO_STD'
    addtags = replicate(create_struct( $
                    'EBV',!Values.F_NAN, $
                    'EBV_TYPE', 'SFD'), n_elements(plugmap))
    plugmap = struct_addtags(plugmap, addtags)
    if not keyword_set(apotags) then plugmap.EBV=plugmap.sfd_ebv
    if keyword_set(calibobj) then begin
        rjce_extintion = 0
        if keyword_set(rjce_extintion) then begin
            if keyword_set(plates) and keyword_set(programname) then begin
                if ((strmatch(programname, '*MWM*', /fold_case) eq 1) $
                    || (strmatch(programname, '*OFFSET*', /fold_case) eq 1)) then begin
                    splog, "Using RJCE extintion"
                    spht = strmatch(plugmap.objtype, stdtype, /FOLD_CASE)
                    ispht = where(spht, nspht)
                    ebv = plugmap[ispht].sfd_ebv
                    EBV_TYPE = plugmap[ispht].EBV_TYPE
                    for i=0, n_elements(plugmap[ispht])-1 do begin
                        dat=(plugmap[ispht])[i].ebv_rjce
                        if (finite(dat) ne 0) and (dat gt 0) and (dat le 1.2*(plugmap[ispht])[i].sfd_ebv) then begin
                                ebv[i] = (plugmap[ispht])[i].ebv_rjce
                                EBV_TYPE[i] = 'RJCE'
                        endif
                    endfor
                    plugmap[ispht].ebv = ebv
                    plugmap[ispht].EBV_TYPE = EBV_TYPE
                    gaiaext = 0
                endif
            endif
        endif

        if keyword_set(MWM_fluxer) then begin
            if keyword_set(plates) then begin
                if keyword_set(programname) then begin
                    if ((strmatch(programname, '*MWM*', /fold_case) eq 1) $
                        || (strmatch(programname, '*OFFSET*', /fold_case) eq 1)) then begin
                        gaiaext = 1
                    endif
                endif
            endif else begin
                euler, ra_field[0], dec_field[0], ll_field, bb_field, 1
                if abs(bb_field) lt 15 then gaiaext = 1
            endelse
        endif
        if keyword_set(no_db) then gaiaext = 0
        if keyword_set(gaiaext) then begin
            splog, "Using dust_3d_map"
            ebv = plugmap.sfd_ebv
            EBV_TYPE = plugmap.EBV_TYPE
            badstdmask = plugmap.badstdmask
            for i=0, n_elements(plugmap)-1 do begin
                    dat=plugmap[i].ebv_gaia
                    if (finite(dat) ne 0) and (dat ne -999) and (dat le plugmap[i].sfd_ebv) then begin
                            ebv[i] = plugmap[i].ebv_gaia
                            EBV_TYPE[i] = 'dust_3d_map'
                    endif ;else begin
                          ;  if strmatch(plugmap[i].objtype, stdtype, /FOLD_CASE) then badstdmask[i] = 1
                    ;endelse
            endfor
            plugmap.badstdmask = badstdmask
            plugmap.ebv = ebv
            plugmap.EBV_TYPE = EBV_TYPE
        endif
    endif

    cartid = (yanny_par_fc(hdr, 'observatory'))[0]
    if strmatch(cartid, '*LCO*',/fold_case) then begin
        sid = plugmap.spectrographid
        sid[where(plugmap.spectrographid eq 2)] = 0
        sid[where(plugmap.spectrographid eq 1)] = 2
        plugmap.spectrographid = sid
    endif
    if keyword_set(plates) or keyword_set(legacy) then begin
        if keyword_set(deredden) then begin
            ; They are just the ratios of new/old,
            redden_corr=[0.822308,0.870815,0.830607,0.813998,0.853955]
            if (keyword_set(deredden)) then begin
                splog, 'Applying reddening vector ', redden_med
                for ifilt=0, 4 do begin
                    if (plateid ge 7572) then begin
                        plugmap.mag[ifilt] = plugmap.mag[ifilt] - redden_med[ifilt]
                    endif else begin
                        plugmap.mag[ifilt] = plugmap.mag[ifilt] - redden_med[ifilt]*redden_corr[ifilt]
                        splog, 'modified extinction co-efficients: ',redden_med[ifilt]*redden_corr[ifilt]
                    endelse
                endfor
            endif
        endif
        if keyword_set(plates) then begin
            nfiber = 1000
            plugmap.spectrographid=1
        endif else begin
            nfiber = n_elements(plugmap)
        endelse
        if (keyword_set(spectrographid)) then begin
            indx = (spectrographid-1)*nfiber/2 + lindgen(nfiber/2)
            plugmap = plugmap[indx]
            fibermask = fibermask[indx]
        endif
    endif else begin
        if (keyword_set(spectrographid)) then begin
            outfibers=[]
            outplugmap=[]
            nt=where(fibermask EQ -100)
            for iflat=0, n_elements(plugfile)-1 do begin
                tmp_fibmask = fibermask[nt[iflat]+1:nt[iflat+1]-1]
                tmp_plugmap = plugmap[(nt[iflat]-iflat):(nt[iflat+1]-(2+iflat))]
        
                indx = where(tmp_plugmap.spectrographid eq spectrographid)
                tmp_plugmap = tmp_plugmap[indx]
                tmp_fibmask = tmp_fibmask[indx]
                
                outfibers = [outfibers, -100, tmp_fibmask]
                outplugmap = [outplugmap, tmp_plugmap]
            endfor
            fibermask = [outfibers,  -100,  -100]
            plugmap = outplugmap
        endif
    endelse
    if keyword_set(apotags) then fibermask=fibermask[where(fibermask ne -100)]
 
    tags_to_delete=['SFD_EBV','ebv_gaia','ebv_rjce']
    foreach tag, tags_to_delete do begin
             if tag_exist(plugmap,tag) then $
                   plugmap = struct_trimtags(plugmap,except_tags=[tag])
    endforeach
 
    plugmap=merge_fmap_datamodel(plugmap, plates=plates, legacy=legacy)
 
    if (not keyword_set(apotags)) then plugmap = clean_fibermap(plugmap)
 
    fibermask=fibermask
;    splog, plugmap
;    splot, fibermask    
;    struct_print, plugmap, filename='fibermap.html', /html
    return, plugmap
end
