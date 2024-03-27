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
;   05-May-2023  Modified by S Morrison to use spFibermap (produced by readfibermaps.py)
;-
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

pro run_readfibermap, spFibermap, spplan=spplan, $
        mjd=mjd, ccd=ccd, plugfile=plugfile, filehdr=filehdr
        
    if keyword_set(apotags) then begin
        flags = ' --clobber --SOS --log'
        flags = flags + ' --topdir /data/boss/sos/'+strtrim(mjd,2)
        flags = flags + ' --confSummary '+plugfile
        flags = flags + ' --ccd '+ccd
        flags = flags + ' --mjd '+mjd
        obs = (yanny_par_fc(filehdr, 'observatory'))[0]
        if strmatch(obs, '*LCO*',/fold_case) then $
            flags = flags + ' --lco'

        cmd = "readfibermaps.py "+ flags
        splog,cmd
        spawn, cmd, dat

    endif else  begin
        flags  = ' --spplan2d '+spplan
        flags  = flags + ' --clobber'
        cmd = "readfibermaps.py "+ flags
        splog,cmd
        spawn, cmd, dat
        splog, dat

    endelse
end

;------------------------------------------------------------------------------


function readplugmap, plugfile, spectrographid, plugdir=plugdir, savdir=savdir, $
    apotags=apotags, deredden=deredden, exptime=exptime, calibobj=calibobj, $
    hdr=hdr, fibermask=fibermask, plates=plates, legacy=legacy, $
    gaiaext=gaiaext, map3d = map3d, MWM_fluxer=MWM_fluxer,$
    nfiles=nfiles, ccd=ccd, cartid=cartid, no_db=no_db, $
    _EXTRA=KeywordsForPhoto
    
    
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
            spFibermap = 'spfibermap-'+fieldid+'-'+mjd+'-'+ccd+'.fits'
        endif else begin
            spFibermap = 'spfibermap-'+fieldid+'-'+mjd+'.fits'
	    spplan     = 'spPlan2d-'+fieldid+'-'+mjd+'.par'
        endelse
        if keyword_set(savdir) then spFibermap=djs_filepath(spFibermap, root_dir=savdir)
        
    endif else begin
        if keyword_set(plugdir) then begin
;splog, 'DEBUG:', djs_filepath(plugfile[0], root_dir=plugdir, subdir='*/*')
            yanny_read, (findfile(djs_filepath(plugfile[0], root_dir=plugdir, subdir='*/*'), count=ct))[0], junk, hdr=filehdr, /anonymous
        endif else begin
;splog, 'DEBUG:', plugfile[0]
            yanny_read, plugfile[0], junk, hdr=filehdr, /anonymous
        endelse
;splog, filehdr
        fieldid = field_to_string(long(yanny_par_fc(filehdr, 'field_id')))
        mjd = (yanny_par_fc(filehdr, 'MJD'))[0]
        if keyword_set(apotags) then begin
            if keyword_set(KeywordsForPhoto) then begin
                junk = where(strmatch(tag_names(KeywordsForPhoto), 'mjd',/fold_case),ct)
                if ct ne 0 then mjd = strtrim(KeywordsForPhoto.MJD,2)
            endif
        endif
        confid = (yanny_par_fc(filehdr, 'configuration_id'))[0]
        if keyword_set(apotags) AND keyword_set(ccd) then begin
            spFibermap = 'spfibermap-'+fieldid+'-'+mjd+'-'+ccd+'.fits'
        endif else begin
            spFibermap = 'spfibermap-'+fieldid+'-'+mjd+'.fits'
            spplan     = 'spPlan2d-'+fieldid+'-'+mjd+'.par'
        endelse
        if keyword_set(savdir) then spFibermap=djs_filepath(spFibermap, root_dir=savdir)
     endelse
    
    if not keyword_set(FILE_TEST(spFibermap)) then begin
        splog, 'Missing spFibermap file: ', spFibermap
        run_readfibermap, spFibermap, spplan=spplan, $
                mjd=mjd, ccd=ccd, plugfile=plugfile, filehdr=filehdr
        if not keyword_set(FILE_TEST(spFibermap)) then $
            message, 'Missing spFibermap file: '+ spFibermap
            
    endif
  
    hdr_struct=MRDFITS(spFibermap, 'SUMMARY', sumhdr,/silent, status=st)

    if st ne 0 then begin
        message, 'Bad spFibermap file: '+ spFibermap
    endif
    
    hdr = ['cut']
    plugmap = []
    fibermask = [-100]
    if keyword_set(plugdir) then plugdir_raw = plugdir
    foreach pf, plugfile do begin
        map_ext = where(strmatch(hdr_struct.EXTNAME, file_basename(pf)+"*", /fold_case),ct)
        if ct eq 0 then begin
            splog,file_basename(pf)+ ' Missing from '+spFibermap
            run_readfibermap, spFibermap, spplan=spplan, $
                mjd=mjd, ccd=ccd, plugfile=plugfile, filehdr=filehdr
            map_ext = where(strmatch(hdr_struct.EXTNAME, file_basename(pf)+"*", /fold_case),ct)
            if ct eq 0 then begin
                message, file_basename(pf)+ ' Missing from '+spFibermap
            endif
        endif
        splog, 'Reading fits fiber map extension for '+file_basename(pf)+' from '+spFibermap
        fibermap = mrdfits(spFibermap, STRUPCASE(file_basename(pf)), fhdr1,/silent)

        if keyword_set(plugmap) then begin
            nflags = max([n_elements(plugmap[0].SDSS5_TARGET_FLAGS),n_elements(fibermap[0].SDSS5_TARGET_FLAGS)])
            plugmap = rename_tags(plugmap, 'SDSS5_TARGET_FLAGS','SDSS5_TARGET_FLAGS_raw')
            fibermap = rename_tags(fibermap, 'SDSS5_TARGET_FLAGS','SDSS5_TARGET_FLAGS_raw')
            plugmap = struct_addtags(plugmap, $
                                     replicate(create_struct('SDSS5_TARGET_FLAGS', BYTARR(nflags)),$
                                              n_elements(plugmap)))
            plugmap.SDSS5_TARGET_FLAGS[0:n_elements(plugmap[0].SDSS5_TARGET_FLAGS_raw)-1] = plugmap.SDSS5_TARGET_FLAGS_raw
            plugmap = struct_trimtags(plugmap, except_tags='SDSS5_TARGET_FLAGS_RAW')
            fibermap = struct_addtags(fibermap, $
                                      replicate(create_struct('SDSS5_TARGET_FLAGS', BYTARR(nflags)),$
                                               n_elements(fibermap)))
           fibermap.SDSS5_TARGET_FLAGS[0:n_elements(fibermap[0].SDSS5_TARGET_FLAGS_raw)-1] = fibermap.SDSS5_TARGET_FLAGS_raw
           fibermap = struct_trimtags(fibermap,except_tags='SDSS5_TARGET_FLAGS_RAW')
        endif
        plugmap = [plugmap, fibermap]
        hdr1 = struct_to_yannyhdr(file_basename(pf),hdr_struct=hdr_struct)
        hdr = [hdr, hdr1,'cut']
        fibermask = [fibermask, fibermap.fibermask, -100]

    endforeach
    hdr = [hdr,'cut', '  ']
    fibermask=[fibermask,-100]
        
    fieldid = (yanny_par_fc(hdr, 'field_id'))[0]
    ra_field=float(yanny_par_fc(hdr, 'raCen'))
    dec_field=float(yanny_par_fc(hdr, 'decCen'))

    if keyword_set(plates) then programname = yanny_par_fc(hdr, 'programname')

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
            if not keyword_set(map3d) then map3d = 'merge3d'
            map3d = STRLOWCASE(map3d)
            splog, "Using ",map3d," dust_3d_map"
            ebv = plugmap.sfd_ebv
            EBV_TYPE = plugmap.EBV_TYPE
            badstdmask = plugmap.badstdmask
            for i=0, n_elements(plugmap)-1 do begin
                    case map3d of
                        'bayestar15': begin
                                     dat=plugmap[i].EBV_BAYESTAR15
                                     tmap3d = 'bayestar15'
                                end
                        'bay15': begin
                                     dat=plugmap[i].EBV_BAYESTAR15
                                     tmap3d = 'bayestar15'
                                end
;                        'edenhofer2023': begin
;                                     dat=plugmap[i].EBV_EDENHOFER2023
;                                     tmap3d = 'edenhofer2023'
;                                end
;                        'eden23': begin
;                                     dat=plugmap[i].EBV_EDENHOFER2023
;                                     tmap3d = 'edenhofer2023'
;                                end
                         'merge3d': begin
                                     dat=plugmap[i].EBV_3D
                                     tmap3d = plugmap[i].EBV_3DSRC
                                 end
                         else: begin
                                     dat=plugmap[i].EBV_BAYESTAR15
                                     tmap3d = 'bayestar15'
                               end
                    endcase
                    if (finite(dat) ne 0) and (dat ne -999) and (dat le plugmap[i].sfd_ebv) then begin
                            ebv[i] = dat
                            EBV_TYPE[i] = tmap3d
                    endif
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
            nt=where(fibermask ne -100,ct)
            if ct ne 0 then fibermask = fibermask[nt]
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
 
    tags_to_delete=['SFD_EBV','EBV_BAYESTAR15','ebv_rjce']
    foreach tag, tags_to_delete do begin
             if tag_exist(plugmap,tag) then $
                   plugmap = struct_trimtags(plugmap,except_tags=[tag])
    endforeach
  
    if (not keyword_set(apotags)) then plugmap = clean_fibermap(plugmap)
 
    fibermask=fibermask
;    splog, plugmap
;    splog, fibermask
; help, plugmap
;    MWRFITS, plugmap, 'fibermap.fits', Status=Status, /silent
;    struct_print, plugmap, filename='fibermap.html', /html
    return, plugmap
end
