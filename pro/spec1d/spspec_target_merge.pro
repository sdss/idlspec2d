

;+
; NAME:
;   spspec_target_merge
;
; PURPOSE:
;   To create spSpec and spFullsky target level coadds (independent of field-mjd) 
;
; CALLING SEQUENCE:
;
; INPUTS:
;   customplan - The spPlanCustom file for the coadd
;
; OPTIONAL KEYWORDS:
;   topdir - the daily coadd base directory
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;
;-
;------------------------------------------------------------------------------
; redefine splog_multi as splog
; (so all functions/procedure called by this procedure that use splog use the parallel logger)

pro splog, noname=noname, prelog=prelog, $
 filename=filename, append=append, close_all=close_all,$
 secondary=secondary, close_secondary=close_secondary, $
 no_stdout=no_stdout, fname=fname, close=close,$
 v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, _EXTRA=extra

    help, calls=calls
    fname = (str_sep(calls[1], ' '))[0] + ': '
    for i=0,n_elements(calls)-4 do fname = ' ' + fname

    if keyword_set(close) then close_all = 1
    splog_multi, noname=noname, prelog=prelog, $
        filename=filename, append=append, close_all=close_all,$
        secondary=secondary, close_secondary=close_secondary, $
        no_stdout=no_stdout, fname=fname, nv=n_params(),$
        v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, _EXTRA=extra
END

;------------------------------------------------------------------------------

pro spspec_target_merge, customplan, topdir=topdir, mjd = mjd

    RESOLVE_ROUTINE,'sdss_maskbits',/EITHER,/SKIP_EXISTING, /quiet
    RESOLVE_ALL, /SKIP_EXISTING, /quiet, /CONTINUE_ON_ERROR
    CPU, TPOOL_NTHREADS = 1


    binsz = 1.0d-4
    wavemin=3.5523d
    wavemax=4.0171d
    camnames = ['B1','R1','B2','R2']
    nord = 3
  
    exclude_mask = (fibermask_bits('NOPLUG') OR $
                    fibermask_bits('BADTRACE') OR $
                    fibermask_bits('BADFLAT') OR $
                    fibermask_bits('BADARC') OR $
                    fibermask_bits('NODATA'))

    allseq = yanny_readone(customplan, 'COADDPLAN', hdr=hdr, /anon)
    if not keyword_set(topdir) then begin
        topdir=filepath(getenv('RUN2D'), root_dir=getenv('BOSS_SPECTRO_REDUX'))
    endif
    logfile = repstr(repstr(customplan,'.par','.log'), 'spPlanCustom','spDiagcomb')

    customplan = fileandpath(customplan, path=custom_dir)
    if not keyword_set(custom_dir) then $
        custom_dir = get_field_dir(topdir,'',yanny_par(hdr,'NAME'),/custom)
    custom  = yanny_par(hdr,'NAME')
    runmjd  = yanny_par(hdr,'CreateMJD')
    targid  = yanny_par(hdr,'TARGID')
    mjd_flag = 0
    if keyword_set(mjd) then begin
        IF n_elements(mjd) gt 1 then begin
            mjd_flag = strtrim(min(mjd),2)+'_'+strtrim(max(mjd),2)
        endif else mjd_flag = strtrim(mjd,2)
        logfile = repstr(logfile, '.log','_'+mjd_flag+'.log')
        if n_elements(mjd) gt 1 then mjd_flag = 0
    endif
    cpbackup, logfile
    splog, filename=logfile
    splog, 'Log file ' + logfile + ' opened ' + systime()

    epoch_combine = allseq.EPOCH_COMBINE
    if keyword_set(mjd) then epoch_combine = [mjd]
    epoch_combine = epoch_combine[UNIQ(epoch_combine, SORT(epoch_combine))]
    foreach ec, epoch_combine, ic do begin
        epseq = allseq[where(allseq.EPOCH_COMBINE eq ec)]
        if not keyword_set(mjd_flag) then begin
            logfile_epoch = repstr(logfile, '.log','_'+strtrim(ec,2)+'.log')
            splog, secondary = logfile_epoch
            splog, 'Starting Epoch Combine Logging to '+logfile_epoch+' at '+systime()
        endif
        splog, 'Building spSpec ',custom,' files for '+strtrim(ec,2)+' ('+strtrim(ic+1,2)+'/'+strtrim(n_elements(epoch_combine),2)+')'
        splog, ''
        nspec = strtrim(n_elements(epseq),2)
        foreach targ, epseq, target_idx do begin
            cid = max(targ.CATALOGID_LIST)
            splog, 'spSpec-'+custom+'-'+strtrim(ec,2)+'-'+strtrim(cid,2)+'.fits ('+strtrim(target_idx+1,2)+'/'+nspec+')'

            mjds = targ.MJD_LIST
            fields = targ.FIELDS_LIST
            fmjds = targ.FMJD_LIST
            nfiles = n_elements(fields)
            nexp = 0
            exptime = 0
            npix = []
            hdrs = []
            exps = []
            obs  = []
            nexpvec = [0,0,0,0]
            spspecfiles = []
            nexps = []
            exptimevec = [0.0,0.0,0.0,0.0]
            valid_tar = 1
            for i = 0, n_elements(targ.FIELDS_LIST)-1 do begin
                if fields[i] eq -1 then continue
                foreach cid, targ.CATALOGID_LIST do begin
                    spspecfile = filepath('spSpec-'+fmjds[i]+'-'+strtrim(cid,2)+'.fits', $
                                            root_dir = get_field_dir(topdir, '', fields[i]),$
                                            subdirectory=['coadd',strtrim(mjds[i],2)])
                    valid_tar = File_test(spspecfile,/READ)
                    if keyword_set(valid_tar) then break
                endforeach
                if not keyword_set(valid_tar) then begin
                    foreach cid, targ.CATALOGID_LIST do begin
                    spspecfile = filepath('spSpec-'+fmjds[i]+'-'+strtrim(cid,2)+'.fits', $
                                            root_dir = get_field_dir(topdir, '', fields[i]),$
                                            subdirectory=['coadd',strtrim(mjds[i],2)])
                        splog, 'Missing specfile ',fileandpath(spspecfile),' for '+targid+':',strtrim(targ.TARGID,2),' SKIPPING'
                    endforeach
                    break
                endif
                print, spspecfile, format = '(%"\t%s\t")'
                spspecfiles = [ spspecfiles, spspecfile]
                
                fibermap = mrdfits(spspecfile, 2,/silent)
                bighdr = headfits(spspecfile, EXTEN=0)
                
                foreach cam, camnames, icam do begin
                    exptimevec[icam] = exptimevec[icam] + sxpar(bighdr, 'EXPT_'+strtrim(cam,2))
                    nexpvec[icam]    = nexpvec[icam]    + sxpar(bighdr, 'NEXP_'+strtrim(cam,2))
                endforeach
                exptime = exptime + sxpar(bighdr, 'EXPTIME')
                nexp = nexp + fibermap.NEXP
                nexps = [nexps, fibermap.NEXP]
                obs  = [obs, sxpar(bighdr,'OBSERVAT')]
                temp = mrdfits(spspecfile, 3,/silent)
                tempflux   = temp.flux
                npix = [npix,  n_elements(tempflux)]
                if i eq 0 then begin
                    temploglam = temp.loglam
                    tempivar   = temp.ivar
                    tempamask  = temp.and_mask
                    tempormask = temp.or_mask
                    tempwdisp  = temp.wdisp
                    tempsky    = temp.sky
                    tempwresl  = temp.wresl
                    hdrplug    = headfits(spspecfile, EXTEN=2)
                endif
            endfor
            if not keyword_set(valid_tar) then continue
       
            ; Clears cards
            rawhdr = bighdr
            for i = 0, n_elements(rawhdr)-1 do begin
                fcard = rawhdr[i]
                fcard = strtrim((strsplit(fcard, '=',/extract))[0],2)
                if strmatch(fcard, 'EXPID*', /fold_case) then begin
                    sxdelpar, bighdr, fcard
                endif
            endfor
       
            npixmax = max(npix)
            flux    = make_array(npixmax,nexp,type=size(tempflux,/type))
            ivar    = make_array(npixmax,nexp,type=size(tempivar,/type))
            llam    = make_array(npixmax,nexp,/DOUBLE)
            amask   = make_array(npixmax,nexp,type=size(tempamask,/type))
            ormask  = make_array(npixmax,nexp,type=size(tempormask,/type))
            wdisp   = make_array(npixmax,nexp,type=size(tempwisp,/type))
            sky     = make_array(npixmax,nexp,type=size(tempsky,/type))
            wresl   = make_array(npixmax,nexp,type=size(tempwresl,/type))

            bighdr = delhdrsec(bighdr, 'EXPOSURE SETTINGS')
            bighdr = clearhdrsec(bighdr, 'TELESCOPE INFO', exclude = ['RA','DEC','RADEG','DECDEG','EQUINOX'])
            bighdr = clearhdrsec(bighdr,'FIELD/PLATE INFO')
            bighdr = clearhdrsec(bighdr,'INSTRUMENT INFO')
            bighdr = delhdrsec(bighdr, 'SPECTROGRAPH STATUS')
            bighdr = delhdrsec(bighdr, 'REDUCTION')
            bighdr = delhdrsec(bighdr, 'VERSIONS')
            bighdr = delhdrsec(bighdr, 'APO WEATHER')
            bighdr = delhdrsec(bighdr, 'LCO WEATHER')
            
            cleankey = ['PRESSURE','WINDD','WINDS','GUSTD','GUSTS','AIRTEMP',$
                        'DEWPOINT','TRUSTEMP','HUMIDITY', 'DUSTA', 'DUSTB',$
                        'WINDD25M','WINDS25M','T_OUT','T_IN','T_PRIM','T_CELL',$
                        'T_FLOOR','T_TRUSS', 'CCDTEMP','LN2TEMP']
            foreach key, cleankey do begin
                sxdelpar, bighdr, key
            endforeach

            j = 0
            l = 0

            ffm = 0
            fmaps = []
            includes = []
            exclude = []
            foreach exp, spspecfiles, i do begin
                temp_fibermap = mrdfits(spspecfiles[i], 2,/silent)
                
                mask = temp_fibermap[0].fibermask AND exclude_mask eq 0
                if mask eq 0 then begin
                    fmaps =    [fmaps, ptr_new(temp_fibermap)]
                    includes = [includes, spspecfiles[i]]
                endif else begin
                    splog, 'dropping'
                
                endelse
            endforeach
            if n_elements(includes) gt 0 then begin
                spspecfiles = includes
                foreach e, exclude do $
                    splog, 'Excluding '+strtrim(e,2)+' due to fibermask'
            endif
            
            foreach exp, spspecfiles, i do begin
                temp_fibermap = *(fmaps[i])
;                temp_fibermap = mrdfits(spspecfiles[i], 2,/silent)
                nexp1  = temp_fibermap.NEXP
                
                for k = 0, nexp1-1 do begin
                    temp = mrdfits(spspecfiles[i], k+3,/silent)
                    hdr = headfits(spspecfiles[i], EXTEN=0)
                    npix = n_elements(temp.flux)
                    flux[0:npix-1,l]   = temp.flux
                    ivar[0:npix-1,l]   = temp.ivar
                    llam[0:npix-1,l]   = dindgen(npix) * binsz+wavemin;temp.loglam
                    amask[0:npix-1,l]  = temp.and_mask
                    ormask[0:npix-1,l] = temp.or_mask
                    wdisp[0:npix-1,l]  = temp.wdisp
                    sky[0:npix-1,l]    = temp.sky
                    wresl[0:npix-1,l]  = temp.wresl
                    l = l + 1
                    hdrs = [hdrs,ptr_new(headfits(spspecfiles[i], EXTEN=k+3))]
                    exps = [exps,ptr_new(temp)]
                endfor

           
                list_cols = ['MOON_DIST', 'MOON_PHASE', 'FIBERID_LIST', 'RA_LIST', 'DEC_LIST', $
                            'DELTA_RA_LIST', 'DELTA_DEC_LIST','EXPTIME','FIRSTCARTON_LIST', $
                            'CARTON_TO_TARGET_PK_LIST', 'ASSIGNED_LIST', 'ON_TARGET_LIST', $
                            'VALID_LIST', 'DECOLLIDED_LIST',  'TOO_LIST','XFOCAL_LIST', 'YFOCAL_LIST', $
                            'TAI_LIST', 'MJDLIST', 'DESIGNS', 'CONFIGS', 'AIRMASS_LIST', $
                            'FIELDSNR2G_LIST', 'FIELDSNR2R_LIST', 'FIELDSNR2I_LIST', $
                            'SEEING20_LIST', 'SEEING50_LIST', 'SEEING80_LIST']
                cart = sxpar(hdr, 'CARTID')
                if isa(cart, /NUMBER) then begin
                    temp_fibermap.ASSIGNED_LIST = strjoin(replicate('1',temp_fibermap.NEXP),' ')
                    temp_fibermap.ON_TARGET_LIST = strjoin(replicate('1',temp_fibermap.NEXP),' ')
                    temp_fibermap.VALID_LIST = strjoin(replicate('1',temp_fibermap.NEXP),' ')
                    temp_fibermap.DECOLLIDED_LIST = strjoin(replicate('0',temp_fibermap.NEXP),' ')
                    temp_fibermap.TOO_LIST = strjoin(replicate('0',temp_fibermap.NEXP),' ')
                    temp_fibermap.DELTA_RA_LIST = strjoin(replicate('0.0',temp_fibermap.NEXP),' ')
                    temp_fibermap.DELTA_DEC_LIST = strjoin(replicate('0.0',temp_fibermap.NEXP),' ')
                endif

                foreach col, list_cols do begin
                    inTaginx  = where(strmatch(tag_names(fibermap), col,/FOLD_CASE) eq 1, ct)
                    outTaginx = where(strmatch(tag_names(temp_fibermap), col,/FOLD_CASE) eq 1, ct1)
                    if (ct eq 0) or (ct1 eq 0) then continue
                    fibermap[0].(outTaginx) = strjoin([fibermap[0].(outTaginx),temp_fibermap[0].(inTaginx)],' ')
                endforeach
                
                ; join fibermask
                ffm = ffm OR temp_fibermap[0].fibermask
                fibermap[0].fibermask = ffm
                
                ;add cards
                foreach fcard, headfits(spspecfiles[i],EXTEN=0) do begin
                    if strmatch(fcard, 'EXPID*', /fold_case) then begin
                        card = (strsplit(fcard,'=',/extract))
                        if n_elements(card) eq 1 then continue
                        cardv = (strsplit(card[1],/extract))[0]
                        sxaddpar, bighdr, string('EXPID',j+1, format='(a5,i3.3)'), cardv, $
                                        ' ID string for exposure '+strtrim(j+1,2), before='EXPTIME'
                        j=j+1
                    endif
                endforeach
            endforeach

            nfinalpix = long((wavemax - wavemin)/binsz)
            finalwave = dindgen(nfinalpix) * binsz + wavemin

            bestandmask= amask[*,*]
            bestormask = ormask[*,*]
            temppixmask = amask[*,*]
            bestresolution = wresl[*,0]
            rm_combine1fiber, llam, flux, ivar, indisp=wdisp, skyflux=sky, $
                            inormask=ormask, inandmask=amask, inresl=wresl, $
                            binsz=binsz, nord=nord, bkptbin=bkptbin, maxsep=maxsep, $
                            maxiter=0, upper=3d6, lower=3d6, maxrej=1, $            ; for _EXTRA DJS_REJECT
                            newloglam=finalwave, newflux=combinedflux, newivar=combinedivar, $
                            finalmask=temppixmask, andmask=bestandmask, ormask=bestormask, $
                            newdisp=bestwdisp,newsky=bestsky, newresl=bestresolution

            Assigned  = fix(strsplit(fibermap.ASSIGNED_LIST,/extract))
            valid     = fix(strsplit(fibermap.VALID_LIST,/extract))
            On_target = fix(strsplit(fibermap.ON_TARGET_LIST,/extract))
            tai      = double(strsplit(fibermap.TAI_LIST, /extract))
            weights  = double(strsplit(fibermap.FIELDSNR2I_LIST, /extract))
            airmass  = double(strsplit(fibermap.AIRMASS_LIST, /extract))
            SEEING20 = double(strsplit(fibermap.SEEING20_LIST,/extract))
            SEEING50 = double(strsplit(fibermap.SEEING50_LIST,/extract))
            SEEING80 = double(strsplit(fibermap.SEEING80_LIST,/extract))
            junk = where(strmatch(TAG_NAMES(fibermap), 'RMSOFF20_LIST',/fold_case), ct)
            if ct then begin
                RMSOFF20 = double(strsplit(fibermap.RMSOFF20_LIST,/extract))
                RMSOFF50 = double(strsplit(fibermap.RMSOFF50_LIST,/extract))
                RMSOFF80 = double(strsplit(fibermap.RMSOFF80_LIST,/extract))
            endif else begin
                RMSOFF20 = 0.d
                RMSOFF50 = 0.d
                RMSOFF80 = 0.d
                ntarget  = n_elements(fibermap)
                rms_tags=replicate(create_struct('RMSOFF20_LIST',' ', 'RMSOFF20',0.d, $
                                                 'RMSOFF50_LIST',' ', 'RMSOFF50',0.d,$
                                                 'RMSOFF80_LIST',' ', 'RMSOFF80',0.d),ntarget)
                fibermap=struct_addtags(fibermap,rms_tags)
            endelse
            firstcarton = strsplit(fibermap.FIRSTCARTON_LIST,/extract)
            cartpk = strsplit(fibermap.CARTON_TO_TARGET_PK_LIST,/extract)

            test=where((Assigned eq 1) AND (valid eq 1) AND (On_target eq 1), ctv)
            EXP_DISP_MED = -1.d
            if ctv ne 0 then begin
                if nexp gt 1 then begin
                    EXP_DISP_MED = abs(median(STDDEV(flux,DIMENSION=2,/double)/combinedflux,/EVEN,/double))
                endif else begin
                    EXP_DISP_MED = 0.d
                endelse
            endif
        

            fibermap.TARGET_INDEX = target_idx+1
            fibermap.NEXP = nexp
            fibermap.EXP_DISP_MED = EXP_DISP_MED
            
            if total(weights,/DOUBLE) gt 0 then weights = DBLARR(n_elements(weights)) +1
            fibermap.MJD_FINAL = total(tai/(24.D*3600.D)*weights,/DOUBLE)/total(weights,/DOUBLE)
            fibermap.AIRMASS   = total(airmass *weights,/DOUBLE)/total(weights,/DOUBLE)
            bighdr = clearhdrcard(bighdr, 'AIRMASS', value = (fibermap[0].AIRMASS))
            fibermap.SEEING20  = total(SEEING20*weights,/DOUBLE)/total(weights,/DOUBLE)
            bighdr = clearhdrcard(bighdr, 'SEEING20', value = fibermap[0].SEEING20)
            fibermap.SEEING50  = total(SEEING50*weights,/DOUBLE)/total(weights,/DOUBLE)
            bighdr = clearhdrcard(bighdr, 'SEEING50', value = fibermap[0].SEEING50)
            fibermap.SEEING80  = total(SEEING80*weights,/DOUBLE)/total(weights,/DOUBLE)
            bighdr = clearhdrcard(bighdr, 'SEEING80', value = fibermap[0].SEEING80)
            fibermap.RMSOFF20  = total(RMSOFF20*weights,/DOUBLE)/total(weights,/DOUBLE)
            bighdr = clearhdrcard(bighdr, 'RMSOFF20', value = fibermap[0].RMSOFF20)
            fibermap.RMSOFF50  = total(RMSOFF50*weights,/DOUBLE)/total(weights,/DOUBLE)
            bighdr = clearhdrcard(bighdr, 'RMSOFF50', value = fibermap[0].RMSOFF50)
            fibermap.RMSOFF80  = total(RMSOFF80*weights,/DOUBLE)/total(weights,/DOUBLE)
            bighdr = clearhdrcard(bighdr, 'RMSOFF80', value = fibermap[0].RMSOFF80)
       
            fibermap.CATALOGID = max(targ.CATALOGID_LIST)
            fibermap.iCATALOGID = LONG64(max(targ.CATALOGID_LIST))
            fibermap=struct_addtags(fibermap, replicate(create_struct('OBS',strjoin(obs[UNIQ(obs, SORT(obs))])),n_elements(fibermap)))
            fibermap.FIRSTCARTON_LIST = strjoin(firstcarton[UNIQ(firstcarton, SORT(firstcarton))], ' ')
            fibermap.CARTON_TO_TARGET_PK_LIST = strjoin(cartpk[UNIQ(cartpk, SORT(cartpk))], ' ')
            cid = max(targ.CATALOGID_LIST)
                        
            coaddname = 'spSpec-'+custom+'-'+strtrim(targ.EPOCH_COMBINE,2)+'-'+strtrim(cid,2)+'.fits'

            coaddname = filepath(coaddname, root_dir=custom_dir, subdirectory=['coadd',strtrim(targ.EPOCH_COMBINE,2)])
            spawn,'mkdir -p '+filepath('coadd',root_dir=custom_dir)
            spawn,'mkdir -p '+filepath(strtrim(targ.EPOCH_COMBINE,2),root_dir=custom_dir,subdirectory='coadd')

            sxaddpar, bighdr, 'VERSCOMB', idlspec2d_version(), $
                    ' Version of idlspec2d for combining multiple spectra', after='VERS2D'
            sxaddpar, bighdr, 'NEXP', nfiles, $
                    ' Number of exposures in this file', before='EXPID01'

            sxaddpar, bighdr, 'EXPTIME', min(exptimevec), $
                    ' Minimum of exposure times for all cameras'
       
            sxaddpar, bighdr, 'OBSERVATORY', strjoin(obs[UNIQ(obs, SORT(obs))],','),$
                    ' Observatory of observations'
            
            if strmatch(strjoin(obs[UNIQ(obs, SORT(obs))],','), 'APO', /FOLD_CASE) then begin
                sxaddpar, bighdr, 'SPEC', 'SP1'
                sxaddpar, bighdr, 'TELESCOP', 'SDSS 2.5-M'
            endif else begin
                if strmatch(strjoin(obs[UNIQ(obs, SORT(obs))],','), 'APO', /FOLD_CASE) then begin
                    sxaddpar, bighdr, 'SPEC', 'SP2'
                    sxaddpar, bighdr, 'TELESCOP', 'LCO Du Pont'
                endif else begin
                    sxaddpar, bighdr, 'SPEC', 'SP1+SP2'
                endelse
            endelse
            sxdelpar, bighdr, 'OBSMODE'
            sxaddpar, bighdr, 'MJD', LONG(max(mjds))
            sxaddpar, bighdr, 'TAI', mean(tai,/double)
            sxaddpar, bighdr, 'DATE-OBS', sxpar(*hdrs[0],'DATE-OBS')

            sxaddpar, bighdr, 'EXPTIME', max(exptimevec)

            foreach cam, ['B1','R1','B2','R2'], idx do begin
                icam = where(strmatch(camnames, cam, /fold_case),ct)
                if ct eq 0 then nexpcam = 0 else nexpcam = nexpvec[icam]
                sxdelpar, bighdr, 'NEXP_'+cam
                sxaddpar, bighdr, 'NEXP_'+cam, nexpcam, $
                        ' '+cam+' camera number of exposures', before='EXPTIME'
            endforeach
            foreach cam, ['B1','R1','B2','R2'], idx do begin
                icam = where(strmatch(camnames, cam, /fold_case),ct)
                if ct eq 0 then exptimecam = 0.0 else exptimecam = exptimevec[icam]
                sxdelpar, bighdr, 'EXPT_'+cam
                sxaddpar, bighdr, 'EXPT_'+cam, exptimecam, $
                        ' '+cam+' camera exposure time (seconds)', before='EXPTIME'
            endforeach

            sxdelpar, bighdr, 'SPCOADD'
            sxaddpar, bighdr, 'SPCOADD', systime(), ' SPCOADD finished', after='EXPTIME'
            sxaddpar, bighdr, 'NWORDER', 2, ' Linear-log10 coefficients'
            sxaddpar, bighdr, 'NWORDER', 2, ' Linear-log10 coefficients'
            sxaddpar, bighdr, 'WFITTYPE', 'LOG-LINEAR', ' Linear-log10 dispersion'
            sxaddpar, bighdr, 'COEFF0', wavemin, ' Central wavelength (log10) of first pixel'
            sxaddpar, bighdr, 'COEFF1', binsz, ' Log10 dispersion per pixel'

            sxdelpar, bighdr, 'EXTNAME'
            merge_spechdrmodel, hdr=bighdr

            sxdelpar, bighdr, 'GUIDER1'
            sxdelpar, bighdr, 'GUIDERN'
            sxdelpar, bighdr, 'NGUIDE'
            sxdelpar, bighdr, 'EXPIDXXX'
            
            finalvalues = create_struct('FLUX',0.0, 'LOGLAM',0.0, 'IVAR',0.0, $
                                        'AND_MASK', long(0), 'OR_MASK', long(0), $
                                        'WDISP', 0.0, 'SKY', 0.0, 'WRESL', 0.0)
            finalvalues = replicate(finalvalues, n_elements(finalwave))
            finalvalues.FLUX = combinedflux
            finalvalues.LOGLAM = finalwave
            finalvalues.IVAR = combinedivar
            finalvalues.AND_MASK = bestandmask
            finalvalues.OR_MASK = bestormask
            finalvalues.WDISP = bestwdisp
            finalvalues.SKY = bestsky
            finalvalues.WRESL = bestresolution


            ; HDU # 0 header
            mwrfits_named, junk_d, coaddname, hdr=bighdr, /create, /silent
       
            ; HDU # 1 header
            mwrfits_named, finalvalues, coaddname, hdr=coadd_hdr, name= 'COADD', desc=' Coadded spectrum', /silent
            coadd_hdr = 0
            ;sxdelpar, coadd_hdr, 'COMMENT'
              
            ; HDU #2 is plugmap
            mwrfits_named, fibermap, coaddname, hdr=hdrplug, name='PLUGMAP', desc=' Plugmap structure', /silent
            sxdelpar, hdrplug, 'COMMENT'
            
            for i = 0, (fibermap.NEXP)-1 do begin
                mwrfits_named, *(exps[i]), coaddname, hdr=*(hdrs[i]), /SILENT
                sxdelpar, hdr, 'COMMENT'
            endfor
        endforeach
        splog, target_idx , ' of ',nspSpec, ' Targets Merged'

        coaddhdr = 0
        coadd = custom
        
        spSpecfiles = findfile(filepath('spSpec*.fits', root_dir=get_field_Dir(topdir, '', custom, /custom),$
                                        subdirectory=['coadd', strtrim(ec,2)]), count = nspSpec)
        
        if nspSpec le 1 then begin
            splog, 'Skipping '+strtrim(ec,2)+' with '+strtrim(nspSpec,2)+' targets'
            continue
        endif
        spSpec2spFullSky, coadd, topdir=topdir, mjd=ec, runmjd=runmjd, outdir=custom_dir, coaddhdr = coaddhdr
        if keyword_set(coaddhdr) then begin


            foreach spf, spSpecfiles, idx do begin
                ;read spSpec file
                splog, idx+1 , ' of ',nspSpec, ' spSpec header Modified'
                fits_open,spf,io,/update    ;Faster to explicity open
                FITS_READ, io, data, hdr, EXTEN_NO=0, /HEADER_ONLY

                SPCOADD_card_idx = where(strmatch(hdr, 'SPCOADD*',/fold_case), ct)
                if ct gt 0 then begin
                    SPCOADD_card_idx = SPCOADD_card_idx[0]
                    hdr = [hdr[0:SPCOADD_card_idx], hdr[-1]]
                endif
                
                
                
                sxaddpar, hdr, "SPEC1_G",  sxpar(coaddhdr,"SPEC1_G"),  "(S/N)^2 for spec 1 at mag 21.20",  before='END'
                sxaddpar, hdr, "FSPEC1_G", sxpar(coaddhdr,"FSPEC1_G"), "Fit (S/N)^2 for spec  1 at mag 21.20",  before='END'
                sxaddpar, hdr, "SN2EXT1G", sxpar(coaddhdr,"SN2EXT1G"), "Extinction corrected (S/N)^2",  before='END'
                sxaddpar, hdr, "FSN2EX1G", sxpar(coaddhdr,"FSN2EX1G"), "Extinction corrected Fit (S/N)^2",  before='END'
                sxaddpar, hdr, "SPEC1_R",  sxpar(coaddhdr,"SPEC1_R"),  "(S/N)^2 for spec  1 at mag 20.20",  before='END'
                sxaddpar, hdr, "FSPEC1_R", sxpar(coaddhdr,"FSPEC1_R"), "Fit (S/N)^2 for spec  1 at mag 20.20",  before='END'
                sxaddpar, hdr, "SN2EXT1R", sxpar(coaddhdr,"SN2EXT1R"), "Extinction corrected (S/N)^2",  before='END'
                sxaddpar, hdr, "FSN2EX1R", sxpar(coaddhdr,"FSN2EX1R"), "Extinction corrected Fit (S/N)^2",  before='END'
                sxaddpar, hdr, "SPEC1_I",  sxpar(coaddhdr,"SPEC1_I"),  "(S/N)^2 for spec  1 at mag 20.20",  before='END'
                sxaddpar, hdr, "FSPEC1_I", sxpar(coaddhdr,"FSPEC1_I"), "Fit (S/N)^2 for spec  1 at mag 20.20",  before='END'
                sxaddpar, hdr, "SN2EXT1I", sxpar(coaddhdr,"SN2EXT1I"), "Extinction corrected (S/N)^2",  before='END'
                sxaddpar, hdr, "FSN2EX1I", sxpar(coaddhdr,"FSN2EX1I"), "Extinction corrected Fit (S/N)^2",  before='END'
                sxaddpar, hdr, "SPEC2_G",  sxpar(coaddhdr,"SPEC2_G"),  "(S/N)^2 for spec  2 at mag 21.20",  before='END'
                sxaddpar, hdr, "FSPEC2_G", sxpar(coaddhdr,"FSPEC2_G"), "Fit (S/N)^2 for spec  2 at mag 21.20",  before='END'
                sxaddpar, hdr, "SN2EXT2G", sxpar(coaddhdr,"SN2EXT2G"), "Extinction corrected (S/N)^2",  before='END'
                sxaddpar, hdr, "FSN2EX2G", sxpar(coaddhdr,"FSN2EX2G"), "Extinction corrected Fit (S/N)^2",  before='END'
                sxaddpar, hdr, "SPEC2_R",  sxpar(coaddhdr,"SPEC2_R"),  "(S/N)^2 for spec  2 at mag 20.20",  before='END'
                sxaddpar, hdr, "FSPEC2_R", sxpar(coaddhdr,"FSPEC2_R"), "Fit (S/N)^2 for spec  2 at mag 20.20",  before='END'
                sxaddpar, hdr, "SN2EXT2R", sxpar(coaddhdr,"SN2EXT2R"), "Extinction corrected (S/N)^2",  before='END'
                sxaddpar, hdr, "FSN2EX2R", sxpar(coaddhdr,"FSN2EX2R"), "Extinction corrected Fit (S/N)^2",  before='END'
                sxaddpar, hdr, "SPEC2_I",  sxpar(coaddhdr,"SPEC2_I"),  "(S/N)^2 for spec  2 at mag 20.20",  before='END'
                sxaddpar, hdr, "FSPEC2_I", sxpar(coaddhdr,"FSPEC2_I"), "Fit (S/N)^2 for spec  2 at mag 20.20",  before='END'
                sxaddpar, hdr, "SN2EXT2I", sxpar(coaddhdr,"SN2EXT2I"), "Extinction corrected (S/N)^2",  before='END'
                sxaddpar, hdr, "FSN2EX2I", sxpar(coaddhdr,"FSN2EX2I"), "Extinction corrected Fit (S/N)^2",  before='END'
                modfits, io, 0, hdr, exten_no = 0               ; update hdr
                
                FITS_CLOSE,io
            endforeach
        endif
        splog, 'Successful completion of spspec_target_merge for Epoch Combine '+strtrim(ec,2)+' at '+systime()
        if not keyword_set(mjd_flag) then splog, /close_secondary

    endforeach
   
    splog, 'Successful completion of spspec_target_merge at ' + systime()
    splog, /close_all
end
