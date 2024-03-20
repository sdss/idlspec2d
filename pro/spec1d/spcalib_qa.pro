; NAME:
;   spcalib_qa
;
; PURPOSE:
;   Compare photometric accuracy of standards
;
; CALLING SEQUENCE:
;   SpCalib_QA, [run2d=, fieldid=, mjd=]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   field       - field to include
;   mjd         - MJD to include
;   run2d       - RUN2D version of reduction
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Depends on the spAll files (either full run2d version or field-mjd version)
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; Function Called:
;   mpfitfun
;   field_to_string
;   djs_filepath
;   mrdfits
;   sdss_flagval
;
; External PROCEDURES CALLED:
;   plot
;   XYOUTS
;   cpbackup
;
; Internal PROCEDURES CALLED:
;   std_hist
;
; REVISION HISTORY:
;   21-June-2022  Written by S. Morrison (UIUC)
;------------------------------------------------------------------------------

pro std_hist, fratio, bs=bs, xmin=xmin, xmax=xmax, filt=filt, fit=fit

    if not keyword_set(bs) then bs=0.015
    if not keyword_set(xmin) then xmin=-.3
    if not keyword_set(xmax) then xmax=-.2
    if not keyword_set(filt) then filt='g'

    logf = alog10(fratio)
    junk = where((logf ge xmin) AND (logf le xmax),ct)
    if ct eq 0 then begin
        splog, "catastrophic failure of Flux calibration"
        fit = [!VALUES.D_NAN,!VALUES.D_NAN]
        return
    endif
    yhist=histogram(logf,locations=xhist,binsize=bs,omin=hmin,omax=hmax,min=xmin, max=xmax)
    p=mpfitfun('gauss1', xhist, yhist, yhist,weights=1.D, [0,0.2,1],/QUIET)
    xhist=[hmin-bs, xhist,hmax+bs/2.0]
    yhist=[0, yhist, 0]
    plot, xhist, yhist, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[0,max(yhist)+2], xticks=5, xminor=10, FONT=0, thick=2.75, psym=10, xtitle='log (synflux/calibflux) ['+filt+']', ytitle='Nstd'
    XYOUTS, xmin+0.02,(max(yhist)*.9), 'mean, sigma!C'+string(p[0],format='(F5.2)')+','+string(p[1],format='(F5.3)'),  CHARSIZE=1.2,  FONT=0
    xfit = findgen((xmax-xmin)/0.001, increment=0.001,start=xmin)
    fit = p
    oplot, xfit, GAUSS1(xfit, p), psym=0, linestyle=1
    return
end

pro SpCalib_QA, run2d=run2d, fieldid=fieldid, mjd=mjd, rerun=rerun, nobkup=nobkup, epoch=epoch
    RESOLVE_ALL, /QUIET, /SKIP_EXISTING, /CONTINUE_ON_ERROR

    if not keyword_set(run2d) then run2d = getenv('RUN2D')
    outname='spCalib_QA-'+run2d
    if keyword_set(fieldid) then begin
        spallfile='spAll-'+field_to_string(fieldid)+'-'+strtrim(mjd,2)+'.fits'
        if keyword_set(epoch) then begin
            spallfile = djs_filepath(spallfile, root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d+'/epoch/spectra/full/'+field_to_string(fieldid)+'/'+strtrim(mjd,2))
            out_csv = djs_filepath(outname+'.csv', root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d+'/epoch')
            outname = outname+'-'+field_to_string(fieldid)+'-'+strtrim(mjd,2)
            outname = djs_filepath(outname+'.ps', root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d+'/'+field_to_string(fieldid)+'/epoch')
        endif else begin
            spallfile = djs_filepath(spallfile, root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d+'/spectra/full/'+field_to_string(fieldid)+'/'+strtrim(mjd,2))
            out_csv = djs_filepath(outname+'.csv', root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d)
            outname = outname+'-'+field_to_string(fieldid)+'-'+strtrim(mjd,2)
            outname = djs_filepath(outname+'.ps', root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d+'/'+field_to_string(fieldid))
        endelse
        logfile = repstr(outname,'.ps','.log')
        
        if not keyword_set(nobkup) then cpbackup, logfile
        splog, filename=logfile
        splog, 'Log file ' + logfile + ' opened ' + systime()

        if not keyword_set(nobkup) then cpbackup, outname
    endif else begin
        spallfile = 'spAll-'+run2d+'.fits'
        
        if keyword_set(epoch) then begin
            spallfile = djs_filepath(spallfile, root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d+'/epoch')
        endif else begin
            spallfile = djs_filepath(spallfile, root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d)
        endelse
        if not keyword_set(rerun) then begin
            if keyword_set(epoch) then begin
                out_csv = djs_filepath(outname+'.csv', root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d+'/epoch')
                outname = djs_filepath(outname+'.ps', root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d+'/epoch')
            endif else begin
                out_csv = djs_filepath(outname+'.csv', root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d)
                outname = djs_filepath(outname+'.ps', root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2d)
            endelse
        endif else begin
            spall_tag=mrdfits(lookforgzip(spallfile),1,/silent)
            fieldids = spall_tag.field
            fieldids = fieldids[uniq(fieldids)]
            foreach fieldid, fieldids do begin
                mjds = spall_tag[where(spall_tag.field eq fieldid)].mjd
                mjds = mjds[uniq(mjds)]
                foreach mjd, mjds do begin
                    SpCalib_QA, run2d=run2d, fieldid=fieldid, mjd=mjd, nobkup=nobkup
                endforeach
            endforeach
            return
        endelse
    endelse
    
    spallfile = lookforgzip(spallfile)
    if file_test(spallfile,/ZERO_LENGTH) or ( not file_Test(spallfile,/read)) then begin
        wait, 5
        if file_test(spallfile,/ZERO_LENGTH) or ( not file_Test(spallfile,/read)) then begin
            splog,'skipping', spallfile
            return
        endif
    endif
    print,spallfile
    spall=mrdfits(spallfile,1,/silent)
    ind = where(strmatch(spall.objtype, 'SPECTROPHOTO_STD',/fold_case) $
            and ((spall.zwarning AND sdss_flagval('ZWARNING', 'UNPLUGGED')) eq 0), ct_std) ; std after removing unplugged ones
    if ct_std ne 0 then begin
        splog,spallfile, n_elements(ind)
        f_ratio = spall[ind].SPECTROSYNFLUX/spall[ind].calibflux
        ind_good = where(f_ratio[1,*] gt 0, ct_std)
    endif else begin
        splog,spallfile
        ct_std = 0
    endelse
    if ct_std ne 0 then begin
        mydevice = !D.NAME
        SET_PLOT, 'ps'
        DEVICE, FILENAME=outname,/LANDSCAPE, /times

        !p.multi=[0,1,3]
        !psym=10
        !y.margin=[4,4]
        !x.ticklen=0.075
        !X.thick=2.0
        !Y.thick=2.0
        !X.charsize=2
        !Y.charsize=2

        std_hist, f_ratio[1, ind_good], bs=0.015, xmin=-.3, xmax=.2, filt='g', fit=fit_g
        std_hist, f_ratio[2, ind_good], bs=0.015, xmin=-.3, xmax=.2, filt='r', fit=fit_r
        std_hist, f_ratio[3, ind_good], bs=0.015, xmin=-.3, xmax=.2, filt='i', fit=fit_i

        if keyword_set(fieldid) then XYOUTS, 0.5,0.98, 'Field='+field_to_string(fieldid)+' MJD='+strtrim(mjd,2),  CHARSIZE=1.2,  FONT=0, alignment=0.5,/NORMAL

        DEVICE, /CLOSE
        ps2pdf, outname
        ; Return plotting to the original device:
        SET_PLOT, mydevice
    endif else begin
        fit_g = [!VALUES.D_NAN,!VALUES.D_NAN]
        fit_r = [!VALUES.D_NAN,!VALUES.D_NAN]
        fit_i = [!VALUES.D_NAN,!VALUES.D_NAN]
    endelse
    if keyword_set(fieldid) then begin
        try = 0
        retry: try = try+1
        if tag_exist(spall,'OBS') then obs = spall[0].OBS else obs='APO'
        if file_test(out_csv) then begin
            outs = create_struct( $
                                'field' , 0L,   $
                                'mjd'   , 0L,   $
                                'obs'   , ' ',  $
                                'g_mean', 0.0d, $
                                'g_sig',  0.0d, $
                                'r_mean', 0.0d, $
                                'r_sig',  0.0d, $
                                'i_mean', 0.0d, $
                                'i_sig',  0.0d, $
                                'n_std',  0)

            nrows = File_Lines(out_csv) - 1
            if nrows eq 0 and try < 3 then begin
                wait,10
                splog, 'retrying ',out_csv
                goto, retry
            endif
            ins = Replicate(outs, nrows)
            temp_in = read_csv(out_csv)
            foreach tag, tag_names(ins), i do ins.(i) = temp_in.(i)
                        
            match = where(ins.field eq fieldid and ins.mjd eq mjd, ct)
            if ct eq 0 then begin
                outs = struct_append( ins, outs)
                match = nrows
            endif else begin
                outs = ins
            endelse
            outs[match].field  = fieldid
            outs[match].mjd    = mjd
            outs[match].obs    = obs
            outs[match].g_mean = fit_g[0]
            outs[match].g_sig  = fit_g[1]
            outs[match].r_mean = fit_r[0]
            outs[match].r_sig  = fit_r[1]
            outs[match].i_mean = fit_i[0]
            outs[match].i_sig  = fit_i[1]
            outs[match].n_std  = ct_std
        endif else begin
            outs = create_struct( $
                                'field' ,  fieldid,  $
                                'mjd'   ,  mjd,      $
                                'obs'   ,  obs,      $
                                'g_mean',  fit_g[0], $
                                'g_sig',   fit_g[1], $
                                'r_mean',  fit_r[0], $
                                'r_sig',   fit_r[1], $
                                'i_mean',  fit_i[0], $
                                'i_sig',   fit_i[1], $
                                'n_std',   ct_std)
        endelse
        WRITE_CSV, out_csv, outs, HEADER=tag_names(outs)
    endif
    splog, 'SpectroPhoto QA Complete'
    if (keyword_set(logfile)) then splog, /close
end
