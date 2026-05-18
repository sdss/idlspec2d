
pro spbuild_traceflat, plan2d=plan2d, mjd=mjd, obs=obs, flat=flat, arc=arc, debug=debug, saveraw=saveraw, hart=hart
    CPU, TPOOL_NTHREADS = 1
    if not keyword_set(plan2d) then plan2d = ['SOS']
    print,plan2d
    foreach thisplan, plan2d do begin
        if not strmatch(plan2d, 'SOS', /fold_case) then begin
            allseq = yanny_readone(thisplan, 'SPEXP', hdr=hdr, /anon)
            print, allseq
            obs =  yanny_par(hdr, 'OBS')
            mjd = yanny_par(hdr, 'MJD')
            if mjd lt 59146 then begin
                legacy = 1
                plates = 0
            endif else begin
                if mjd lt 59564 then begin
                    legacy = 0
                    plates = 1
                endif else begin
                    legacy = 0
                    plates = 0
                endelse
            endelse

            if strmatch(obs,'LCO',/fold_Case) then rawdir=getenv('BOSS_SPECTRO_DATA_S') $
                    else rawdir=getenv('BOSS_SPECTRO_DATA_N')
                
            if keyword_set(legacy) then ccds=['b1','b2','r1','r2'] else begin
                if strmatch(obs,'LCO',/fold_Case) then ccds=['b2','r2'] else ccds=['b1','r1']
            endelse
            logfile_base = repstr(repstr(plan2d, 'spPlan2d', 'spTrace'),'.par','')
            logfile = logfile_base+'.log'
            psfile  = logfile_base+'.ps'

            stime0 = systime(1)
            cpbackup, logfile
            splog, filename=logfile
            splog, 'Log file '+ logfile + ' opened ' + systime()
            cpbackup, psfile
            set_plot, 'ps'
            device, filename=psfile, /color
            splog, 'Plot file ' + psfile
        endif else begin
            rawdir = '/data/spectro'
            ccds = ['SOS']
        endelse
        splog, 'CCDs', ccds

        rawdir = filepath(strtrim(mjd,2), root_dir=rawdir)
        config_dir = filepath('', root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
        ;ecalibfile = findopfile('opECalib*par', mjd, config_dir, /abort_notfound, silent=silent)
        
        fieldid = long(yanny_par(hdr, 'fieldname'))
        foreach ccd, ccds, icam do begin
            if not strmatch(plan2d, 'SOS', /fold_case) then begin
                j = where(allseq.flavor EQ 'TRACEFLAT' $
                          AND allseq.name[icam] NE 'UNKNOWN', nflat )
                flatname = allseq[j].name[icam]
                if not keyword_set(hart) then begin
                    k = where(allseq.flavor EQ 'TRACEARC' $
                              AND allseq.name[icam] NE 'UNKNOWN', narc )
                endif else begin
                    k = where(allseq.flavor EQ 'TRACEARC' $
                              AND allseq.name[icam] NE 'UNKNOWN', narc )
                    if narc eq 0 then begin
                        k = where(allseq.flavor EQ 'HART' $
                              AND allseq.name[icam] NE 'UNKNOWN', narc )
                    endif
                endelse
                arcname = allseq[k].name[icam]
                flatinfoname = 'spTraceFlat-'+ccd+'-'
                arcinfoname = 'spTraceArc-'+ccd+'-'
            endif else begin
                flatinfoname = repstr(repstr( flat, 'sdr', 'spTraceFlat'), '.fit', '.fits')
                arcinfoname  = repstr(repstr( arc,  'sdr', 'spTraceArc'),  '.fit', '.fits')
            endelse

            if FILE_TEST(filepath(flatname+'.gz', root_dir=rawdir)) then begin
                junk = mrdfits(filepath(flatname+'.gz', root_dir=rawdir),0, framehdr, /silent)
            endif else junk = mrdfits(filepath(flatname, root_dir=rawdir),0, framehdr, /silent)
            cartid=sxpar(framehdr, 'cartid')
            plottitle = repstr(flatinfoname,'.fits','.ps')

            plottitle = ' FIELDID='+strtrim(yanny_par(hdr, 'fieldname'),2) +' '+ ' MJD='+strtrim(mjd,2)+' '
            if narc eq 1 then begin
                spcalib, flatname, arcname, cartid=cartid, indir=rawdir, ecalibfile=ecalibfile, $
                        plottitle=plottitle, flatinfoname=flatinfoname, arcinfoname=arcinfoname,$
                        plates = plates, legacy=legacy, timesep=0, saveraw = saveraw, debug=debug, /buildTraceFlat
            endif else begin
                foreach arcname1, arcname, iarc do begin
                    splog, 'Running ', arcname1, ' ', arcinfoname
                    spcalib, flatname, arcname1, cartid=cartid, indir=rawdir, ecalibfile=ecalibfile, $
                            plottitle=plottitle, flatinfoname=flatinfoname, arcinfoname=arcinfoname,$
                            plates = plates, legacy=legacy, timesep=0, saveraw = saveraw, debug=debug, /buildTraceFlat
                endforeach
            endelse
        endforeach
        device, /close
        set_plot, 'x'
        ps2pdf, psfile
        splog, 'Total time for run_spcalib = ', systime(1)-stime0, ' seconds', format='(a,f6.0,a)'
        splog, 'Successful completion of run_spcalib at ' + systime()
        splog, /close
    endforeach
end
