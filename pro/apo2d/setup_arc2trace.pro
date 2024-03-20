;+
; NAME:
;   setup_arc2trace
;
; PURPOSE:
;   Setup files for arc_to_trace
;
; CALLING SEQUENCE:
;   setup_arc2trace, tsetfile, arcfile, indir, outdir
;
; INPUTS:
;   tsetfile - SOS traceset file
;   arcfile  - raw arc lamp calibration frame
;   indir    - location of raw exposure files
;   outdir   - current sosdir
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   buildplan
;   sdsshead
;   yanny_read
;   MRDFITS
;   traceset2xy
;   create_struct
;   write_spflat
;   yanny_write
;   djs_filepath
;
; REVISION HISTORY:
;   Written S.Morrison Nov 15, 2022
;-
;------------------------------------------------------------------------------
function buildplan, filename, cam, flavor, indir=indir, spexp=spexp
    if keyword_set(indir) then $
        filename = djs_filepath(filename, root_dir=indir)
    hdr = sdsshead(filename, do_lock=do_lock)
    filename_short = FILE_BASENAME(filename)
    filename_short = repstr(filename_short,'.gz','')
    
    CASE cam OF
        'b1': name = [filename_short,repstr(filename_short, 'b1','r1')]
        'r1': name = [repstr(filename_short, 'r1','b1'),filename_short]
        'b2': name = [filename_short,repstr(filename_short, 'b2','r2')]
        'r2': name = [repstr(filename_short, 'r2','b2'),filename_short]
    endcase
    
    spexp1 = { $
            confid : strtrim(sxpar(hdr, 'CONFID'),2), $
            fieldid: strtrim(string(sxpar(hdr, 'FIELDID'),f='(i6.6)'),2), $
            mjd    : LONG(sxpar(hdr, 'MJD')), $
            flavor : string(flavor), $
            exptime: float(sxpar(hdr, 'EXPTIME')), $
            name    :name }
    if keyword_set(spexp) then spexp = [spexp, spexp1] else spexp = spexp1
    return, spexp
end


pro setup_arc2trace, tsetfile, fflatfile, arcfile, indir, outdir, mjd, cam, fieldid

    obs = GETENV("OBSERVATORY")

    ; search for plan file
    planfile = 'spPlanTrace-'+strtrim(mjd,2)+'_'+obs+'.par'
    fullplanfile = djs_filepath(planfile, root_dir=outdir, subdirectory=['trace',strtrim(mjd,2)])
    file_dir = djs_filepath(strtrim(mjd,2), root_dir =outdir, subdirectory=['trace'] )
    if not file_test(file_dir,/DIRECTORY) then FILE_MKDIR, file_dir
    
    
    while(djs_lockfile(fullplanfile+'.'+cam) eq 0) do wait, 5
    i = 0
    while(djs_lockfile(fullplanfile) eq 0) do begin
        i = i+1
        if strmatch(cam, 'b?') then begin
            altcam = repstr(cam,'b','r')
            af = fullplanfile+'.'+altcam
            if i ge 4 then begin
                if (file_test(af,/dangling_symlink)+file_test(af,/symlink)) gt 0 then begin
                    splog, 'file locked by both cam process for > 20s'
                    break
                endif
            endif
        endif
        wait, 5
    endwhile
    plan = file_search(fullplanfile)
    splog, 'test0'

    tsetsplit = strsplit((strsplit((file_basename(tsetfile)),'.',/extract))[0],'-',/extract)
    ;tsetsplit = strsplit((strsplit(tsetfile,'.',/extract))[0],'-',/extract)
    flatfile = 'sdR-'+tsetsplit[-1]+'-'+tsetsplit[-2]+'.fit'
    
    if keyword_set(plan) then begin
        arcflavor = 'arc'
        flatflavor = 'flat'
        yanny_read, plan, pp, hdr=hdr
        spexp = *pp
        names = spexp.name

    endif else begin
        arcflavor  = 'TRACEARC'
        flatflavor = 'TRACEFLAT'
        names = ['']
    endelse
    
    splog, 'test1'
    junk = where(strmatch(names, repstr(FILE_BASENAME(flatfile),'.gz','')), ct)
    if ct eq 0 then begin
        spexp = buildplan(flatfile, cam, flatflavor, indir=indir, spexp=spexp)
    endif else begin
        flatfile = djs_filepath(flatfile, root_dir=indir)
    endelse

    junk = where(strmatch(names, repstr(FILE_BASENAME(arcfile),'.gz','')), ct)
    if ct eq 0 then $
        spexp = buildplan(arcfile,  cam, arcflavor,  spexp=spexp)
    
    splog, 'test2'

    flatinfoname = 'spTraceFlat-'+cam+'-'
    if not file_test(djs_filepath(flatinfoname+'*', root_dir=file_dir) )then begin
            splog, 'test3a'
        ; convert tset to spflat format
        flathdr = sdsshead(flatfile, do_lock=do_lock)

        tset = MRDFITS(tsetfile,2)
        fibermask = MRDFITS(tsetfile,5)
        traceset2xy, tset, ycen, xsol

        ftemp = create_struct( name='FLAT_STRUCT', $
                              'NAME', FILE_BASENAME(flatfile), $
                              'IARC', 0, 'PROFTYPE', 1, $
                              'MEDWIDTH', fltarr(4), 'FIBERMASK', ptr_new(fibermask), $
                              'TSET', ptr_new(tset), 'XSOL', ptr_new(xsol), $
                              'WIDTHSET', ptr_new(0), 'FFLAT', ptr_new(0), $
                              'SUPERFLATSET', ptr_new(0), 'HDR', ptr_new(flathdr))
        flatstruct = replicate(ftemp, 1)
        
        
        toutdir = djs_filepath('', root_dir=outdir, subdirectory=['trace',strtrim(mjd,2)])
        FILE_MKDIR, toutdir
        write_spflat, flatinfoname, 0, flatstruct, flathdr, [FILE_BASENAME(arcfile)], 0, 0, 0, $
                      outdir=toutdir
    endif else begin
            splog, 'test3b'
        flatinfoname1 = 'spFlat-'+cam+'-'
        ; convert tset to spflat format
        flathdr = sdsshead(flatfile, do_lock=do_lock)
            splog, 'test3b1'
        splog, tsetfile
        tset = MRDFITS(tsetfile,2)
                    splog, 'test3b2'
        fibermask = MRDFITS(tsetfile,5)
                splog, 'test3b3'
        traceset2xy, tset, ycen, xsol
                splog, 'test3b4'
        splog, fflatfile
        fflat = MRDFITS(fflatfile, 0)
                splog, 'test3b5'

        ftemp = create_struct( name='FLAT_STRUCT', $
                              'NAME', FILE_BASENAME(flatfile), $
                              'IARC', 0, 'PROFTYPE', 1, $
                              'MEDWIDTH', fltarr(4), 'FIBERMASK', ptr_new(fibermask), $
                              'TSET', ptr_new(tset), 'XSOL', ptr_new(xsol), $
                              'WIDTHSET', ptr_new(0), 'FFLAT', ptr_new(fflat), $
                              'SUPERFLATSET', ptr_new(0), 'HDR', ptr_new(flathdr))
        flatstruct = replicate(ftemp, 1)
                splog, 'test3b6'

        toutdir = djs_filepath('', root_dir=outdir, subdirectory=[fieldid])
        ;toutdir = djs_filepath('', root_dir=outdir, subdirectory=[tsetsplit[-3]])
        FILE_MKDIR, toutdir
        write_spflat, flatinfoname1, 0, flatstruct, flathdr, [FILE_BASENAME(arcfile)], 0, 0, 0, $
                      outdir=toutdir
    

    endelse
    splog, 'test4'

    

    hdr = ''
    hdr = [hdr, "MJD      " + strtrim(mjd,2)        + "  # Modified Julian Date"]
    hdr = [hdr, "OBS      " + obs                   + "  # Observatory"]
    hdr = [hdr, "RUN2D    " + GETENV("IDLSPEC2D_VER") + "  # 2D reduction name"]

    yanny_write, fullplanfile, ptr_new(spexp), hdr=hdr, stnames='SPEXP'
    djs_unlockfile, fullplanfile
    djs_unlockfile, fullplanfile+'.'+cam
end
;------------------------------------------------------------------------------
