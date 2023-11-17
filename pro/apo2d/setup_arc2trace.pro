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
        'b1': name = [filename_short,'']
        'r1': name = ['',filename_short]
        'b2': name = [filename_short,'']
        'r2': name = ['',filename_short]
    endcase
    spexp1 = { $
            confid : strtrim(sxpar(hdr, 'CONFID'),2), $
            fieldid: strtrim(sxpar(hdr, 'FIELDID'),2), $
            mjd    : LONG(sxpar(hdr, 'MJD')), $
            flavor : string(flavor), $
            exptime: float(sxpar(hdr, 'EXPTIME')), $
            name    :name }
    if keyword_set(spexp) then spexp = [spexp, spexp1] else spexp = spexp1
    return, spexp
end


pro setup_arc2trace, tsetfile, arcfile, indir, outdir, mjd, cam

    obs = GETENV("OBSERVATORY")

    ; search for plan file
    planfile = 'spPlanTrace-'+strtrim(mjd,2)+'_'+obs+'_'+cam+'.par'
    fullplanfile = djs_filepath(planfile, root_dir=outdir, subdirectory=['trace',strtrim(mjd,2)])
    while(djs_lockfile(fullplanfile) eq 0) do wait, 5
    plan = file_search(fullplanfile)
    
    tsetsplit = strsplit((strsplit(tsetfile,'.',/extract))[0],'-',/extract)
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
    
    junk = where(strmatch(names, FILE_BASENAME(flatfile)), ct)
    if ct eq 0 then $
        spexp = buildplan(flatfile, cam, flatflavor, indir=indir, spexp=spexp)

    junk = where(strmatch(names, FILE_BASENAME(arcfile)), ct)
    if ct eq 0 then $
        spexp = buildplan(arcfile,  cam, arcflavor,  spexp=spexp)
    


    
    if not keyword_set(plan) then begin
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
        
        flatinfoname = 'spTraceFlat-'+cam+'-'
        toutdir = djs_filepath('', root_dir=outdir, subdirectory=['trace',strtrim(mjd,2)])
        FILE_MKDIR, toutdir
        write_spflat, flatinfoname, 0, flatstruct, flathdr, [FILE_BASENAME(arcfile)], 0, 0, 0, $
                      outdir=toutdir
    endif

    hdr = ''
    hdr = [hdr, "MJD      " + strtrim(mjd,2)        + "  # Modified Julian Date"]
    hdr = [hdr, "OBS      " + obs                   + "  # Observatory"]
    hdr = [hdr, "RUN2D    " + GETENV("IDLSPEC2D_VER") + "  # 2D reduction name"]

    planfile = 'spPlanTrace-'+strtrim(mjd,2)+'_'+obs+'_'+cam+'.par'
    fullplanfile = djs_filepath(planfile, root_dir=outdir, subdirectory=['trace',strtrim(mjd,2)])

    yanny_write, fullplanfile, ptr_new(spexp), hdr=hdr, stnames='SPEXP'
    djs_unlockfile, fullplanfile
end
;------------------------------------------------------------------------------
