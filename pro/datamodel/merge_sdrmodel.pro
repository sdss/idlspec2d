




pro merge_sdrmodel, data=data, hdr=hdr

    hdr_model = yanny_readone(filepath('sdr_hdr.par',root_dir=getenv('IDLSPEC2D_DIR'), subdir='datamodel'))

    dropcards = ['CAMROW','CAMCOL','VERSIDL','VERSUTIL','VERSREAD','VERSLOG','VERSFLAT','TWOPHASE','QUALITY']


    if isa(sxpar(hdr,'CARTID'), /NUMBER) then begin
        FPS = 0
        obs='APO'
        platetype = sxpar(hdr, 'PLATETYP', count=nhdr)
        if nhdr eq 0 then sxaddpar, hdr, 'PLATETYP', 'BOSS',  before='END'
        
    endif else begin
        FPS = 1
        if strmatch(sxpar(hdr,'CARTID'), '*FPS-S*', /fold_case) then obs='LCO' else obs='APO'
        sxaddpar, hdr, 'PLATETYP', 'BHM&MWM',  before='END'
    endelse

    cam = strtrim(SXPAR(hdr, 'CAMERAS'),2)
    sdr = fix(cam.Extract('[1-9]{1}'))
    

    if keyword_set(FPS) then begin
        if sdr eq 2 then LCO = 1
    endif else LCO = 0

    mergedHdr = ['']

    if keyword_set(data) then $
        mkhdr, mergedHdr, data $
    else  mkhdr, mergedHdr, ''
    sxdelpar, mergedHdr, ['COMMENT','DATE']
    sxdelpar, mergedHdr, ['CAMROW','CAMCOL','VERSIDL','VERSUTIL','VERSREAD','VERSLOG','VERSFLAT','TWOPHASE','QUALITY']

    
    ; Copy from raw to output header
 
    group = 'Primary'
    for i=0, n_elements(hdr_model)-1 do begin
        value = 0
        if not strmatch(hdr_model[i].Group, '*'+group+'*', /fold_case) then begin
            group = hdr_model[i].Group
            if not strmatch(group, 'Secondary', /fold_case) then $
                mergedHdr = [mergedHdr[0:-2], '', '        '+ STRUPCASE(group),mergedHdr[-1]]
        endif
        
        if strmatch(group,'Primary', /fold_case) then continue
        
        ct = 0
        if sdr eq 1 then begin
            if strlen(hdr_model[i].CopySP1) gt 0 then Value = SXPAR(hdr, hdr_model[i].CopySP1, count = ct)
        endif else begin
            if strlen(hdr_model[i].CopySP2) gt 0 then Value = SXPAR(hdr, hdr_model[i].CopySP2, count = ct)
        endelse
        if ct eq 0 then Value = SXPAR(hdr, hdr_model[i].card, count = ct)
        if ct eq 0 then begin
            if keyword_set(LCO) then begin
                Value = hdr_model[i].LCO_Default
            endif else Value = hdr_model[i].APO_Default
        endif
        
        
        if isa(Value, 'byte') then begin
            if Value eq 1 then Value = 'T' else Value = 'F'
        endif
        sxaddpar, mergedHdr, hdr_model[i].card, Value, hdr_model[i].Comment, before='END'
    endfor
    
    ; Check datamodel for missing cards and report
    for i = 0, n_elements(hdr)-1 do begin
        if n_elements(strsplit(hdr[i],'=',/extract)) eq 1 then continue
        card = strtrim((strsplit(hdr[i],'=',/extract))[0],2)
        junk = where(strmatch(hdr_model.Card, card,/fold_case), ct)
        if ct > 0 then continue

        if sdr eq 1 then begin
            junk = where(strmatch(hdr_model.CopySP1, card,/fold_case), ct)
            if ct > 0 then continue
        endif

        if sdr eq 2 then begin
            junk = where(strmatch(hdr_model.CopySP2, card,/fold_case), ct)
            if ct > 0 then continue
        endif
        
        junk = where(strmatch(['CHECKSUM','DATASUM','GAIN','RDNOISE', 'END'], card, /fold_case), ct)
        if ct > 0 then continue
        
        junk = where(strmatch(dropcards,card,/fold_case), ct)
        if ct eq 0 then $
            splog, card+' missing from datamodel - Not propagated to outputs'
    endfor
    hdr = mergedHdr
    hdr = [hdr[0:-2], '', '        '+'PIPELINE OUTPUTS',mergedHdr[-1]]
    
end
