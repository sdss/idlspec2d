
pro merge_spechdrmodel, hdr=hdr, drop=drop

    hdr_model = yanny_readone(filepath('spec_hdr.par',root_dir=getenv('IDLSPEC2D_DIR'), subdir='datamodel'))
    dropcards = ['CCD','CCDID','CCDTYPE','FLAVOR','INTSTART','INTEND',$
                 'OFFRA','OFFDEC','OFFPA','GSEEING','HA','ROTTYPE','ROTPOS',$
                 'FOCUS','SCALE','OBJSYS','BOREOFFX','BOREOFFY',$
                 'ARCOFFX','ARCOFFY','CALOFFX','CALOFFY','CALOFFR',$
                 'GUIDOFFX','GUIDOFFY','GUIDOFFR',$
                 'M2PISTON','M2XTILT','M2YTILT','M2XTRAN','M2YTRAN','M2ZROT',$
                 'M1PISTON','M1XTILT','M1YTILT','M1XTRAN','M1YTRAN','M1ZROT',$
                 'HEAR','HARTMANN','UNAME','DEWPOINT','DUSTA','DUSTB','GUSTS',$
                 'HUMIDITY','PRESSURE','WINDD','WINDS','PLUGFILE','TILEID',$
                 'WTIME','DEWDEP','DUSTC','DUSTD','HUMIDOUT',$
                 'TEMP01','TEMP02','TEMP03','TEMP04','XCHI2',$
                 'XCHI2MAX','XCHI2MIN','SFLATMIN','TRACFLAT'$
                 ]
    if keyword_set(drop) then dropcards = [dropcards, drop]
    mergedHdr = [' ']
    
    mkhdr, mergedHdr, ''
    
    endidx = where(strmatch(mergedHdr, 'END',/fold_case),ct)
    mergedHdr = mergedHdr[0:endidx[0]-1]
    ; Copy from spFrame to output header
 
    group = 'Primary'
    for i=0, n_elements(hdr_model)-1 do begin
        value = 0
        if not strmatch(hdr_model[i].Group, '*'+group+'*', /fold_case) then begin
            group = hdr_model[i].Group
            if not strmatch(group, 'Secondary', /fold_case) then begin
                mergedHdr = [mergedHdr[0:-2], '', '        '+ STRUPCASE(group),mergedHdr[-1]]
            endif
        endif
        
        if strmatch(group,'Primary', /fold_case) then continue
        
        
        if strmatch(hdr_model[i].card, 'EXPIDXX', /fold_case) then begin
            idx = where(strmatch(hdr, 'EXPID*', /fold_case), ct)
            if ct gt 0 then begin
                foreach j, idx do begin
                    tcard = strtrim(strsplit(hdr[j],'=',/extract),2)
                    comment = strsplit(hdr[j],'/',/extract)
                    if n_elements(comment) ge 2 then comment = comment[1] else comment = ''
                    sxaddpar, mergedHdr, tcard, sxpar(hdr, tcard), comment, before='END'
                endforeach
            endif
            continue
        endif
        
        
        loopset= 0
        foreach cc, ['SIGBSXXX','CENBSXXX','AVGBSXXX','STDBSXXX',$
                     'SIGASXXX','CENASXXX','AVGASXXX','STDASXXX'] do begin
            if strmatch(hdr_model[i].card, cc, /fold_case) then begin

                idx = where(strmatch(hdr, STRMID(cc,0,strlen(cc)-3)+'*',/fold_case),ct)
                if ct gt 0 then begin
                    foreach j, idx do begin
                        tcard = strtrim((strsplit(hdr[j],'=',/extract))[0],2)
                        comment = (strsplit(hdr[j],'/',/extract))[-1]
                        if n_elements(comment) ge 2 then comment = comment[1] else comment = ''
                        sxaddpar, mergedHdr, tcard, sxpar(hdr, tcard), comment, before='END'
                    endforeach
                endif
                loopset =1
                break
            endif
        endforeach
        if keyword_set(loopset) then continue
        
        
        
        ct = 0
        if ct eq 0 then Value = SXPAR(hdr, hdr_model[i].card, count = ct)
        if ct eq 0 then Value = hdr_model[i].Default
                
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
        
        junk = where(strmatch(['CHECKSUM','DATASUM','GAIN','RDNOISE', 'END'], card, /fold_case), ct)
        if ct > 0 then continue
        
        if strmatch(card,'EXPID*',/fold_case) then continue

        loopset= 0
        foreach cc, ['SIGBSHXX','CENBSHXX','AVGBSHXX','STDBSHXX',$
                     'SIGASHXX','CENASHXX','AVGASHXX','STDASHXX'] do begin
            if strmatch(card, STRMID(cc,0,strlen(cc)-3)+'*',/fold_case) then begin
                loopset = 1
                break
            endif
        endforeach
        if keyword_set(loopset) then continue
        
        if strmatch(card,'SN2*', /fold_case) then begin
            comment = strsplit(hdr[i],'/',/extract)
            if n_elements(comment) ge 2 then comment = comment[1] else comment = ''
            sxaddpar, mergedHdr, card, sxpar(hdr, card), comment, before='END'
            continue
        endif
 
        if strmatch(card,'SNC*', /fold_case) then begin
            comment = strsplit(hdr[i],'/',/extract)
            if n_elements(comment) ge 2 then comment = comment[1] else comment = ''
            sxaddpar, mergedHdr, card, sxpar(hdr, card), comment, before='END'
            continue
        endif
        
        if keyword_set(dropcards) then begin
            junk = where(strmatch(dropcards,card,/fold_case), ct)
        endif
        if ct eq 0 then begin
            splog, card+' missing from datamodel - Not propagated to outputs'
            print, sxpar(hdr, card)
        endif
    endfor
    hdr = mergedHdr
    
    
end
