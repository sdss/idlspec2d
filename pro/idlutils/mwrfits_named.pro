pro mwrfits_named, data, filename, hdr=hdr, name=name, desc=desc, $
                 _EXTRA = _EXTRA
                 
    auto_keys = ['XTENSION','BITPIX','NAXIS','NAXIS1','NAXIS1','PCOUNT','GCOUNT','TFIELDS','SIMPLE','EXTEND']
    if keyword_set(hdr) then begin
        hdr_tmp = hdr[*]
        sxdelpar, hdr_tmp, auto_keys

        keep =[]
        foreach card, hdr_tmp, idx do begin
            if strmatch(card, 'COMMENT*', /FOLD_CASE) eq 1 then begin
                if strmatch(strtrim(card, 2),'COMMENT', /fold_case) eq 1 then continue
                if strmatch(card, 'Comment*\*\*\* Column names \*\*\**', /fold_case)  eq 1 then continue
                if strmatch(card, 'Comment*\*\*\* Column formats \*\*\**', /fold_case)  eq 1 then continue
                if strmatch(card, 'Comment*\*\*\* End of mandatory fields \*\*\**', /fold_case) eq 1 then continue
            endif
            keep = [keep, idx]
        endforeach
        if n_elements(keep) gt 0 then hdr_tmp = hdr_tmp[keep] else hdr_tmp = 0
    endif
    if keyword_set(name) then begin
        name = STRUPCASE(name)
        if keyword_set(desc) then $
            sxaddpar, hdr_tmp, 'EXTNAME', name, desc $
        else sxaddpar, hdr_tmp, 'EXTNAME', name

    endif
    
    data = fix_empty_strings(data)
    mwrfits, data, filename, hdr_tmp, _EXTRA=_EXTRA
    hdr_tmp = 0
    return
end
