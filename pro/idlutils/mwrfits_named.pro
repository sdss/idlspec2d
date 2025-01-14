pro mwrfits_named, data, filename, hdr=hdr, name=name, desc=desc, $
                 _EXTRA = _EXTRA
                 
    auto_keys = ['XTENSION','BITPIX','NAXIS','NAXIS1','NAXIS1','PCOUNT','GCOUNT','TFIELDS','SIMPLE','EXTEND']
    if keyword_set(hdr) then begin
        hdr_tmp = hdr[*]
        sxdelpar, hdr_tmp, auto_keys
    endif
    if keyword_set(name) then begin
        name = STRUPCASE(name)
        if keyword_set(desc) then $
            sxaddpar, hdr_tmp, 'EXTNAME', name, desc $
        else sxaddpar, hdr_tmp, 'EXTNAME', name

    endif
    mwrfits, data, filename, hdr_tmp, _EXTRA=_EXTRA
    hdr_tmp = 0
    return
end
