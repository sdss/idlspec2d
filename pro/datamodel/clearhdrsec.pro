function clearhdrsec, hdr, sec, exclude = exclude
    start = where(strmatch(hdr, '*'+sec+'*', /fold_case), ct)
    if ct eq 0 then return, hdr
    if ct gt 1 then begin
        start = where(strmatch(hdr, '*'+STRUPCASE(sec)+'*'), ct)
        if ct eq 0 then return, hdr
    endif
    foreach card, hdr[start+1:*] do begin
        if strmatch(card, ' *', /fold_case) then return, hdr
        card = strtrim((strsplit(card,'=',/extract))[0],2)
        if keyword_set(exclude) then begin
            junk = where(strmatch(exclude, card,/fold_case),ct)
            if ct eq 1 then continue
        endif
        hdr = clearhdrcard(hdr, card)
    endforeach
    return, hdr
end
