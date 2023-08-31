function delhdrsec, hdr, sec
    ;if strmatch('REDUCTION',sec, /fold_case) then print, hdr
    start = where(strmatch(hdr, '*'+sec+'*', /fold_case), ct)
    if ct eq 0 then return, hdr
    if ct gt 1 then begin
        start = where(strmatch(hdr, '*'+STRUPCASE(sec)+'*'), ct)
        if ct eq 0 then return, hdr
    endif
    rawhdr = hdr
    foreach card, rawhdr[start+1:*] do begin
        if strmatch(card, ' *', /fold_case) then break
        card = strtrim((strsplit(card,'=',/extract))[0],2)
        ;splog, sec+':', card
        sxdelpar, hdr, card
    endforeach
    remove,[start, start-1], hdr
    return, hdr
end
