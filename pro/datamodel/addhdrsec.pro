
function addhdrsec, hdr, sec, key, value, comment, before = before
    temphdr = [""]
    sxaddpar, temphdr, key, value, comment
    if keyword_set(before) then begin
        start = where(strmatch(hdr, '*'+before+'*', /fold_case), ct)
        hdr = [hdr[0:start-2], '', '        '+ STRUPCASE(sec), temphdr[0], hdr[start-1:-1]]
    endif else begin
        hdr = [hdr[0:-2], '', '        '+ STRUPCASE(sec),temphdr[0], hdr[-1]]
    endelse

    return, hdr
end
