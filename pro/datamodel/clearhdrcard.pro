
function clearhdrcard, hdr, card, value=value
    idx = where(strmatch(hdr, card+'*',/fold_case), ct)
    if ct eq 0 then return, hdr
    if ct gt 1 then begin
        foreach i, idx do begin
            if strmatch(strtrim((strsplit(hdr[i],'=',/extract))[0],2),card, /fold_case) then begin
                idx = i
                break
            endif
        endforeach
    endif
    comment = strsplit(hdr[idx], '/', /extract)
    if n_elements(comment) eq 2 then comment = comment[1] else comment = ''

    if strlen(strtrim(hdr[idx+1],2 )) eq 0 then begin ; strmatch(hdr[idx+1], ' *') then begin
        loc = strtrim((strsplit(hdr[idx-1],'=',/extract))[0],2)
        sxdelpar, hdr, card
        if keyword_set(value) then begin
            sxaddpar, hdr, card, value, comment, /NULL, after=loc
        endif else begin
            sxaddpar, hdr, card, !Values.F_NAN, comment, /NULL, after=loc
        endelse
    endif else begin
        if strlen(strtrim(hdr[idx-1],2 )) ne 0 then begin
        ; if not strmatch(hdr[idx-1], ' *') then begin
            loc = strtrim((strsplit(hdr[idx+1],'=',/extract))[0],2)
            sxdelpar, hdr, card
            if keyword_set(value) then begin
                sxaddpar, hdr, card, value, comment, /NULL, before=loc
            endif else begin
                sxaddpar, hdr, card, !Values.F_NAN, comment, /NULL, before=loc
            endelse
        endif else begin
            splog,card
        endelse
    endelse
    return, hdr
end
