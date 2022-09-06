function yanny_par_fc, hdr, key, count=count, indx=indx
    val = yanny_par(hdr, key, count=count, indx=indx)
    if not keyword_set(val) then begin
        idx = where(strmatch(hdr, key+'*', /fold_case), ct)
        if ct ne 0 then begin
                foreach k, hdr[idx] do begin
                        k = (strsplit(k,' ',/extract))[0]
                        if strmatch(k, key, /fold_case) then begin
                                key = k
                                break
                        endif
                endforeach
        endif
    endif
    return, yanny_par(hdr, key, count=count, indx=indx)
end
