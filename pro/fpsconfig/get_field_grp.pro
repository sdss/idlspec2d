function get_field_grp, field, custom=custom
    if keyword_set(custom) then begin
        if n_elements(strsplit(field,'_',/extract)) eq 1 then begin
            return, field
        endif else begin
            return, strjoin((strsplit(field,'_',/extract))[0:-2],'_')
        endelse
    endif else begin
        if valid_num(string(field)) then begin
            return, strtrim(string(long(field)/1000,f='(i3.3)'),2)+'XXX'
        endif else begin
            return, field
        endelse
    endelse

end
