function struct_to_yannyhdr,extname, hdr_struct=hdr_struct
    indx=where(strmatch(hdr_struct.EXTNAME, extname, /fold_case), ct)
    this_hdr_struct = hdr_struct[indx]
    keys = tag_names(this_hdr_struct)
    nhead=n_elements(keys)
    yanny_hdr=strarr(nhead)
    
    for i=0, nhead-1 do begin
        if strlen(this_hdr_struct.(i)) eq 0 then continue
        key = strtrim(keys[i],2)
        if key eq 'EXTNAME' then extname = this_hdr_struct.(i)
        if key eq 'PLUGDIR' then PLUGDIR = this_hdr_struct.(i)
    endfor
    for i=0, nhead-1 do begin
        if strlen(this_hdr_struct.(i)) eq 0 then continue
        key = strtrim(keys[i],2)
        if key eq 'EXTNAME' then continue
        if key eq 'PLUGDIR' then continue
        ;key = restore_key(key, plugdir+extname)
        yanny_hdr[i] = key+' '+this_hdr_struct.(i)
    endfor
    return, yanny_hdr
end

