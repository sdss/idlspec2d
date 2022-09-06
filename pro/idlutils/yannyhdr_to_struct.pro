
function yannyhdr_to_struct, hdr, extname, plugdir=plugdir, hdr_struct=hdr_struct, cartid=cartid

    if not keyword_set(plugdir) then plugdir = file_dirname(extname)
    extname = file_basename(extname)
    nhead=n_elements(hdr)
    this_hdr_struct = create_struct('EXTNAME', extname, 'PLUGDIR', plugdir)
    for i=0, nhead-1 do begin
        if strmid(hdr[i],0,1) EQ "#" then continue
        if strmid(hdr[i],0,8) EQ "EVILSCAN" then continue
        if strlen(hdr[i]) eq 0 then continue
        key = (str_sep( strtrim(hdr[i],2), ' '))[0]
        if strlen(key) eq 0 then continue
        ct = 0
        key_flag = 0
        check = where(strmatch(tag_names(this_hdr_struct), key, /fold_case), ct)
        while ct ne 0 do begin
            key_flag = 1
            check = where(strmatch(tag_names(this_hdr_struct), key+'_'+strtrim(key_flag,2), /fold_case), ct)
        endwhile
        if keyword_set(key_flag) then key = key+'_'+strtrim(key_flag,2)
        this_hdr_struct = struct_addtags(this_hdr_struct, create_struct(key, strjoin(yanny_par_fc(hdr, key),' ')))
    endfor

    if not keyword_set(cartid) then begin
        if strtrim(yanny_par_fc(hdr,'observatory'),2) eq 'APO' then begin
            this_hdr_struct = struct_addtags(this_hdr_struct, create_struct('cartridgeId', 'FPS-N'))
        endif else this_hdr_struct = struct_addtags(this_hdr_struct, create_struct('cartridgeId', 'FPS-S'))
    endif else begin
        check = where(strmatch(tag_names(this_hdr_struct), 'cartridgeId', /fold_case), ct)
        if ct eq 0 then this_hdr_struct = struct_addtags(this_hdr_struct, create_struct('cartridgeId', cartid))
    endelse

    if keyword_set(hdr_struct) then begin

        hdr_struct_tags = tag_names(hdr_struct)
        this_hdr_struct_tags = tag_names(this_hdr_struct)

        hdr_struct_add_tags = []
        this_hdr_struct_add_tags =[]
        foreach tag, hdr_struct_tags do begin
            junk = where(this_hdr_struct_tags eq tag, ct)
            if ct eq 0 then this_hdr_struct_add_tags=[this_hdr_struct_add_tags,tag]
        endforeach
        foreach tag, this_hdr_struct_tags do begin
            if (tag eq 'EXTNAME') or (tag eq 'PLUGDIR') then continue
            junk = where(hdr_struct_tags eq tag, ct)
            if ct eq 0 then hdr_struct_add_tags=[hdr_struct_add_tags,tag]
        endforeach

        empty_hdr_struct = hdr_struct[0]
        struct_assign, {junk:!values.f_nan}, empty_hdr_struct
        empty_this_struct = replicate(this_hdr_struct[0], n_elements(hdr_struct))
        struct_assign, {junk:!values.f_nan}, empty_this_struct

        if n_elements(hdr_struct_add_tags) gt 0 then $
            hdr_struct = struct_addtags(hdr_struct,struct_selecttags(empty_this_struct, select_tags=hdr_struct_add_tags))
        if n_elements(this_hdr_struct_add_tags) gt 0 then $
            this_hdr_struct = struct_addtags(this_hdr_struct,struct_selecttags(empty_hdr_struct, select_tags=this_hdr_struct_add_tags))
  
        this_hdr_struct = struct_selecttags(this_hdr_struct, select_tags=tag_names(hdr_struct))
        hdr_struct=struct_append(hdr_struct, this_hdr_struct)

    endif else hdr_struct = this_hdr_struct
    return, hdr_struct
end
