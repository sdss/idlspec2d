
function get_field_dir, topdir, run2d, field, custom=custom
    if keyword_set(custom) then begin
        field = strtrim(field,2)
        fieldgrp = get_field_grp(field,/custom)
            
        dir_ = djs_filepath(field, root_dir=topdir,$
                            subdirectory=[run2d, 'fields',fieldgrp])
    endif else begin
        if strcmp(field,'*') then begin
            zfield = field_to_string(0)
            zfieldgrp = repstr(get_field_grp(0),'0','?')
            zfield = repstr(zfield,'0','?')
            dir_ = djs_filepath(zfield, root_dir=topdir, $
                                 subdirectory=[run2d,'fields',zfieldgrp])
        endif else begin
            dir_ = djs_filepath(field_to_string(field), root_dir=topdir, $
                                 subdirectory=[run2d,'fields',get_field_grp(field)])
        endelse
    endelse
    return, dir_
end
