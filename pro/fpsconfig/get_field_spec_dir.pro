

function get_field_spec_dir, topdir, run2d, field, mjd, epoch=epoch, $
                             lite=lite, custom_name=custom_name
    if keyword_set(lite) then stype = 'lite' else stype = 'full'
    if keyword_set(epoch) then begin
        subdir = [run2d,'spectra','epoch',stype,get_field_grp(field),$
                 field_to_string(field)]
    endif else begin
        if keyword_set(custom_name) then begin
            subdir = [run2d,'spectra',get_field_grp(field), stype, $
                      get_field_grp(field), field_to_string(field)]
        endif else begin
            subdir = [run2d,'spectra','daily', stype,get_field_grp(field),$
field_to_string(field)]
        endelse
    endelse
    return, djs_filepath(strtrim(mjd,2), root_dir=topdir, subdirectory=subdir)
end

