
function qso_sortpeaks, full_list, tt, linelist=linelist
 
    if NOT keyword_set(linelist) then linelist = 'mgII.lines'

    readcol, linelist, linelabel, restwave, f_value, damping, $
              format='a,d,f,f'

    label = linelabel+string(format='(i4)',restwave)
    logwave = alog10(restwave)
    pixel_sep = logwave*1.0e4

    mg1 = (where(fix(restwave) EQ 2796))[0]
    if mg1 EQ -1 then return, 0

    pixel_sep = pixel_sep - pixel_sep[mg1]
    nsep = n_elements(pixel_sep)

    nspec = n_elements(tt)
    all_mgii = 0

    for i=0, nspec -1 do begin
    
      any = where(full_list.fiber EQ i, nthis)

      if nthis GE 2 then begin

        results = look_mgii(full_list[any], tt[i], i,  pixel_sep, $
          restwave, f_value, label)

        if keyword_set(results) then begin
          if NOT keyword_set(all_mgii) then all_mgii = results $
          else all_mgii = struct_append(all_mgii, results)
        endif
      endif
    endfor

   return, all_mgii
end
