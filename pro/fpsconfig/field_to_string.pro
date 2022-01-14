function field_to_string, field

  if n_elements(field) gt 1 then begin
    sconfig = strarr(n_elements(field))
    for i=0, n_elements(field)-1 do sfield[i] = field_to_string(field[i])
    return, sfield
  endif
    return, strtrim(string(field,f='(i6.6)'),2) 
  ;if field lt 10000 then $
  ;  return, strtrim(string(field,f='(i6.6)'),2) $
  ;else $
  ;  return, strtrim(field,2)

end
