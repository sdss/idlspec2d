function field_to_string, field

  if n_elements(field) gt 1 then begin
    sfield = strarr(n_elements(field))
    for i=0, n_elements(field)-1 do sfield[i] = field_to_string(field[i])
    return, sfield
  endif
  if long(field) lt 0 then field = long(0) 
  return, strtrim(string(field,f='(i6.6)'),2) 


end
