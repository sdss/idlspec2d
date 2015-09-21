function plate_to_string, plate

if n_elements(plate) gt 1 then begin
    splate = strarr(n_elements(plate))
    for i=0, n_elements(plate)-1 do splate[i] = plate_to_string(plate[i])
    return, splate
endif

if plate lt 10000 then $
	return, strtrim(string(plate,f='(i4.4)'),2) $
else $
	return, strtrim(string(plate,f='(i6.6)'),2)

end

