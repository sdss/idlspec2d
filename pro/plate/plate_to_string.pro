function plate_to_string, plate

if plate lt 10000 then $
	return, strtrim(string(plate,f='(i4.4)'),2) $
else $
	return, strtrim(string(plate,f='(i6.6)'),2)

end

