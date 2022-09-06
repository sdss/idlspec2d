
pro get_field_type, fieldid=fieldid, mjd=mjd, legacy=legacy, plates=plates, fps=fps
    plates=0
    legacy=0
    fps=0
    if keyword_set(fieldid) then begin
        if long(fieldid) lt 15000 then begin
            legacy = 1
        endif else begin
            if long(fieldid) lt 16000 then begin
                plates = 1
            endif else fps=1
        endelse
    endif else begin
        if keyword_set(mjd) then begin
            if long(mjd) lt 59030 then begin
                legacy = 1
            endif else begin
                if long(mjd) lt 59550 then begin
                    plates = 1
                endif else fps=1
            endelse
        endif
    endelse
end
