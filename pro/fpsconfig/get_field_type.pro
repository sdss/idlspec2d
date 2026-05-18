
pro get_field_type, fieldid=fieldid, mjd=mjd, legacy=legacy, plates=plates, fps=fps
    COMMON generations_data_block, GENERATIONS_DATA

    IF ~ISA(JSON_DATA) THEN BEGIN
        GENERATIONS_DATA = JSON_PARSE(filepath('generations.json',root_dir=getenv('IDLSPEC2D_DIR'), $
                                                          subdir=['python','boss_drp','etc']),/TOSTRUCT)
    ENDIF

    plates=0
    legacy=0
    fps=0
    if keyword_set(fieldid) then begin
        if long(fieldid) lt GENERATIONS_DATA.legacy[0].field_range[1] then begin
            legacy = 1
        endif else begin
            if long(fieldid) lt GENERATIONS_DATA.plates[0].field_range[1] then begin
                plates = 1
            endif else fps=1
        endelse
    endif else begin
        ; TODO: Handle obs if difference are needed post SDSSV

        if keyword_set(mjd) then begin
            if long(mjd) lt GENERATIONS_DATA.legacy[0].mjd_range[0].apo[1] then begin
                legacy = 1
            endif else begin
                if long(mjd) lt GENERATIONS_DATA.legacy[0].mjd_range[0].apo[1] then begin
                    plates = 1
                endif else fps=1
            endelse
        endif
    endelse
end





; Parse JSON into IDL structure
data = JSON_PARSE(json_text)