FUNCTION fix_empty_strings, data, PLACEHOLDER=placeholder, VERBOSE=verbose
    COMPILE_OPT idl2

    ; Default placeholder is a single space
    IF N_ELEMENTS(placeholder) EQ 0 THEN placeholder = ' '

    t = TYPENAME(data)

    CASE t OF

        ;--------------------------------------------------------
        ; STRUCT or array of STRUCTS → recurse into each field
        ;--------------------------------------------------------
        'STRUCT': BEGIN
            n = N_ELEMENTS(data)
            tags = TAG_NAMES(data)
            FOR ie = 0, n-1 DO BEGIN
                FOR i = 0, N_ELEMENTS(tags)-1 DO BEGIN
                    newval = fix_empty_strings(data[ie].(i), $
                                               PLACEHOLDER=placeholder, $
                                               VERBOSE=verbose)
                    data[ie].(i) = newval
                ENDFOR
            ENDFOR
        END

        ;--------------------------------------------------------
        ; STRING or string array → replace all-zero-length entries
        ;--------------------------------------------------------
        'STRING': BEGIN
            lens = STRLEN(data)
            IF TOTAL(lens) EQ 0 THEN BEGIN
                data[*] = placeholder
                IF KEYWORD_SET(verbose) THEN $
                    PRINT, 'Replaced empty string array with placeholder "' + placeholder + '"'
            ENDIF ELSE BEGIN
                ; handle partial empties too
                empty_idx = WHERE(lens EQ 0, count)
                IF count GT 0 THEN BEGIN
                    data[empty_idx] = placeholder
                    IF KEYWORD_SET(verbose) THEN $
                        PRINT, 'Replaced ', count, ' empty string(s) with "' + placeholder + '"'
                ENDIF
            ENDELSE
        END

        ;--------------------------------------------------------
        ; All other data types → no change
        ;--------------------------------------------------------
        ELSE: BEGIN
            ; no modification
        END
    ENDCASE

    RETURN, data
END

