;------------------------------------------------------------------------------
FUNCTION apo_color2hex, colorname
    CASE STRUPCASE(STRTRIM(colorname,2)) OF
        'RED': hexname = '#FF0000'
        'YELLOW': hexname = '#FFFF00'
        ELSE: hexname = 'black'
    ENDCASE
    RETURN, hexname
END
;------------------------------------------------------------------------------
