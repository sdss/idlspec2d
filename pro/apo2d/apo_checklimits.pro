;+
; NAME:
;   apo_checklimits()
;
; PURPOSE:
;   Convert output FITS file from APOREDUCE to HTML format.
;
; CALLING SEQUENCE:
;   markstring = apo_checklimits(flavor, field, camera, value, [ /html ] )
;
; INPUTS:
;   flavor     - FLAVOR to match in the opLimits file.
;   field      - FIELD to match in the opLimits file.
;   camera     - CAMERA to match in the opLimits file.
;   value      - Value to test in the opLimits file.  If this is a
;                string value, then one matches to STRVAL in that file.
;                Otherwise, test for values within [LOVALUE,HIVALUE].
;
; OPTIONAL INPUTS:
;   html       - If set, then convert the color name in MARKSTRING into
;                an HTML string.  For example, a return value of 'red'
;                becomes '<B><FONT COLOR="#FF0000">'.
;
; OUTPUT:
;   markstring - Return the COLOR from the opLimits file.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/opLimits.par
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
FUNCTION apo_checklimits, flavor, field, camera, value, html=html
    COMMON apo_limits, numlimits, textlimits
    markstring = ''
    IF ~KEYWORD_SET(value) THEN RETURN, markstring
    ;
    ; Read this Yanny file only the first time this routine is called,
    ; then save the limits in a common block.
    ;
    IF ~KEYWORD_SET(numlimits) THEN BEGIN
        limitfile = FILEPATH('opLimits.par', ROOT_DIR=GETENV('IDLSPEC2D_DIR'), $
            SUBDIRECTORY='examples')
        yanny_read, limitfile, pdata, stnames=stnames
        numlimits = *pdata[(WHERE(stnames EQ 'SPECLIMIT'))[0]]
        textlimits = *pdata[(WHERE(stnames EQ 'TEXTLIMIT'))[0]]
        yanny_free, pdata
    ENDIF
    IF (size(value,/tname) EQ 'STRING') THEN BEGIN
        ;
        ; Case of text limit
        ;
        FOR ilim=0, N_ELEMENTS(textlimits)-1 DO BEGIN
            IF (STRMATCH(field, textlimits[ilim].field) $
                && STRMATCH(camera, textlimits[ilim].camera) $
                && STRMATCH(flavor, textlimits[ilim].flavor) $
                && STRMATCH(STRTRIM(value,2), textlimits[ilim].strval)) THEN BEGIN
                markstring = textlimits[ilim].color
                IF KEYWORD_SET(html) THEN markstring = '<span style="color:' $
                    + apo_color2hex(markstring) + ';font-weight:bold;">'
            ENDIF
        ENDFOR
    ENDIF ELSE BEGIN
        ;
        ; Case of floating-point limit
        ;
        FOR ilim=0, N_ELEMENTS(numlimits)-1 DO BEGIN
            IF (STRMATCH(field, numlimits[ilim].field) $
                && strmatch(camera, numlimits[ilim].camera) $
                && strmatch(flavor, numlimits[ilim].flavor) $
                && value GE numlimits[ilim].lovalue $
                && value LE numlimits[ilim].hivalue) THEN BEGIN
                markstring = numlimits[ilim].color
                IF KEYWORD_SET(html) THEN markstring = '<span style="color:' $
                    + apo_color2hex(markstring) + ';font-weight:bold;">'
            ENDIF
        ENDFOR
    ENDELSE
    RETURN, markstring
END
;------------------------------------------------------------------------------
