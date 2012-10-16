;+
; NAME:
;  sdss_specobjid
; PURPOSE:
;  convert spectrum identifiers to CAS-style SPECID
; CALLING SEQUENCE:
;  specobjid= sdss_specobjid(plate, fiber, mjd, rerun [, line=, index=])
; INPUTS:
;  plate - [N] plate id
;  fiberid - [N] fiber number
;  mjd - [N] mjd
;  rerun - [N] reruns of object(s)
; OPTIONAL INPUTS:
;  line - a line measure, with this code number (1-indexed)
;  index - an index measure, with this code number (1-indexed)
; OUTPUTS:
;  specobjid - [N] resulting CAS-style specobjid(s)
; COMMENTS:
;  Cannot specify both "line" and "index" as non-zero
;  Bits are assigned in specobjid (a ulong64 #)
;    50-63 plate id number (14 bits)
;    38-49 fiber id number  (12 bits)
;    24-37 MJD (date) plate was observed minus 50000 (14 bits)
;    10-23 rerun (14 bits)
;    0-9   line/redshift/index: 0 for SpecObj, else number of
;          spectroscopic line (SpecLine) or index (SpecLineIndex)
;          (10 bits)
;  This definition differs from that in DR7.
;  MJD must be > 50000
;  Rerun can be an integer, like 26, or a string of the form 'vN_M_P',
;    where N, M and P are integers, with the restriction 5<=N<=6,
;    0<=M<=99, and 0<=P<=99. This is understood to be the RUN2D value
;    for a spectrum. In the latter case, the 14 bits for the rerun are
;    filled with (N-1)*10000+M*100+P
; BUGS:
;  No provision is made for v7 or greater software versions, or for
;   cases where the RUN1D value is relevant (not either '' or equal to
;   RUN2D).
; REVISION HISTORY:
;  Written by MRB, NYU 24-10-2008
; VERSION:
;  $Id: sdss_specobjid.pro 130412 2012-02-16 21:22:13Z blanton $
;------------------------------------------------------------------------------
FUNCTION sdss_specobjid, plate, fiber, mjd, in_rerun, line=line, index=index
    nplate=N_ELEMENTS(plate)
    IF nplate EQ 0 THEN MESSAGE, 'Must specify PLATE.'
    IF (KEYWORD_SET(line) GT 0 AND KEYWORD_SET(index) GT 0) THEN $
        MESSAGE, 'Cannot set both LINE and INDEX!'
    IF ~KEYWORD_SET(line) THEN line= REPLICATE(0L, nplate)
    IF ~KEYWORD_SET(index) THEN index= REPLICATE(0L, nplate)
    IF (N_ELEMENTS(fiber) NE nplate OR $
        N_ELEMENTS(mjd) NE nplate OR $
        N_ELEMENTS(in_rerun) NE nplate) THEN $
        MESSAGE, 'PLATE, FIBER, MJD, RERUN all need same size'
    mjdm5= mjd-50000L
    codenumber=0L
    ;
    ; Determine whether rerun is an integer
    ;
    IF SIZE(in_rerun,/TYPE) EQ 7 THEN BEGIN
        IF STREGEX(STRTRIM(in_rerun,2),'^[0-9]+$',/BOOLEAN) THEN BEGIN
            ;
            ; OK, even though it's a string, it is still a pure integer
            ;
            rerun = LONG(in_rerun)
        ENDIF ELSE BEGIN
            ;
            ; Here's where we have to deal with mapping rerun strings
            ; to integers.
            ;
           words= STREGEX(STRTRIM(in_rerun,2),'^v([0-9]+)_([0-9]+)_([0-9]+)', /SUB, $
                          /EXTRACT)
           If(N_ELEMENTS(words) ne 4) THEN $
              MESSAGE, 'Rerun name should be of form vN_M_P'
           rerun= (long(words[1])-5L)*10000L+ $
                  (long(words[2])*100L)+ $
                  (long(words[3]))
        ENDELSE
    ENDIF ELSE BEGIN
       rerun= in_rerun
    ENDELSE
    IF KEYWORD_SET(line) THEN codenumber= line
    IF KEYWORD_SET(index) THEN codenumber= index
    ibad= WHERE( plate LT 0 OR plate GE 2L^14 OR $
        fiber LT 0 OR fiber GT 2L^12 OR $
        mjdm5 LT 0 OR mjdm5 GE 2L^14 OR $
        rerun LT 0 OR rerun GE 2L^14 OR $
        codenumber LT 0 OR codenumber GE 2L^10, nbad)
    IF (nbad GT 0) THEN MESSAGE, 'Inputs out of bounds!'
    ;
    ; convert to ulong64 and clip
    ;
    ul2=ULONG64(2)
    ulplate= ULONG64(plate) AND (ul2^14-1)
    ulfiber= ULONG64(fiber) AND (ul2^12-1)
    ulmjdm5= ULONG64(mjdm5) AND (ul2^14-1)
    ulrerun= ULONG64(rerun) AND (ul2^14-1)
    ulcodenumber= ULONG64(codenumber) AND (ul2^10-1)
    ;
    ; load bits into right parts of objid
    ;
    specobjid= ULON64ARR(nplate)
    specobjid= specobjid OR $
        ISHFT(ulplate,50) OR $
        ISHFT(ulfiber,38) OR $
        ISHFT(ulmjdm5,24) OR $
        ISHFT(ulrerun,10) OR $
        ISHFT(ulcodenumber,0)
    RETURN, specobjid
END
