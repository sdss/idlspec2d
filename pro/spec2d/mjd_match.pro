;+
; NAME:
;   mjd_match
;
; PURPOSE:
;   Decide if any of these MJD's are within the bounds specified by
;   MJD,MJSTART,MJEND.
;
; CALLING SEQUENCE:
;   qmjd = mjd_match(mjdlist, [mjd=, mjstart=, mjend=])
;
; INPUTS:
;   mjdlist    - One MJD or an array of MJD's.
;
; OPTIONAL INPUTS:
;   mjd        - Match against these specific MJD's.
;   mjstart    - Starting MJD.
;   mjend      - Ending MJD.
;
; OUTPUT:
;   qmjd       - Set to 1 if any of the MJD's in MJDLIST are within
;                the bounds specified by MJD,MJSTART,MJEND.
;                Otherwise, set to 0.
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   22-Aug-2001  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function mjd_match, mjdlist, mjd=mjd, mjstart=mjstart, mjend=mjend

   qmjd = 1B
   ; If MJD keyword set, then QMJD=1 only if any of the
   ; elements of MJDLIST1 agree with those in MJD.
   if (keyword_set(mjd)) then $
    qmjd = qmjd AND ((where( $
     mjdlist # (lonarr(n_elements(mjd))+1) - $
     (lonarr(n_elements(mjdlist))+1) # mjd EQ 0))[0] NE -1)

   ; If MJSTART keyword set, then QMJD=1 only if any of the
   ; elements of MJDLIST1 is >= MJSTART.
   if (keyword_set(mjstart)) then $
    qmjd = qmjd AND ((where(mjdlist GE mjstart))[0] NE -1)

   ; If MJEND keyword set, then QMJD=1 only if any of the
   ; elements of MJDLIST1 is <= MJEND.
   if (keyword_set(mjend)) then $
    qmjd = qmjd AND ((where(mjdlist LE mjend))[0] NE -1)

   return, qmjd
end
;------------------------------------------------------------------------------
