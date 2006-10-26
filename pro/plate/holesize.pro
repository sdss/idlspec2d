;+
; NAME:
;   holesize
; PURPOSE:
;   Returns size in deg of a hole in an SDSS plate
; CALLING SEQUENCE:
;   sz= holesize([/guide])
; OPTIONAL KEYWORDS:
;   guide    - Return value for a guide hole
; COMMENTS:
;   Returns 55./3600. for a normal hole, 175/3600. for a guide
; REVISION HISTORY:
;   26-Oct-2006  MRB, NYU
;-
;------------------------------------------------------------------------------
function holesize, guide=guide

holesize=55./3600.
if(keyword_set(guide)) then holesize=175./3600.

return, holesize

end

