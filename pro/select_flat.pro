;+
; NAME:
;   select_flat
;
; PURPOSE:
;   Select best flat from a flat-field structure produced with SPCALIB.
;
; CALLING SEQUENCE:
;   bestflat = select_flat( flatstruct, tai )
;
; INPUTS:
;   flatstruct - Structure array with extracted flat calibration information
;   tai        - Time of observation for object frame (from header card)
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   bestflat   - Structure array with extracted flat calibration information
;                for best flat selected for the time specified by TAI.
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
; REVISION HISTORY:
;   25-Jan-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

function select_flat, flatstruct, tai

   indx = where(flatstruct.qbad EQ 0, ct)
   if (ct EQ 0) then begin
      splog, 'No good flats'
      return, 0
   endif else if (ct EQ 1) then begin
      splog, 'One good flat: ' + flatstruct[indx[0]].name
      return, flatstruct[indx[0]]
   endif

   timediff = flatstruct[indx].tai - tai

   ; For now, selected the nearest good flat
   junk = min(abs(timediff), ii)
   ibest = indx[ii]

   return, flatstruct[ibest]
end
;------------------------------------------------------------------------------
