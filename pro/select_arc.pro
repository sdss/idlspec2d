;+
; NAME:
;   select_arc
;
; PURPOSE:
;   Select best arc from an arc structure produced with SPCALIB.
;
; CALLING SEQUENCE:
;   bestarc = select_arc ( arcstruct )
;
; INPUTS:
;   arcstruct  - Structure array with extracted arc calibration information
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   bestarc    - Structure with extracted arc calibration information
;                for that arc selected from ARCSTRUCT
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
;   djsig()
;
; REVISION HISTORY:
;   25-Jan-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

function select_arc, arcstruct

   indx = where(arcstruct.qbad EQ 0, ct)
   if (ct EQ 0) then begin
      splog, 'No good arcs'
      return, 0
   endif else if (ct EQ 1) then begin
      splog, 'One good arc: ' + arcstruct[indx[0]].name
      return, arcstruct[indx[0]]
   endif

   ; Select the arcs with the most matched arc lines; trim INDX to that list
   nmax = max( arcstruct[indx].nmatch )
   indx = indx[ where(arcstruct[indx].nmatch EQ nmax) ]

   ; Of that list, find the one with the arc with the smallest residuals
   ibest = 0
   bestsig = 0
   for ii=0, N_elements(indx)-1 do begin
      sig = djsig( *arcstruct[indx[ii]].xdif_tset )
      if (ii EQ 0 OR sig LT bestsig) then begin
         ibest = indx[ii]
         bestsig = sig
      endif
   endfor

   return, arcstruct[ibest]
end
;------------------------------------------------------------------------------
