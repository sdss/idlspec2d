;+
; NAME:
;   plate_rotate
;
; PURPOSE:
;   Rotate (RA, dec) positions for one plate center to another plate center
;
; CALLING SEQUENCE:
;   plate_rotate, racen1, deccen1, racen2, deccen2, ra1, dec1, ra2, dec2
;
; INPUTS:
;   racen1, deccen1  - Center of pointing 1
;   racen2, deccen2  - Center of pointing 2
;   ra1, dec1        - positions on plate 1
;
; OUTPUTS:
;   ra2, dec2        - positions on plate 2
;
; COMMENTS:
;   This routine determines the rotation that maps (racen1, deccen1)
;    onto (racen2, deccen2), and applies it to (ra1, dec1) (which may
;    be vectors.  By convention, North is "up" on both plates 1 & 2. 
;
;   This routine is useful for the bright stars survey, where one must
;    drill holes for several pointings in the same plate. 
;
; BUGS:
;
; REVISION HISTORY:
;   2001-Sep-07  Written by D. Finkbeiner, Princeton.
;
;-
;------------------------------------------------------------------------------
pro plate_rotate, racen1, deccen1, racen2, deccen2, ra1, dec1, ra2, dec2

  if (n_params() EQ 0) then begin 
     print, 'plate_rotate, racen1, deccen1, racen2, deccen2, ra1, dec1, ra2, dec2'
     return
  endif

; convert input ra1, dec1 to unit vectors
  uv = ll2uv([[ra1-racen1], [dec1]], /double)  ; array [N, 3]

  dtor = !dpi/180.d
  sd = sin((deccen2-deccen1)*dtor)
  cd = cos((deccen2-deccen1)*dtor)
  
  x = uv[*, 0]*cd - uv[*, 2]*sd
  y = uv[*, 1]
  z = uv[*, 0]*sd + uv[*, 2]*cd
  uv2 = [[x], [y], [z]]
      
  radec = uv2ll(uv2) ; array [N, 2]
      
  ra2 = radec[*, 0]+racen2
  dec2 = radec[*, 1]

; wrap RA
  ra2 = (ra2+720.d) MOD 360

  return
end

pro test

  racen1 = 45
  deccen1 = 90.

  racen2 = 135
  deccen2 = 90.
  
  ra1  = randomn(iseed, 4)+racen1
  dec1 = randomn(iseed, 4)+deccen1
  plot,ra1,dec1,ps=1,/yno
  oplot, [racen1], [deccen1], ps=6, syms=3

  plate_rotate, racen1, deccen1, racen2, deccen2, ra1, dec1, ra2, dec2

  plot,ra2,dec2,ps=1,/yno
  oplot, [racen2], [deccen2], ps=6, syms=3

  return
end

