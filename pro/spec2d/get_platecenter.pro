;+
; NAME:
;   get_platecenter
;
; PURPOSE:
;   Return our best guess of the observed plate center.
;
; CALLING SEQUENCE:
;   get_platecenter, objhdr, plugmap, ra, dec
;
; INPUTS:
;   objhdr     - Header cards from object image
;   plugsort   - Plugmap structure w/just one element per fiber.
;
; OUTPUTS:
;   ra, dec    - the center position of the observed fibers. Crude.
;
; NOTES: 
;   See idlspec2d/7201 for some background. Basically, we need to
;   get the observed plate position from the fiber positions: neither
;   the exposure header nor the plugmap keywords always hold the right
;   values. [I wonder how the telescope is pointed...]
;
;   There can be several pointings per plate, and the distribution of
;   observed fibers over a plate can be pretty strange. There are
;   definitely more correct answers, but I'm starting with the middle
;   of the observed ra & dec ranges.
;
; BUGS: Galactic cap plates can get a totally wrong RA.
;
pro get_platecenter, plugsort, ra, dec
    ras = plugsort.ra
    if (max(plugsort.ra) - min(plugsort.ra) GT 180) then $
       ras[where(ras LT 180)] += 360
    ra = (min(ras) + max(ras)) / 2
    if ra GE 360 then ra -= 360
    
    dec = (min(plugsort.dec) + max(plugsort.dec)) / 2
 end

