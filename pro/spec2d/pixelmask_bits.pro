;+
; NAME:
;   pixelmask_bits
;
; PURPOSE:
;   Return mask value corresponding to a mask condition for PIXELMASK.
;
; CALLING SEQUENCE:
;   mask = pixelmask_bits(bitlabel)
;
; INPUTS:
;   bitlabel   - String name specifying bit
;
; OUTPUTS:
;   mask       - Signed long set to 2^BIT, where BIT is the bit specified by
;                BITLABEL, or set to 0 if no label is matched
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mask = pixelmask_bits('FULLREJECT') 
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   23-Jan-2000 Written by S. Burles, Chicago
;   27-Jan-2000 Changed from signed int to signed long
;-
;------------------------------------------------------------------------------
function pixelmask_bits, bitlabel

     pixelbits = ['NEARBADPIXEL', $      ; Bad pixel within 3 pixels of trace
                  'LOWFLAT', $           ; Flat field less than 0.5
                  'FULLREJECT',  $       ; Pixel fully rejected in extraction
                  'PARTIALREJECT',  $    ; Some pixels rejected in extraction
                  'SCATTEREDLIGHT',   $  ; Scattered light significant
                  'CROSSTALK',$          ; Cross-talk significant
                  'NOSKY',$              ; No sky subtraction
                  'SKYLEVEL',$           ; Sky background is > 10*flux
                  'NODATA',$             ; No data available in combine B-spline
                  'COMBINEREJ']          ; Rejected in combine B-spline

     ss = strpos(pixelbits,strupcase(bitlabel))

     match = where(ss NE -1, nmatch)

     if (nmatch NE 1) then return, 0

     return, 2L^match[0]
end

