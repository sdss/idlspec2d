;+
; NAME:
;   pixelmask_bits
;
; PURPOSE:
;   Return an integer with a single bit set which matches bitlabel
;
; CALLING SEQUENCE:
;   pixelmask_bits(bitlabel)
;
; INPUTS:
;   bitlabel   - String to match to corresponding bit
;
; OUTPUTS:
;   One integer with 1 bit set 
;   returns 0 if no label is matched
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   pixelmask[i] = pixelmask_bits('FULLREJECT') 
;
; BUGS:
;   ?  We may want to return unsigned long instead
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   23-Jan-2000 Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
function pixelmask_bits, bitlabel

     pixelbits = ['NEARBADPIXEL', $      ; Bad pixel within 3 pixels of trace
                  'LOWFLAT', $           ; Flat field less than 0.5
                  'FULLREJECT',  $       ; Pixel fully rejected in extraction
                  'PARTIALREJECT',  $    ; Some pixels rejected in extraction
                  'SCATTEREDLIGHT',   $  ; Scattered light significant
                  'CROSSTALK',$          ; Cross-talk significant
                  'SKYLEVEL']           ; Sky background is > 10*flux
                    

     ss = strpos(pixelbits,strupcase(bitlabel))
	
     match = where(ss NE -1, nmatch)
     
     if (nmatch NE 1) then return, 0b

     return, 2^match[0]
end

	
