;+
; NAME:
;   fibermask_bits
;
; PURPOSE:
;   Return a byte with a single bit set which matches bitlabel
;
; CALLING SEQUENCE:
;   fibermask_bits(bitlabel)
;
; INPUTS:
;   bitlabel   - String to match to corresponding bit
;
; OUTPUTS:
;   One byte with 1 bit set 
;   returns 0 if no label is matched
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fibermask[i] = fibermask_bits('BADARC') 
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   23-Jan-2000 Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
function fibermask_bits, bitlabel

     fiberbits = ['NOPLUG', $            ; Not in plugmap file
                     'BADTRACE', $       ; Bad trace in trace320crude
                     'BADFLAT',  $       ; Low counts in fflat
                     'BADARC',   $       ; Bad arc solution
                     'MANYBADCOLUMNS',$  ; >10% bad columns 
                     'MANYREJECTED',  $  ; >10% rejected in extraction
                     'LARGESHIFT']       ; Large shift between flat and object?
                    

     ss = strpos(fiberbits,strupcase(bitlabel))
	
     match = where(ss NE -1, nmatch)
     
     if (nmatch NE 1) then return, 0b

     return, 2b^match[0]
end

	
