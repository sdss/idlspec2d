;+
; NAME:
;   fibermask_bits
;
; PURPOSE:
;   Return mask value corresponding to a mask condition for FIBERMASK.
;
; CALLING SEQUENCE:
;   mask = fibermask_bits(label)
;
; INPUTS:
;   label      - String name specifying bit
;
; OUTPUTS:
;   mask       - Signed long set to 2^BIT, with BIT specified by LABEL.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mask = fibermask_bits('NOPLUG') 
;
; BUGS:
;
; PROCEDURES CALLED:
;   pixelmask_bits()
;
; REVISION HISTORY:
;   23-Jan-2000 Written by S. Burles, Chicago
;   14-Jul-2000 Combine with PIXELMASK_BITS() (DJS).
;-
;------------------------------------------------------------------------------
function fibermask_bits, label

   return, pixelmask_bits(label)
end
;------------------------------------------------------------------------------
