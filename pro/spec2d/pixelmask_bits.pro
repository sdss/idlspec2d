;+
; NAME:
;   pixelmask_bits
;
; PURPOSE:
;   Return mask value corresponding to a mask condition for either FIBERMASK
;   or PIXELMASK.
;
; CALLING SEQUENCE:
;   mask = pixelmask_bits(label)
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
;   This routine is also called by FIBERMASK_BITS().
;
; EXAMPLES:
;   mask = pixelmask_bits('FULLREJECT') 
;
; BUGS:
;
; PROCEDURES CALLED:
;   sdss_flagval()
;
; DATA FILES:
;
; REVISION HISTORY:
;   23-Jan-2000 Written by S. Burles, Chicago
;   27-Jan-2000 Changed from signed int to signed long.
;   14-Jul-2000 Combine with FIBERMASK_BITS(), and use data file in /etc
;               subdirectory (DJS).
;   02-Apr-2002 Now use the even more generalized call to SDSS_FLAGVAL().
;-
;------------------------------------------------------------------------------
function pixelmask_bits, label

   return, sdss_flagval('SPPIXMASK', label)
end
;------------------------------------------------------------------------------
