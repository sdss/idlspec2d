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
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/spMaskbits.par
;
; REVISION HISTORY:
;   23-Jan-2000 Written by S. Burles, Chicago
;   27-Jan-2000 Changed from signed int to signed long.
;   14-Jul-2000 Combine with FIBERMASK_BITS(), and use data file in /etc
;               subdirectory (DJS).
;-
;------------------------------------------------------------------------------
function pixelmask_bits, label

   ; Declare a common block so that the mask names are remembered between calls.
   common com_maskbits, maskbits

   if (NOT keyword_set(maskbits)) then begin
      maskfile = filepath('spMaskbits.par', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
      if (NOT keyword_set(maskfile)) then $
       message, 'File with mask bits not found'
      yanny_read, maskfile, pdat
      maskbits = *pdat[0]
      yanny_free, pdat
   endif

   imatch = where(strupcase(label) EQ strupcase(maskbits.label), ct)
   if (ct NE 1) then $
    message, 'Error matching bit name: ' + label

   return, 2L^(maskbits[imatch[0]].bit)
end
;------------------------------------------------------------------------------
