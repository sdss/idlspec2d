;+
; NAME:
;   sshiftvec
;
; PURPOSE:
;   Shift vector or image (line-at-a-time) using a damped sinc function.
;
; CALLING SEQUENCE:
;   simage = sshiftvec( fimage, shift, [ sincrad=sincrad, dampfac=dampfac, $
;    eps=eps ] )
;
; INPUTS:
;   fimage     - Input vector (1D) or image (2D)
;   shift      - Distance to shift
;
; OPTIONAL KEYWORDS:
;   sincrad    - Half-width of sinc convolution kernal; default to 10 pixels
;   dampfac    - Damping factor for gaussian; default to 3.25
;   eps        - Smallest fractional shift allowed; default to 1.0e-5
;
; OUTPUTS:
;   simage     - Shifted image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine is based upon the IDL routine SSHIFT from Marc Buie,
;   which in turn is based upon the Zodiac routine shiftc/sshift.
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   Dynamic link to sshiftvec.c
;
; REVISION HISTORY:
;   14-Apr-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function sshiftvec, fimage, shift, sincrad=sincrad, dampfac=dampfac, eps=eps

   ; Need 2 parameters
   if (N_params() LT 2) then begin
      print, 'Syntax - simage = sshiftvec( fimage, shift, [ sincrad=sincrad, $'
      print, ' dampfac=dampfac, eps=eps ] )'
      return, -1
   endif

   if (NOT keyword_set(sincrad)) then sincrad = 10L
   if (NOT keyword_set(dampfac)) then dampfac = 3.25
   if (NOT keyword_set(eps)) then eps = 1.0e-5

   szf = size(fimage)
   szs = size(shift)
   if (szf[0] EQ 1) then begin
      nx = szf[1]
      ny = 1L
      if (szs[0] EQ 0) then shiftvec = shift $
       else message, 'Invalid dimensions for SHIFT'
   endif else if (szf[0] EQ 2) then begin
      nx = szf[1]
      ny = szf[2]
      if (szs[0] EQ 0) then shiftvec = replicate(shift,ny) $
       else if (szs[0] EQ 1 AND szs[1] EQ ny) then shiftvec = shift $
       else message, 'Invalid dimensions for SHIFT'
   endif else begin
      message, 'Invalid dimensions for FIMAGE'
   endelse

   ibad = where(shiftvec GT nx, nbad)
   if (nbad GT 0) then message, 'Shift too large'

   simage = float(0 * fimage)
   soname = filepath('libspec2d.'+idlutils_so_ext(), $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='lib')
   result = call_external(soname, 'sshiftvec', $
    nx, ny, float(fimage), float(shiftvec), simage, $
    long(sincrad), float(dampfac), float(eps))

   return, simage
end
;------------------------------------------------------------------------------
