;------------------------------------------------------------------------------
;+
; NAME:
;   fiber_rollcall
;
; PURPOSE:
;   Print the "roll call" of how many plug-map bits are set.
;
; CALLING SEQUENCE:
;   fiber_rollcall, andmask, loglam
;
; INPUTS:
;   andmask    - Bit mask
;   loglam     - Wavelength mapping in log10(Ang)
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   splog
;   sdss_flagname()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   17-Feb-2004  Written by D. Schlegel & S. Burles; copied from PLATESN.
;-
;------------------------------------------------------------------------------
pro fiber_rollcall, andmask, loglam

   dims = size(andmask, /dimens)
   nfiber = dims[1]

   ;----------
   ; Print roll call of bad fibers and bad pixels.
   ; Assume that all mask bits under bit #16 are for the entire fiber.

   bitlabel = sdss_flagname('SPPIXMASK', 2UL^32-1, /silent)
   bitnum = where(bitlabel NE '', nlabel)
   bitlabel = bitlabel[bitnum]

   ncam = 4
   camnames = ['b1','r1','b2','r2']
   camwave1 = [3500, 6000, 3500, 6000]
   camwave2 = [6000, 9500, 6000, 9500]
   specidimg = bytarr(dims) + 1B
   specidimg[*,nfiber/2:nfiber-1] = 2B

   gwave = where(loglam GT alog10(4000) AND loglam LT alog10(5500))
   rwave = where(loglam GT alog10(5600) AND loglam LT alog10(6900))
   iwave = where(loglam GT alog10(6910) AND loglam LT alog10(8500))

   rollcall = fltarr(ncam)

   splog, ' '
   splog, camnames, format='(24x,4a7)'
   splog, format='(24x,4(" ------"))'
   for ilabel=0, nlabel-1 do begin
      for icam=0, ncam-1 do begin
         specid = fix( strmid(camnames[icam],1) )
         thiscam = strmid(camnames[icam],0,1)
         fib1 = (specid-1) * (nfiber/2)
         fib2 = fib1 + (nfiber/2) - 1

         if (bitnum[ilabel] LT 16) then begin
            ; CASE: Fiber mask
            ; Look at only the g-band or i-band wavelengths, so that
            ; we ignore any overlap in wavelength between cameras.
            if (thiscam EQ 'b') then indx = gwave $
             else indx = iwave
            qmask = (andmask[indx,fib1:fib2] AND 2L^bitnum[ilabel]) NE 0
            rollcall[icam] = total( total(qmask,1) NE 0 )
         endif else begin
            ; CASE: Pixel mask
            ; Include overlap in wavelength between cameras, calling
            ; any wavelengths < 6000 Ang the blue camera, and
            ; any wavelengths > 6000 Ang the red camera.
            indx = where(loglam GT alog10(camwave1[icam]) $
             AND loglam LT alog10(camwave2[icam]), nindx)
            qmask = (andmask[indx,fib1:fib2] AND 2L^bitnum[ilabel]) NE 0
            rollcall[icam] = 100.0 * total(qmask NE 0) $
             / (nindx * (fib2-fib1+1))
         endelse
      endfor
      if (bitnum[ilabel] LT 16) then begin
         splog, 'N(fiber) '+bitlabel[ilabel]+'               ', $
          long(rollcall), format='(a24,4i7)'
      endif else begin
         splog, '%(pixel) '+bitlabel[ilabel]+'               ', $
          rollcall, format='(a24,4f7.2)'
      endelse
   endfor

   return
end
;------------------------------------------------------------------------------
