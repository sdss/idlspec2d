;+
; NAME:
;   redmonster
;
; PURPOSE:
;   Search for the Red Monster (contiguous region of bad sky residuals),
;   and set a pixelmask bit.
;
; CALLING SEQUENCE:
;   redmonster, relloglam, relchi2, [ objloglam, filtsz=, thresh=, pixelmask= ]
;
; INPUTS:
;   relloglam  - Log10(Angstroms) for RELCHI2 vector
;   relchi2    - Relative chi^2 vector from sky fiber residuals
;
; OPTIONAL INPUTS:
;   objlogam   - Log10(Angstroms) image for object frame [NPIX,NFIBER];
;                required input for modifying PIXELMASK
;   filtsz     - Filter size for looking for Red Monster; default to 25 pix,
;                but not more than the number of pixels in RELLOGLAM
;   thresh     - Treshold in relative chi^2 for identifying REDMONSTER;
;                in the pixel mask; default to 4.0.
;   pixelmask  - If this and OBJLOGLAM are specified, then add the REDMONSTER
;                bit to this mask [NPIX,NFIBER]
;
; OUTPUTS:
;   pixelmask  - (Modified.)
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   pixelmask_bits()
;   splog
;
; REVISION HISTORY:
;   14-Mar-2001  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro redmonster, relloglam, relchi2, objloglam, filtsz=filtsz1, $
 thresh=thresh, pixelmask=pixelmask

   if (keyword_set(filtsz1)) then filtsz = filtsz1 $
    else filtsz = 25
   if (NOT keyword_set(thresh)) then thresh = 4.0

   nbin = n_elements(relloglam)
   filtsz = filtsz < (nbin - 1L)
   filtwd = (filtsz-1)/2 ; Half-width of filter

   ;----------
   ; Filter the relative chi^2 vector to be the lowest quartile value
   ; within a moving filter of size FILTSZ.

   filtchi2 = fltarr(nbin)
   for ibin=0, nbin-1 do begin
      i0 = (ibin - filtwd) > 0L
      i1 = (ibin + filtwd) < (nbin - 1)
      isort = sort(relchi2[i0:i1])
      filtchi2[ibin] = (relchi2[i0:i1])[isort[(i1-i0+1)/4]]
   endfor

   ;----------
   ; Find all (filtered) chi^2 values above a threshold, and grow
   ; this by convolving with a boxcar of size FILTSZ.

   badmask = filtchi2 GE thresh
   badmask = convol(long(badmask), lonarr(filtsz)+1L, /center, /edge_truncate) $
    NE 0

   ;----------
   ; Identify each contiguous group of bad pixels, then trigger a warning
   ; message and add a bit to the pixel mask.

   while (total(badmask) NE 0) do begin
      i0 = (where(badmask EQ 1))[0]
      i1 = (where(badmask[i0:nbin-1] EQ 0))[0]
      if (i1 EQ -1) then i1 = nbin - 1 $
       else i1 = i1 + i0
      peakchi = sqrt(max(relchi2[i0:i1]))
      splog, string(fix(10^relloglam[i0]), fix(10^relloglam[i1]), peakchi, $
       format='("WARNING: Bad sky residuals at ", i5, " to ", i5, " Ang; Peak chi=", f6.2)')
      if ( 10^relloglam[i0] GT 6250 AND 10^relloglam[i0] LT 6500 $
       AND 10^relloglam[i1] GT 6500 AND 10^relloglam[i1] LT 6750) then begin
         splog, string(fix(10^relloglam[i0]), fix(10^relloglam[i1]), peakchi, $
          format='("WARNING: Red Monster at ", i5, " to ", i5, " Ang; Peak chi=", f6.2)')
      endif
      badmask[i0:i1] = 0

      if (keyword_set(objloglam) AND keyword_set(pixelmask)) then begin
         indx = where(objloglam GE relloglam[i0] AND objloglam LE relloglam[i1])
         if (indx[0] NE -1) then $
          pixelmask[indx] = pixelmask[indx] OR pixelmask_bits('REDMONSTER')
      endif
   endwhile

   return
end
;------------------------------------------------------------------------------
