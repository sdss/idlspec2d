;+
; NAME:
;   reject_flat
;
; PURPOSE:
;   Decide whether a flat is bad.
;
; CALLING SEQUENCE:
;   qbad = reject_flat(img, [ hdr, nsatrow=, fbadpix= ] )
;
; INPUTS:
;   img        - Raw flat-field image
;
; OPTIONAL INPUTS:
;   hdr        - Header for image
;   nsatrow    - Returned from SDSSPROC()
;   fbadpix    - Returned from SDSSPROC()
;
; OUTPUTS:
;   qbad       - Return 1 if a flat is bad, 0 if good.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
; Decide if this flat is bad:
;   Reject if more than 2% of the pixels are marked as bad.
;   Reject if more than 100 rows are saturated.
;   Reject if the 80-th percentile is less than 1000 electrons.
;   Reject if the flat-field screens are not closed, which should appear
;     as FFS = '1 1 1 1 1 1 1 1' in the header.
;   Reject if the flat-field lamps are not turned on, which should appear
;     as FF = '1 1 1 1' in the header.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   25-Jan-2001  Written by D. Schlegel, Princeton.
;                This code is copied out of SPCALIB.
;-
;------------------------------------------------------------------------------
function reject_flat, img, hdr, nsatrow=nsatrow, fbadpix=fbadpix

   qbad = 0

   if (keyword_set(hdr)) then begin
      ffs = sxpar(hdr, 'FFS')
      if (keyword_set(ffs)) then begin
         if (strtrim(ffs) NE '1 1 1 1 1 1 1 1') then begin
            qbad = 1
            splog, 'WARNING: Reject flat: Flat-field screens not closed!'
         endif
      endif

      ff = sxpar(hdr, 'FF')
      if (keyword_set(ff)) then begin
         if (strtrim(ff,2) NE '1 1 1 1') then begin
            qbad = 1
            splog, 'WARNING: Reject flat: Flat-field lamps not turned on!'
         endif
      endif
   endif

   if (keyword_set(fbadpix)) then begin
      if (fbadpix GT 0.02) then begin
         qbad = 1
         splog, 'WARNING: Reject flat: ' $
          + string(format='(i3)', fix(fbadpix*100)) + '% bad pixels'
      endif
   endif

   if (keyword_set(nsatrow)) then begin
      if (nsatrow GT 100) then begin
         qbad = 1
         splog, 'WARNING: Reject flat: ' $
          + string(format='(i4)', nsatrow) + ' saturated rows'
      endif
   endif

   isort = sort(img)
   percent80 = img[ isort[ 0.80 * n_elements(img) ] ]
   if (percent80 LT 1000.) then begin
      qbad = 1
      splog, 'WARNING: Reject flat as too faint: 80-th-percentile =' $
       + string(percent80)
   endif

   return, qbad
end
;------------------------------------------------------------------------------
