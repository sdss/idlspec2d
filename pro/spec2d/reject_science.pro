;+
; NAME:
;   reject_science
;
; PURPOSE:
;   Decide whether a science exposure is bad.
;
; CALLING SEQUENCE:
;   qbad = reject_science(img, [ hdr, nsatrow=, fbadpix= ] )
;
; INPUTS:
;   img        - Raw science image
;
; OPTIONAL INPUTS:
;   hdr        - Header for image (unused)
;   nsatrow    - Returned from SDSSPROC()
;   fbadpix    - Returned from SDSSPROC()
;
; OUTPUTS:
;   qbad       - Return 1 if a science exposure is bad, 0 if good.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Decide if this science exposure is bad:
;     Reject if more than 10% of the pixels are marked as bad.
;     Reject if more than 100 rows are saturated.
;     Reject if the 25-th percentile is more than 1000 electrons.
;       This percentile should be very low (of order 10 electrons), since
;       it will be the counts between fibers.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   15-Jun-2001  Written by D. Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function reject_science, img, hdr, nsatrow=nsatrow, fbadpix=fbadpix

   qbad = 0

   if (keyword_set(fbadpix)) then begin
      if (fbadpix GT 0.10) then begin
         qbad = 1
         splog, 'WARNING: Reject science: ' $
          + string(format='(i3)', fix(fbadpix*100)) + '% bad pixels'
      endif
   endif

   if (keyword_set(nsatrow)) then begin
      if (nsatrow GT 100) then begin
         qbad = 1
         splog, 'WARNING: Reject science: ' $
          + string(format='(i4)', nsatrow) + ' saturated rows'
      endif
   endif

   isort = sort(img)
   percent80 = img[ isort[ 0.25 * n_elements(img) ] ]
   if (percent80 GT 1000.) then begin
      qbad = 1
      splog, 'WARNING: Reject science as too bright: 25-th-percentile =' $
       + string(percent80)
   endif

   return, qbad
end
;------------------------------------------------------------------------------
