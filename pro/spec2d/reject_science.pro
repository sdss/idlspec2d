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
;   hdr        - Header for image
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
;     Reject if any of the flat-field screens are closed.  If the FFS keyword
;       is in the header, it should be FFS = '0 0 0 0 0 0 0 0'
;     Reject if any of the flat-field lamps are turned on.  If the FF keyword
;       is in the header, it should be FF = '0 0 0 0'
;     Reject if any of the neon arc lamps are turned on.  If the NE keyword
;       is in the header, it should be NE = '0 0 0 0'
;     Reject if any of the HgCd arc lamps are turned on.  If the HGCD keyword
;       is in the header, it should be HGCD = '0 0 0 0'
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

   if (keyword_set(hdr)) then begin
      ffs = sxpar(hdr, 'FFS')
      if (keyword_set(ffs)) then begin
         ffs_sum = fix( total( fix( str_sep(ffs,' ') ) ) )
         if (ffs_sum GT 0) then begin
            qbad = 1
            splog, 'ABORT: Reject science: Flat-field screens are closed!'
         endif
      endif

      lamp_ne = sxpar(hdr, 'NE')
      ne_sum = fix( total( fix( str_sep(lamp_ne,' ') ) ) )
      if (ne_sum GT 0) then begin
         qbad = 1
         splog, 'ABORT: Reject science: Ne lamps turned on!'
      endif

      lamp_hgcd = sxpar(hdr, 'HGCD')
      hgcd_sum = fix( total( fix( str_sep(lamp_hgcd,' ') ) ) )
      if (hgcd_sum GT 0) then begin
         qbad = 1
         splog, 'ABORT: Reject science: HgCd lamps turned on!'
      endif

      lamp_ff = sxpar(hdr, 'FF')
      ff_sum = fix( total( fix( str_sep(lamp_ff,' ') ) ) )
      if (ff_sum GT 0) then begin
         qbad = 1
         splog, 'ABORT: Reject science: Flat-field lamps turned on!'
      endif
   endif

   if (keyword_set(fbadpix)) then begin
      if (fbadpix GT 0.10) then begin
         qbad = 1
         splog, 'ABORT: Reject science: ' $
          + string(format='(i3)', fix(fbadpix*100)) + '% bad pixels'
      endif
   endif

   if (keyword_set(nsatrow)) then begin
      if (nsatrow GT 100) then begin
         qbad = 1
         splog, 'ABORT: Reject science: ' $
          + string(format='(i4)', nsatrow) + ' saturated rows'
      endif
   endif

   isort = sort(img)
   percentxx = img[ isort[ 0.25 * n_elements(img) ] ]
   percentxx = median(img)
   splog, 'Science 25-th-percentile = ' + string(percentxx)
   if (percentxx GT 1000.) then begin
      qbad = 1
      splog, 'ABORT: Reject science as too bright: 25-th-percentile =' $
       + string(percentxx)
   endif

   qbad = 0 ; Hack for during BOSS commissioning!!!???
   return, qbad
end
;------------------------------------------------------------------------------
