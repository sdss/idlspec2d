;
; Simple procedure to parse header and return best estimates of time-stamps
;
; Check header for keywords detailing exposure time-stamps
;
; If TAI-END and TAI-BEG exist, define TAI-MID = (TAI-BEG + TAI-END)/2
; if not, define TAI-END = TAI - 60.
;                TAI-BEG = TAI-END - EXPTIME
;                TAI-MID = (TAI-BEG + TAI-END)/2
;
; S. Burles, 28 Jan 2002

pro get_tai, hdr, tai_beg, tai_mid, tai_end

   exptime = sxpar(hdr, 'EXPTIME')
   tai     = sxpar(hdr, 'TAI')
   tai_beg = sxpar(hdr, 'TAI-BEG')
   tai_end = sxpar(hdr, 'TAI-END')

   if (tai_beg GT 4.0d9 AND tai_end GT 4.0d9 AND $
    tai_end - tai_beg GE 0) then begin
      tai_mid = (tai_end + tai_beg)/2.0
      exptime_guess = tai_end - tai_beg
      ; The following condition is possible if the exposure was paused for
      ; more than 2 minutes.
      if (exptime_guess GT exptime + 120) then $
       splog, "Warning: (TAI-END) - (TAI-BEG) > EXPTIME + 120 sec"
      ; The following should never happen.
      if (exptime_guess LT exptime + 20) then $
       splog, "Warning: (TAI-END) - (TAI-BEG) < EXPTIME + 20 sec"
   endif else if (tai GT 4.0d9) then begin
      tai_end = tai - 60.0   ; average buffer for readout 
      tai_beg = tai_end - exptime
      tai_mid = (tai_beg + tai_end)/2.0
   endif else begin
      tai_mid = 0
      splog, 'WARNING: TAI,TAI-BEG,TAI-END all appear to be invalid!'
   endelse

   return
end
