;+
; NAME:
;   reject_arc
;
; PURPOSE:
;   Decide whether an arc is bad.
;
; CALLING SEQUENCE:
;   qbad = reject_arc(img, [ hdr, nsatrow=, fbadpix= ] )
;
; INPUTS:
;   img        - Raw arc image
;
; OPTIONAL INPUTS:
;   hdr        - Header for image
;   nsatrow    - Returned from SDSSPROC()
;   fbadpix    - Returned from SDSSPROC()
;
; OUTPUTS:
;   qbad       - Return 1 if a arc is bad, 0 if good.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Decide if this arc is bad:
;     Reject if more than 2% of the pixels are marked as bad.
;     Reject if more than 100 rows are saturated.
;     Reject if the flat-field screens are not closed, which should appear
;       as FFS = '1 1 1 1 1 1 1 1' in the header.
;     Reject if the arc lamps are not turned on, which should appear
;       as NE = '1 1 1 1' and HGCD = '1 1 1 1' in the header.
;       Accept if either set of lamps is properly turned on.
;
;   Note that the FFS, FF, NE and HGCD keywords only appear at MJD >= 51629.
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
function reject_arc, img, hdr, nsatrow=nsatrow, fbadpix=fbadpix

   qbad = 0

   if (keyword_set(hdr)) then begin
      ffs = sxpar(hdr, 'FFS')
      if (keyword_set(ffs)) then begin
         ffs_sum = fix( total( fix( str_sep(ffs,' ') ) ) )
         if (ffs_sum LT 8) then begin
            qbad = 1
            splog, 'WARNING: Reject arc: Flat-field screens not closed!'
         endif
      endif

      lamp_ne = sxpar(hdr, 'NE')
      lamp_hgcd = sxpar(hdr, 'HGCD')
      lamp_hear = sxpar(hdr, 'HEAR')

      if strmatch(string(sxpar(hdr,'CARTID')), '*FPS-S*', /fold_case) then obs='LCO' else obs='APO'
      
      if keyword_set(lamp_ne) then begin
          ne_sum = fix( total( fix( str_sep(lamp_ne,' ') ) ) )
          ne_max = n_elements(str_sep(lamp_ne,' '))
      endif else begin
          ne_sum = 0
      endelse
      
      if keyword_set(lamp_hgcd) then begin
          hgcd_sum = fix( total( fix( str_sep(lamp_hgcd,' ') ) ) )$
          hgcd_max = n_elements(str_sep(lamp_hgcd,' '))
      endif else begin
          hgcd_sum = 0
      endelse
      
      if keyword_set(lamp_hear) then begin
          hear_sum = fix( total( fix( str_sep(lamp_hear,' ') ) ) )$
          hear_max = n_elements(str_sep(lamp_hear,' '))
      endif else begin
          hear_sum = 0
      endelse
    
      if strmatch(obs, 'APO', /fold_case) then begin
         if (ne_sum LT 4) then $
          splog, 'WARNING: ' + strtrim(ne_sum,2) +   '/'+strtrim(ne_max)+     ' Ne lamps are off'
         if (hgcd_sum LT 4) then $
          splog, 'WARNING: ' + strtrim(hgcd_sum,2) + '/'+strtrim(hgcd_max)+ ' HgCd lamps are off'

         if (ne_sum LT 4 AND hgcd_sum LT 4) then begin
            qbad = 1
            splog, 'WARNING: Reject arc: Neither Ne nor HgCd lamps are off!'
         endif
      endif else begin
         if (ne_sum LT 4) then $
            splog, 'WARNING: ' + strtrim(ne_sum,2) +   '/'+strtrim(ne_max)+     ' Ne lamps are off'
         if (hear_sum LT 4) then $
            splog, 'WARNING: ' + strtrim(hear_sum,2) + '/'+strtrim(hear_max)+ ' HeAr lamps are off'

         if (ne_sum LT 4 AND hear_sum LT 4) then begin
            qbad = 1
            splog, 'WARNING: Reject arc: Neither Ne nor HeAr lamps turned on!'
         endif
      endelse
      
   endif

   if (keyword_set(fbadpix)) then begin
      if (fbadpix GT 0.02) then begin
         qbad = 1
         splog, 'WARNING: Reject arc: ' $
          + string(format='(i3)', fix(fbadpix*100)) + '% bad pixels'
      endif
   endif

   if (keyword_set(nsatrow)) then begin
      if (nsatrow GT 100) then begin
         qbad = 1
         splog, 'WARNING: Reject arc: ' $
          + string(format='(i4)', nsatrow) + ' saturated rows'
      endif
   endif

   return, qbad
end
;------------------------------------------------------------------------------
