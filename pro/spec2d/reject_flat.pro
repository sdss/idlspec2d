;+
; NAME:
;   reject_flat
;
; PURPOSE:
;   Decide whether a flat is bad.
;
; CALLING SEQUENCE:
;   qbad = reject_flat(img, [ hdr, nsatrow=, fbadpix=, percent80thresh= ] )
;
; INPUTS:
;   img        - Raw flat-field image
;
; OPTIONAL INPUTS:
;   hdr        - Header for image
;   nsatrow    - Returned from SDSSPROC()
;   fbadpix    - Returned from SDSSPROC()
;                percent80thresh - threshold for the number of pixels
;                                  that can fall below 80% flux
;
; OUTPUTS:
;   qbad       - Return 1 if a flat is bad, 0 if good.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Decide if this flat is bad:
;     Reject if more than 2% of the pixels are marked as bad.
;     Reject if more than 100 rows are saturated.
;     Reject if the 80-th percentile is less than 1000 electrons.
;     Reject if the flat-field screens are not closed, which should appear
;       as FFS = '1 1 1 1 1 1 1 1' in the header.
;     Reject if the flat-field lamps are not turned on, which should appear
;       as FF = '1 1 1 1' in the header.
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
function reject_flat, img, hdr, nsatrow=nsatrow, fbadpix=fbadpix, $
        percent80thresh=percent80thresh, noreject=noreject

  if (NOT keyword_set(percent80thresh)) then percent80thresh=400.

   qbad = 0

   if (keyword_set(hdr)) then begin
      ffs = sxpar(hdr, 'FFS')
      if (keyword_set(ffs)) then begin
         ffs_sum = fix( total( fix( str_sep(ffs,' ') ) ) )
         if (ffs_sum LT 8) then begin
            qbad = 1
            splog, 'WARNING: Reject flat: Flat-field screens not closed!'
         endif
      endif

      ff = sxpar(hdr, 'FF')
      if (keyword_set(ff)) then begin
         ff_sum = fix( total( fix( str_sep(ff,' ') ) ) )
         if (ff_sum LT 4) then begin
            qbad = 1
            splog, 'WARNING: Reject flat: Flat-field lamps not turned on!'
         endif
      endif
      if not strmatch(sxpar(hdr,'HARTMANN'), 'out*', /fold_case) then begin
         splog,'WARNING: Hartmann doors closed'
      endif
      

      lamp_ne = sxpar(hdr, 'NE')
      lamp_hgcd = sxpar(hdr, 'HGCD')
      lamp_hear = sxpar(hdr, 'HEAR')

      if strmatch(string(sxpar(hdr,'CARTID')), '*FPS-S*', /fold_case) then obs='LCO' else obs='APO'
      
      if keyword_set(lamp_ne) then begin
          ne_sum = fix( total( fix( str_sep(lamp_ne,' ') ) ) )
          ne_max = n_elements(str_sep(strtrim(lamp_ne,2),' '))
      endif else begin
          ne_sum = 0
      endelse
      
      if keyword_set(lamp_hgcd) then begin
          hgcd_sum = fix( total( fix( str_sep(lamp_hgcd,' ') ) ) )
          hgcd_max = n_elements(str_sep(strtrim(lamp_hgcd,2),' '))
      endif else begin
          hgcd_sum = 0
      endelse
      
      if keyword_set(lamp_hear) then begin
          hear_sum = fix( total( fix( str_sep(lamp_hear,' ') ) ) )
          hear_max = n_elements(str_sep(strtrim(lamp_hear,2),' '))
      endif else begin
          hear_sum = 0
      endelse
      if strmatch(obs, 'APO', /fold_case) then begin
        if (ne_sum gt 0) then begin
            splog, 'WARNING: Reject Flat: ' + strtrim(ne_sum,2) +   '/'+strtrim(ne_max,2)+     ' Ne lamps are On'
            qbad = 1
        endif
        if (hgcd_sum gt 0) then begin
            splog, 'WARNING: Reject Flat: ' + strtrim(hgcd_sum,2) + '/'+strtrim(hgcd_max,2)+ ' HgCd lamps are On'
            qbad = 1
        endif
      endif else begin
        if (ne_sum gt 0) then begin
            splog, 'WARNING: Reject Flat: ' + strtrim(ne_sum,2) +   '/'+strtrim(ne_max,2)+     ' Ne lamps are On'
            qbad =1
        endif
        if (hear_sum gt 0) then begin
            splog, 'WARNING: Reject Flat: ' + strtrim(hear_sum,2) + '/'+strtrim(hear_max,2)+ ' HeAr lamps are On'
            qbad = 1
        endif
      endelse
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
   if (percent80 LT percent80thresh) then begin
       qbad = 1
       splog, 'WARNING: Reject flat as too faint: 80-th-percentile =' $
                    + string(percent80)
   endif else begin
        if (percent80 LT 1.3*percent80thresh) then begin
            splog, 'WARNING: Flat is borderline faint (but adequate): 80-th-percentile =' +string(percent80)
        endif else splog, 'flat 80-th-pecentile ='+string(percent80)
   endelse
   if keyword_set(noreject) then begin
       if keyword_set(qbad) then splog, 'WARNING: OVERRIDING BAD FLAT FLAG'
       qbad = 0
   endif
   return, qbad
end
;------------------------------------------------------------------------------
