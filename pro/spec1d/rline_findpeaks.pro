;+
; NAME:
;   rline_findpeaks
;
; PURPOSE:
;  Given a two dimensional array for Equivalent Width (EW) and
;  weights (EW inverse variance), this function returns a structure
;  containing peak positions, heights, and significance level (SN).
;  
;  Only postive peaks are detected, so transform the array is negative
;  peaks are sought.
;
; CALLING SEQUENCE:
;   peak_list = rline_findpeaks(ew, ewinv, npeak=, threshold=, minsep=)
;
; INPUTS:
;   ew          - Positive equivalent width array
;   ewinv       - Associated inverse variance
;
; OPTIONAL INPUTS:
;   npeak       - Maximum number of peaks to locate   (default 20)
;   threshold   - Minimum significance level (sigma) to return (default 5.0)
;   minsep      - Minimum pixel separation to between peaks   (default 2.5)
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   peak_list   - Array of structures with peak information 
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   find_npeaks
;
; DATA FILES:
;
; REVISION HISTORY:
;  05-Jan-2001  Written by S. Burles, Fermiland
;  17-Dec-2001  Commented, updated to work with new find_npeaks.pro
;-
;------------------------------------------------------------------------------
function rline_findpeaks, ew, ewinv, npeak = npeak, threshold=threshold, $
     minsep=minsep

    npix  = (size(ew,/dim))[0]
    nspec = (size(ew,/dim))[1]
    xtab = lindgen(npix)

    if NOT keyword_set(npeak) then npeak = 20
    if NOT keyword_set(threshold) then threshold=5.0
    if NOT keyword_set(minsep) then minsep=2.5
   
    sn = ew * sqrt(ewinv)

    ttemp = { fiber : -1L, x : -1.0, y: -1.0, sn : 0.0, xerr : 0.0 , $
              class : 'UNK'}

    full_list = 0


    for i=0, nspec -1 do begin

       xpeak = find_npeaks(sn[*,i], nfind=npeak, ypeak=ypeak, xerr=xerr, $
                           minsep=minsep)

       high = where(ypeak GE threshold, nhigh)
       if nhigh GT 1 then begin
         mask = lonarr(nhigh) + 1
         for j=nhigh-1,1,-1 do begin
           identical = where(abs(xpeak[high[j]] - xpeak[high[0:j-1]]) LT minsep)
           if identical[0] NE -1 then mask[j] = 0
         endfor
         goodmask = where(mask, nhigh)
         if nhigh GT 0 then high=high[goodmask]
       endif 

       if nhigh GT 0 then begin
         linterp, xtab, ew[*,i], xpeak[high], y
         tt = replicate(ttemp, nhigh)
         tt.fiber = i
         tt.x = xpeak[high]
         tt.y = y
         tt.sn = ypeak[high]
         tt.xerr = xerr[high]
          
         if keyword_set(full_list) then full_list = [full_list, tt] $
         else full_list = tt

       endif
    endfor

    return, full_list
end

    

    

