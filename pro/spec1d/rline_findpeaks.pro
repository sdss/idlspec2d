function rline_findpeaks, ew, ewinv, npeak = npeak, threshold=threshold

    npix  = (size(ew,/dim))[0]
    nspec = (size(ew,/dim))[1]
    xtab = lindgen(npix)

    if NOT keyword_set(npeak) then npeak = 20
    if NOT keyword_set(threshold) then threshold=5.0
   
    sn = ew * sqrt(ewinv)

    ttemp = { fiber : -1L, x : -1.0, y: -1.0, sn : 0.0, xerr : 0.0 , $
              class : 'UNK'}

    full_list = 0


    for i=0, nspec -1 do begin

       xpeak = find_npeaks(sn[*,i], npeak=npeak, ypeak=ypeak, xerr=xerr)

       high = where(ypeak GE threshold, nhigh)

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

    

    

