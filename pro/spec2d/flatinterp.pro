; This was called by EXTRACT_OBJECT, but isn't used now.
function flatinterp, fflat, minflat, nsmooth=nsmooth

     if NOT keyword_set(nsmooth) then nsmooth=15
     smoothfflat = fflat
     nfiber = (size(fflat))[2]

     for i=0,nfiber - 1 do begin
       bad = where(fflat[*,i] LE minflat, nbad)
       good = where(fflat[*,i] GT minflat, ngood)
       if nbad GT 0 AND ngood GT 10 then begin
         interp = fflat[*,i]
         interp[good] = smooth(fflat[good,i],nsmooth<ngood,/edge_truncate)
         mask = [0,fflat[*,i] LE minflat,0]
         i2 = [interp[good[0]], interp, interp[good[ngood-1]]] 
         smoothfflat[bad,i] = (djs_maskinterp(i2,mask))[bad+1]
       endif
     endfor

return, smoothfflat
end
        
     
       
       
     
   

