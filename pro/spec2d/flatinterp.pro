function flatinterp, fflat, minflat, nsmooth=nsmooth

     if NOT keyword_set(nsmooth) then nsmooth=15
     smoothfflat = fflat
     nfiber = (size(fflat))[2]

     for i=0,nfiber - 1 do begin
       bad = where(fflat[*,i] LE minflat, nbad)
       good = where(fflat[*,i] GT minflat, ngood)
       if nbad GT 0 AND ngood GT 10 then begin
         interp = fflat[*,i]
         interp[good] = smooth(fflat[good,i],nsmooth)
         smoothfflat[bad,i] = (djs_maskinterp(interp,fflat[*,i] LE minflat))[bad]
       endif
     endfor

return, smoothfflat
end
        
     
       
       
     
   

