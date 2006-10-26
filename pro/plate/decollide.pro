;+
; NAME:
;   random_decollide
; PURPOSE:
;   Take a set of RAs and DECs and find a random, decollided set
; CALLING SEQUENCE:
;   random_decollide, 
; INPUTS:
; REVISION HISTORY:
;   26-Oct-2006  Written by MRB, NYU
;-
;------------------------------------------------------------------------------
function random_decollide, ra, dec, holesize=holesize, seed=seed

iran=shuffle_indx(n_elements(ra), seed=seed)
ing=spheregroup(ra[iran], dec[iran], holesize, firstg=firstg, $
                nextg=nextg, multg=multg)
ng=max(ing)+1L

holeok=bytarr(n_elements(ra))
for i=0L, ng-1L do begin
    if(multg[i] eq 1) then begin
        holeok[firstg[i]]=1
    endif else begin
        indx=lonarr(multg[i])-1L
        ok=bytarr(multg[i])
        k=0L
        j=firstg[i]
        while(j ne -1L) do begin
            indx[k]=j
            j=nextg[j]
        endwhile
        ok[0]=1L
        for k=1L, multg[i]-1L do begin
            iok=where(ok)
            spherematch, ra[indx[iok]], dec[indx[iok]], $
              ra[indx[k]], dec[indx[k]], holesize, m1, m2, d12
            if(m2[0] eq -1) then ok[k]=1
        endfor
        iok=where(ok)
        holeok[indx[iok]]=1
    endelse
endfor

return, idec

end
