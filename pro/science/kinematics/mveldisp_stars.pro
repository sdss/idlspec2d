Pro mveldisp_stars, galflux, galsig, galwave, starflux, starsig, starwave, $
sample, czmin, czmax, imin, imax, fullresult

goodmain=where((sample.z_primtarget EQ 64 OR sample.z_primtarget EQ 96) AND sample.z LT 0.3)

if (n_elements(goodmain) LT imax) then imax=n_elements(goodmain)-1

imin1=0
imax1=n_elements(goodmain)-1

for i=imin1,imax1 do begin & $
print,'object',i, sample[goodmain[i]].plate, sample[goodmain[i]].mjd, $
sample[goodmain[i]].fiberid & $
mveldisp, galflux[*,goodmain[i]],galsig[*,goodmain[i]],galwave[*,goodmain[i]],starflux[*,*], starsig[*,*], starwave[*,*], result, redshifts=sample[goodmain[i]].z,czmin=czmin, czmax= czmax, klo_cut=0.016, khi_cut=0.23, maxsig=6,sigmastep=0.4,/noquot  & $

if (size(fullresult,/tname) EQ 'STRUCT') then begin & $
fullresult = [fullresult, result] & $
endif else begin & $
fullresult = result & $
endelse & $

endfor

namefull=string('veldisp_MAIN'+string(imin,format='(i4.4)')+$
                 '_'+string(imax,format='(i4.4)')+'.fit')
mwrfits, fullresult, namefull

return
end