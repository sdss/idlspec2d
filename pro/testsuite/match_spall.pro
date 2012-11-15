;+
;
; NAME:
;  match_spall
;
; PURPOSE:
;  Match two spall structure rows on common
;  plate-mjd-fiberid, making the assumption
;  that they are ganged sensibly by Nfib.
;
; USAGE:
;  match_spall, sp1, sp2, m1, m2 [, nfib=nfib]
;
; INPUTS:
;  sp1, sp2: two spAll structures with the tags
;    PLATE, MJD, and FIBERID at a minimum
;  nfib (optional): number of fibers per plate,
;    default to 1000.
;
; OUTPUTS:
;  m1, m2: matching indices such that sp1[m1] is
;  a row-by-row match to sp2[m2].
;
; WRITTEN:
;  bolton@utah 2012nov
;
;-

pro match_spall, sp1, sp2, m1, m2, nfib=nfib

if (not keyword_set(nfib)) then nfib = 1000L

ns1 = n_elements(sp1)
ns2 = n_elements(sp2)

nfib = 1000L
np1 = ns1 / nfib
np2 = ns2 / nfib

idx1 = nfib * lindgen(np1)
idx2 = nfib * lindgen(np2)

pl1 = sp1[idx1].plate
mj1 = sp1[idx1].mjd
key1 = string(pl1, format='(i4.4)') + '-' + string(mj1, format='(i5.5)')

pl2 = sp2[idx2].plate
mj2 = sp2[idx2].mjd
key2 = string(pl2, format='(i4.4)') + '-' + string(mj2, format='(i5.5)')

key_uniq = [key1, key2]
key_uniq = key_uniq[sort(key_uniq)]
key_uniq = key_uniq[uniq(key_uniq)]

n_uniq = n_elements(key_uniq)

; Find those that are in both arrays:
inboth = replicate(0B, n_uniq)
for i = 0L, n_uniq - 1 do inboth[i] = (total(key1 eq key_uniq[i]) gt 0) and (total(key2 eq key_uniq[i]) gt 0) 

key_both = key_uniq[where(inboth)]
n_common = n_elements(key_both)

m1 = replicate(-1L, nfib * n_common)
m2 = replicate(-1L, nfib * n_common)

for i = 0L, n_common - 1 do begin & $
   wh1 = (where(key1 eq key_both[i]))[0] & $
   m1[nfib*i:nfib*(i+1)-1] = nfib * wh1 + lindgen(nfib) & $
   wh2 = (where(key2 eq key_both[i]))[0] & $
   m2[nfib*i:nfib*(i+1)-1] = nfib * wh2 + lindgen(nfib) & $
endfor

return
end
