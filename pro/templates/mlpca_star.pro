;+
;
; NAME: mlpca_star
;
; PURPOSE:
;  Generate ML-PCA templates for stars
;
; bolton@utah 2010may
;
;-


pro mlpca_star, infile

infile = 'mlpcainput_star.txt'

readcol, infile, plate, mjd, fiberid, z, subc, format='L,L,L,F,A'

subc = strtrim(subc)

; Exclude CVs:
wh = where(subc ne 'CV')
plate = plate[wh]
mjd = mjd[wh]
fiberid = fiberid[wh]
z = z[wh]
subc = subc[wh]

nobj = n_elements(plate)

readspec, plate, fiberid, mjd=mjd, flux=flux, invvar=invvar, $
 loglam=loglam, zans=zans, /align, /silent

; Get everything aligned in rest-frame wavelength:
align_rest, in_loglam=loglam, in_flux=flux, in_invvar=invvar, z=z, $
 out_loglam=rest_loglam, out_flux=rest_flux, out_invvar=rest_ivar

; Require most spectra to be present for good pixels.
; This should be done in a more assiduous subtype-by-subtype manner,
; and thus there are conceivable cases where this might blow up in
; our face down the line.

ngoodspec = total(rest_ivar gt 0., 2)
mingoodspec = round (0.8 * nobj)
p_lo = min(where(ngoodspec ge mingoodspec))
p_hi = max(where(ngoodspec ge mingoodspec))
rest_flux = rest_flux[p_lo:p_hi,*]
rest_ivar = rest_ivar[p_lo:p_hi,*]
rest_loglam = rest_loglam[p_lo:p_hi]

; Number of unique subtypes:
u_subc = subc[sort(subc)]
u_subc = u_subc[uniq(u_subc)]

ntype = n_elements(u_subc)
npix = n_elements(rest_loglam)

startempset = dblarr(npix, ntype)
model = 0.d0 * rest_flux

; Loop over sub-types:
for i = 0L, ntype-1 do begin & $
   wh = where(subc eq u_subc[i]) & $
   startempset[*,i] = bolton_mlpca(objflux=rest_flux[*,wh], objivar=rest_ivar[*,wh], coeffs=coeffs, nbasis=1) & $
   model[*,wh] = startempset[*,i] # coeffs & $
endfor

; Inspect:
;i = -1L
;i++ & splot, rest_loglam, rest_flux[*,i] * (rest_ivar[*,i] gt 0) & soplot, rest_loglam, model[*,i], color=3
i=-1L
i++ & splot, 10.^rest_loglam, startempset[*,i] & print, u_subc[i]

   ;----------
   ; Construct header for output file

get_juldate, jd
mjdstr = STRING(LONG(jd-2400000L), FORMAT='(I5)')
outfile = 'spMLpcaStar-' + mjdstr + '.fits'
;outfile = 'spEigenStar-' + mjdstr + '.fits'
sxaddpar, hdr, 'OBJECT', 'STAR'
sxaddpar, hdr, 'COEFF0', rest_loglam[0]
sxaddpar, hdr, 'COEFF1', rest_loglam[1] - rest_loglam[0]
; Add a space to the name below, so that 'F' appears as a string and
; not as a logical.
for i=0, ntype-1 do $
 sxaddpar, hdr, 'NAME'+strtrim(string(i),2), u_subc[i]+' '

   ;----------
   ; Write output file

   mwrfits, float(startempset), outfile, hdr, /create


return
end
