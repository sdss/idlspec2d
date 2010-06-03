;+
;
; NAME: mlpca_cvstar
;
; PURPOSE:
;  Generate ML-PCA templates for CV stars
;
; bolton@utah 2010may
;
;-


pro mlpca_cvstar, infile

infile = 'mlpcainput_star.txt'

readcol, infile, plate, mjd, fiberid, z, subc, format='L,L,L,F,A'

subc = strtrim(subc)

; Exclude non-CVs:
wh = where(subc eq 'CV')
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


cvstartempset = bolton_mlpca(objflux=rest_flux, objivar=rest_ivar, coeffs=coeffs, nbasis=2)

   ;----------
   ; Construct header for output file

get_juldate, jd
mjdstr = STRING(LONG(jd-2400000L), FORMAT='(I5)')
outfile = 'spMLpcaCVstar-' + mjdstr + '.fits'
sxaddpar, hdr, 'OBJECT', 'STAR_CV'
sxaddpar, hdr, 'COEFF0', rest_loglam[0]
sxaddpar, hdr, 'COEFF1', rest_loglam[1] - rest_loglam[0]


   ;----------
   ; Write output file

   mwrfits, float(cvstartempset), outfile, hdr, /create


return
end
