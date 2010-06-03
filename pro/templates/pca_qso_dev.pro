;+
;
; NAME: pca_qso_dev
;
; PURPOSE: Developmental PCA routine for BOSS QSOs
;
; WRITTEN: bolton@utah 2010may
;
;-

pro pca_qso_dev, infile

; infile = 'mlpcainput_qso_01.txt'
; infile = 'mlpcainput_qso_02.txt'
infile = 'mlpcainput_qso_combined.txt'

readcol, infile, plate, mjd, fiberid, z, format='L,L,L,F', comment='#'

nobj = n_elements(plate)

readspec, plate, fiberid, mjd=mjd, flux=flux, invvar=invvar, $
 loglam=loglam, zans=zans, /align, /silent

; res_qso = pca_solve(flux, invvar, loglam, z, nkeep=5, acoeff=acoeff)

; pca_solve Doing ths straight away gives a chunky output...

; The smoothing that we're using:
flux_sm2 = flux
invvar_sm2 = invvar
nsmooth = 31
npix = n_elements(loglam)
for i = 0L, nobj-1 do begin & $
  numer = smooth(flux[*,i] * (invvar[*,i] gt 0), nsmooth) & $
  denom = smooth(float(invvar[*,i] gt 0), nsmooth) > 1. & $
  flux_sm2[*,i] = numer / denom & $
  gmin = min(where(invvar[*,i] gt 0.)) + nsmooth/2 & $
  gmax = max(where(invvar[*,i] gt 0.)) - nsmooth/2 & $
  invvar_sm2[0:gmin,i] = 0. & $
  invvar_sm2[gmax:npix-1,i] = 0. & $
endfor

flux_sm = flux
bluefitlim = alog10(((1. + z) * 800.) > 3600.)
;bluefitlim = replicate(alog10(3700.), nobj)
;redfitlim = replicate(alog10(9700.), nobj)
for i = 0L, nobj-1 do begin & $
  pmax = max(where(loglam lt bluefitlim[i], nwh)) & $
  if (nwh gt 0) then flux_sm[0:pmax,i] = flux_sm2[0:pmax,i] & $
;  pmin = min(where(loglam gt redfitlim[i], nwh)) & $
;  if (nwh gt 0) then flux_sm[pmin:npix-1,i] = flux_sm2[pmin:npix-1,i] & $
endfor

; Get everything aligned in rest-frame wavelength:
align_rest, in_loglam=loglam, in_flux=flux_sm, in_invvar=invvar_sm2, z=z, $
 out_loglam=rest_loglam, out_flux=rest_flux, out_invvar=rest_ivar

; Make sure we have enough to work with:
;ngoodspec = total(rest_ivar gt 0., 2)
;mingoodspec = 10
;p_lo = min(where(ngoodspec ge mingoodspec))
;p_hi = max(where(ngoodspec ge mingoodspec))
;rest_flux = rest_flux[p_lo:p_hi,*]
;rest_ivar = rest_ivar[p_lo:p_hi,*]
;rest_loglam = rest_loglam[p_lo:p_hi]

; Determine practical convergence as when total chi^2 (in the well-fit
; scenario) changes by less than, say, 0.1 between the two step of
; the MLPCA iteration process.
;ndof = total(rest_ivar gt 0)
;convtest = 0.1 / float(ndof)

;qsotempset = bolton_mlpca(objflux=rest_flux, objivar=rest_ivar, coeffs=coeffs, convtest=convtest, maxiter=10)

res_qso = pca_solve(rest_flux, rest_ivar, nkeep=5, acoeff=acoeff)

;; Generate model & inspect, if running via command line:
;model = qsotempset # coeffs
;model = res_qso # acoeff
;i = -1L
;i++ & splot, rest_loglam, rest_flux[*,i] * (rest_ivar[*,i] gt 0) & soplot, rest_loglam, model[*,i], color=3

; Copied from pca_gal:
get_juldate, jd
mjdstr = STRING(LONG(jd-2400000L), FORMAT='(I5)')
outfile = 'spMLpcaQSO-' + mjdstr + '-01.fits'
;outfile = 'spEigenQSO-' + mjdstr + '.fits'
sxaddpar, hdr, 'OBJECT', 'QSO'
sxaddpar, hdr, 'COEFF0', rest_loglam[0]
sxaddpar, hdr, 'COEFF1', rest_loglam[1] - rest_loglam[0]
sxaddpar, hdr, 'IDLUTILS', idlutils_version(), 'Version of idlutils'
sxaddpar, hdr, 'SPEC2D', idlspec2d_version(), 'Version of idlspec2d'
;mwrfits, float(qsotempset), outfile, hdr, /create
mwrfits, float(res_qso), outfile, hdr, /create

return
end
