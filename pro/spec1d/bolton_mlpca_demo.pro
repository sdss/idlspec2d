;
; bolton_mlpca_demo
;
; Demonstrator script for Bolton's ML-PCA EigenGal implementation.
;
; bolton@utah 2010may07
;

; Require the following, because that's what I've vetted this on:
setenv, 'RUN2D=v5_4_9'
setenv, 'RUN1D=v5_4_9'

; The uber-deep plate:
plate = 3851
mjd = 55298

; Get all sensible galaxies:
readspec, plate, mjd=mjd, zans=zans
nobj = n_elements(zans)
useit = replicate(1B, nobj)
; Things to mask following inspection, as detailed above:
maskid = [119,149,151,152,248,315,426,431,433]-1
useit[maskid] = 0B
idx = where(useit and (strtrim(zans.objtype,2) eq 'GALAXY') $
 and (strtrim(zans.class,2) eq 'GALAXY'), ngal)
zans = zans[idx]

readspec, plate, zans.fiberid, mjd=mjd, flux=flux, synflux=synflux, $
 loglam=loglam, zans=zans, wave=wave, invvar=invvar, tsobj=tsobj, /align

; Get everthing aligned in rest-frame wavelength:
npix = (size(flux))[1]
dloglam = loglam[1] - loglam[0]
pixshift = round(alog10(1. + zans.z) / dloglam)
mp = max(pixshift)
totalpix = npix + mp
rest_flux = replicate(0., totalpix, ngal)
rest_ivar = replicate(0., totalpix, ngal)
rest_synflux = replicate(0., totalpix, ngal)

for i = 0L, ngal-1 do begin & $
  rest_flux[mp-pixshift[i]:mp-pixshift[i]+npix-1,i] = flux[*,i] & $
  rest_ivar[mp-pixshift[i]:mp-pixshift[i]+npix-1,i] = invvar[*,i] & $
  rest_synflux[mp-pixshift[i]:mp-pixshift[i]+npix-1,i] = synflux[*,i] & $
endfor

; Make the rest-frame loglam baseline:
loglam0_rest = loglam[0] - mp * dloglam
loglam_rest = loglam0_rest + dloglam * findgen(totalpix)

; Impose a minimum chi^2 cut on the extreme ends
; (somewhat arbitrarily set at "one sigma" of info per galaxy):
chi2vec = total(rest_flux^2 * rest_ivar, 2)
chi2filt = median(chi2vec, 15)
chi2min = float(ngal)
irange = minmax(where(chi2filt ge chi2min))
rest_flux = rest_flux[irange[0]:irange[1],*]
rest_ivar = rest_ivar[irange[0]:irange[1],*]
rest_synflux = rest_synflux[irange[0]:irange[1],*]
loglam_rest = loglam_rest[irange[0]:irange[1]]
totalpix = n_elements(loglam_rest)

; Cut it down again -- not enough galaxies covering
; the highest rest wavelengths in this sample.
; Including those sparsely covered ranges makes the
; algorithm go unstable.
;
; We can fix this with inclusion of more low-z galaxies.
p_lo = 0
p_hi = 4799
rest_flux = rest_flux[p_lo:p_hi,*]
rest_synflux = rest_synflux[p_lo:p_hi,*]
rest_ivar = rest_ivar[p_lo:p_hi,*]
loglam_rest = loglam_rest[p_lo:p_hi]

galtempset = bolton_mlpca(objflux=rest_flux, objivar=rest_ivar, coeffs=coeffs)

; Check these out:
modflux = galtempset # coeffs
i = -1L

; Up-arrow repeat:
i++ & splot, 10.^loglam_rest, rest_flux[*,i], yrange=[-1.,7.] & soplot, 10.^loglam_rest, modflux[*,i], color=3

; Copied from pca_gal:
get_juldate, jd
mjdstr = STRING(LONG(jd-2400000L), FORMAT='(I5)')
outfile = 'spMlpcaGal-' + mjdstr + '.fits'
sxaddpar, hdr, 'OBJECT', 'GALAXY'
sxaddpar, hdr, 'COEFF0', 10.^loglam_rest[0]
sxaddpar, hdr, 'COEFF1', dloglam
sxaddpar, hdr, 'IDLUTILS', idlutils_version(), 'Version of idlutils'
sxaddpar, hdr, 'SPEC2D', idlspec2d_version(), 'Version of idlspec2d'
mwrfits, galtempset, outfile, hdr, /create

