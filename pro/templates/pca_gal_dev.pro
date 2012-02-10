;+
;
; NAME:
;  pca_gal_dev
;
; PURPOSE:
;  Developmental PCA of BOSS Galaxies
;
; WRITTEN:
;  bolton@utah 2010may
;
;-

pro mlpca_gal_dev, infile

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
maskid = [119,149,151,152,248,315,426,431,433, $
 213, 293, 168]-1
useit[maskid] = 0B
idx = where(useit and (strtrim(zans.objtype,2) eq 'GALAXY') $
 and (strtrim(zans.class,2) eq 'GALAXY'), ngal)
zans = zans[idx]

; Get some low-z galaxies as well:

infile = 'mlpcainput_gal_01.txt'

readcol, infile, fplate, fmjd, ffiberid, fz, format='L,L,L,F', comment='#'


i_lo = 0
i_hi = 125


fplate = fplate[i_lo:i_hi]
fmjd = fmjd[i_lo:i_hi]
ffiberid = ffiberid[i_lo:i_hi]
fz = fz[i_lo:i_hi]

plate = [fplate, zans.plate]
fiberid = [ffiberid, zans.fiberid]
mjd = [fmjd, zans.mjd]
z = [fz, zans.z]


nobj = n_elements(plate)

readspec, plate, fiberid, mjd=mjd, flux=flux, invvar=invvar, $
 loglam=loglam, zans=zans, /align, /silent

; B-spline to anything short of 3700Ang observed or 3400Ang rest:
; (this didn't really work right)
;idx = 10L * lindgen(50) + 5L

;plate = plate[idx]
;mjd = mjd[idx]
;fiberid = fiberid[idx]
;z = z[idx]

;@junkbuffer.pro

;dloglam = loglam[1] - loglam[0]
;bluefitlim = alog10(((1. + z) * 3400.) > 3700.)
;bkspace = 30 * dloglam
;flux_sm = flux
;for i = 0L, nobj-1 do begin & $
;  pmax = max(where(loglam lt bluefitlim[i], nwh)) & $
;  if (nwh gt 0) then begin & $
;    sset = bspline_iterfit(loglam[0:pmax], flux[0:pmax,i], nord=3, $
;     invvar=invvar[0:pmax,i], bkspace=bkspace, yfit=yfit, /silent) & $
;    flux_sm[0:pmax,i] = yfit & $
;  endif & $
;endfor

; Observed and rest-frame smoothing:

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
bluefitlim = alog10(((1. + z) * 3400.) > 3700.)
;bluefitlim = replicate(alog10(3700.), nobj)
redfitlim = replicate(alog10(9700.), nobj)
for i = 0L, nobj-1 do begin & $
  pmax = max(where(loglam lt bluefitlim[i], nwh)) & $
  if (nwh gt 0) then flux_sm[0:pmax,i] = flux_sm2[0:pmax,i] & $
  pmin = min(where(loglam gt redfitlim[i], nwh)) & $
  if (nwh gt 0) then flux_sm[pmin:npix-1,i] = flux_sm2[pmin:npix-1,i] & $
endfor


; Get everything aligned in rest-frame wavelength:
align_rest, in_loglam=loglam, in_flux=flux_sm, in_invvar=invvar_sm2, z=z, $
 out_loglam=rest_loglam, out_flux=rest_flux, out_invvar=rest_ivar


ngoodspec = total(rest_ivar gt 0., 2)
mingoodspec = 20
p_lo = min(where(ngoodspec ge mingoodspec))
p_hi = max(where(ngoodspec ge mingoodspec))
rest_flux = rest_flux[p_lo:p_hi,*]
rest_ivar = rest_ivar[p_lo:p_hi,*]
rest_loglam = rest_loglam[p_lo:p_hi]

; Cap the signal to noise so we don't get bossed around by low-z
; galaxies:

;snrmax = 10.
;snr = abs(rest_flux) * sqrt(rest_ivar)
;snrfact = snr / snrmax
;rest_ivar = rest_ivar / (snrfact > 1.)^2

; Determine practical convergence as when total chi^2 (in the well-fit
; scenario) changes by less than, say, 0.1 between the two steps of
; the MLPCA iteration process.
;ndof = total(rest_ivar gt 0)
;convtest = 0.1 / float(ndof)

; Somewhat heuristic:
convtest = 1.d-6

; Do the ML-PCA analysis:
galtempset = bolton_mlpca(objflux=rest_flux, objivar=rest_ivar, coeffs=coeffs, convtest=convtest, maxiter=25)

; Pad it out to 1800 angstroms rest with constant replication:
dloglam = rest_loglam[1] - rest_loglam[0]
n_new = round((rest_loglam[0] - alog10(1800.)) / dloglam)

new_loglam = [reverse(-dloglam * (dindgen(n_new) + 1.) + rest_loglam[0]), rest_loglam]


n_mean = 15L
mean_vec = total(galtempset[0:n_mean-1,*], 1) / float(n_mean)
new_tempset = [replicate(1.d0, n_new) # mean_vec, galtempset]

;res_gal = pca_solve(rest_flux, rest_ivar, acoeff=acoeff, nkeep=5)

; The models that come out of pca_solve don't seem to fit as well in
; many cases as those from bolton_mlpca.  Is this an under-convergence issue?


;; Generate model & inspect, if running via command line:
;model = galtempset # coeffs
;model = res_gal # acoeff
;i = -1L
;i++ & splot, rest_loglam, rest_flux[*,i] * (rest_ivar[*,i] gt 0) & soplot, rest_loglam, model[*,i], color=3

; Copied from pca_gal:
get_juldate, jd
mjdstr = STRING(LONG(jd-2400000L), FORMAT='(I5)')
outfile = 'spMLpcaGal-' + mjdstr + '-01.fits'
;outfile = 'spEigenGal-' + mjdstr + '.fits'
sxaddpar, hdr, 'OBJECT', 'GALAXY'
;sxaddpar, hdr, 'COEFF0', rest_loglam[0]
;sxaddpar, hdr, 'COEFF1', rest_loglam[1] - rest_loglam[0]
sxaddpar, hdr, 'COEFF0', new_loglam[0]
sxaddpar, hdr, 'COEFF1', new_loglam[1] - new_loglam[0]
sxaddpar, hdr, 'IDLUTILS', idlutils_version(), 'Version of idlutils'
sxaddpar, hdr, 'SPEC2D', idlspec2d_version(), 'Version of idlspec2d'
;mwrfits, float(galtempset), outfile, hdr, /create
mwrfits, float(new_tempset), outfile, hdr, /create
;mwrfits, float(res_gal), outfile, hdr, /create


return
end
