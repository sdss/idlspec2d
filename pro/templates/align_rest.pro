;+
;
; NAME: align_rest
;
; PURPOSE:
;  Align spectra in rest wavelength given redshift.
;
; COMMENTS:
;  Assumes constant log-lam binning, rounds to nearest pixel.
;
; WRITTEN:
;  bolton@utah 2010may
;
;-

pro align_rest, in_loglam=in_loglam, in_flux=in_flux, in_invvar=in_invvar, z=z, $
 out_loglam=out_loglam, out_flux=out_flux, out_invvar=out_invvar

npix = (n_elements(in_loglam))
nobj = (size(in_flux))[2]

dloglam = in_loglam[1] - in_loglam[0]
pixshift = round(alog10(1.d0 + z) / dloglam)
maxp = max(pixshift)
minp = min(pixshift)
totalpix = npix + maxp - minp

out_flux = make_array(totalpix, nobj, type=size(in_flux, /type), value=0)
out_invvar = make_array(totalpix, nobj, type=size(in_invvar, /type), value=0)

for i = 0L, nobj-1 do begin
  out_flux[maxp-pixshift[i]:maxp-pixshift[i]+npix-1,i] = in_flux[*,i]
  out_invvar[maxp-pixshift[i]:maxp-pixshift[i]+npix-1,i] = in_invvar[*,i]
endfor

; Make the rest-frame loglam baseline:
loglam0_rest = in_loglam[0] - maxp * dloglam
out_loglam = loglam0_rest + dloglam * findgen(totalpix)

return
end
