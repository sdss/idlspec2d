; given a spectrum, dereddenning it using the SFD map and Milky Way extinction curve (CCM)
; to correct for Galactic extinction
; NB, the flux errors should also be renormalized
; 

pro dereddening_spec, wave, flux, err=err, ra=ra, dec=dec, $
    dered_flux=dered_flux, dered_err=dered_err, eb_v=eb_v, A_lam=A_lam

rv = 3.1  ; Rv for Milky Way

A_lam = mgb_extinction(ra,dec,rv,wave,0, eb_v=eb_v)
scale = 10.0D^(0.4*A_lam)
dered_flux = flux*scale

if keyword_set(err) then dered_err = err*scale

return
end
