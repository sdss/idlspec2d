;+
; NAME:
;   k-correct
; PURPOSE:
;   Take 5-band photometry and redshift and return k-corrected magnitude.
; COMMENTS:
; CALLING SEQUENCE:
;   restphot= k-correct(obsphot,z)
; INPUTS:
;   obsphot   - 5xN array of SDSS photometry in observed frame
;   z         - N-vector of redshifts
; OPTIONAL INPUTS:
;   obserr    - 5xN array of uncertainties in obsphot
; KEYWORDS:
; OUTPUTS:
;   restphot  - 5xN array of SDSS photometry in rest frame
; OPTIONAL OUTPUTS:
;   resterr   - 5xN array of propagated uncertainties in restphot
; PROCEDURES CALLED:
; BUGS:
;   MAKES NO USE OF ERRORS.
;   GIVES CRAP WHEN EXTRAPOLATING.
;   This uses an incredibly stupid interpolation scheme.
; REVISION HISTORY:
;   2000-Jun-28  Written by Hogg (IAS)
;-
function k_correct, obsphot,z,obserr=obserr,resterr=resterr

; get bandpass information
  bandpassinfo, indgen(5),index=band,name=name,wave=wave,zero=zero, $
                redindex=ri,blueindex=bi
  nband= 5

; get array sizes
  nz= n_elements(z)

; start loop over redshifts
  loglam= alog10(wave)
  logflux= alog10(zero)-0.4*obsphot
  emit= dblarr(5,nz)
  obse= dblarr(5,nz)
  logfluxk= dblarr(5,nz)
  restphot= dblarr(5,nz)
  for i=0,nz-1 do begin

; make log lambda points
    emit[band,i]= loglam[band]-alog10(1.0+z[i]) ; emitted
    obse[band,i]= loglam[band]                  ; observed

; do spline interpolation
    y1= spl_init(emit[band],logflux[band,i])
    logfluxk[band,i]= spl_interp(emit[band],logflux[band,i],y1,obse[band,i])

; convert back to magnitudes
    restphot[band,i]= 2.5*alog10(zero)-2.5*logfluxk[band,i]
  endfor

  return, restphot
end

