;
;+
; NAME:
;	medfilt.pro

;
; PURPOSE:
;	Median filter a spectrum.
;
;
; CALLING SEQUENCE:
;
;	newflux = MEDFILT(Flux, Ind_filter, Nfilt)
;
;
; INPUTS:
;	Flux:       spectrum
;       Ind_filter: pixels which require median filtering (allows proper
;                   filtering at ends where available)
;       Nfilt:      size of filter in pixels
;
; OUTPUTS:
;	This function returns the median filtered spectrum.
;	
; NOTES: 
;       Use Ind_filter = findgen(n_elements(flux)) to median filter
;       whole input spectrum.
;
; MODIFICATION HISTORY: Vivienne Wild, vw@ast.cam.ac.uk, 05/01/05
;
;****************************************************************************************

FUNCTION MEDFILT, Flux, Ind_filter, Nfilt

nhalf = nfilt/2
nbin = n_elements(ind_filter)
npix = n_elements(flux)

med_spec = fltarr(nbin)
for j=0,nbin-1 do begin

    start = j-nhalf+ind_filter[0]
    finish = j+nhalf+ind_filter[0]

    if start lt 0 then start = 0
    if finish gt npix then finish=npix

    med_spec[j] = median(flux[start:finish-1],/even)
endfor

return, med_spec

END
