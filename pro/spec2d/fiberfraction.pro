;+
; NAME:
;   fiberfraction
;
; PURPOSE:
;   Calculates the fraction of light from an object with the specified Gaussian
;   fwhm that enters a fiber of the specified diameter when the object's centroid
;   is offset from the fiber center by the specified amount. All inputs should
;   be specified in the same units. The default diameter is the BOSS fiber size
;   in arcsecs, so requires that fwhm and offset also be specified in arcsecs.
;   Accuracy gets worse with increasing diameter/fwhm, but should be at least 0.5%
;   for fwhm > diameter/4 (or fwhm > 0.5 arcsec for BOSS fibers).
;
; CALLING SEQUENCE:
;  result = fiberfraction(psffwhm, offsets, fibersize)
;
; INPUTS:
;  psffwhm - The 2D Gaussian PSF FWHM.
;  offsets - A scalar or vector of offsets from the fiber center.
;  fibersize - The diameter of the fiber size.
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;  Returns an object matching dimension of the offsets input variable of the fraction
;  of light that enters a fiber.
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; REVISION HISTORY:
;   25-July-2013  Written by Daniel Margala (dmargala@uci.edu), UC Irvine.
;-
;------------------------------------------------------------------------------

function fiberfraction, fwhm, offsets, diameter
    if (fwhm le 0) then begin
        splog, "Invalid fwhm <= 0."
        return, -1
    endif
    index = where(offsets lt 0, count)
    if (count ne 0) then begin
        splog, "Invalid offset < 0."
        return, -1
    endif
    if (diameter le 0) then begin
        splog, "Invalid diameter <= 0."
        return, -1
    endif
    if (diameter gt 4.0*fwhm) then begin
        splog, "diameter > 4*fwhm not implemented."
        return, -1
    endif

    ; Calculate dimensionless ratios
    sigma = fwhm/2.3548200450309493820 ; constant is 2*sqrt(2*log(2))
    t = offsets/sigma
    ss = diameter/(2.*sigma)
    tSqby2 = t*t/2.
    ssSqby2 = ss*ss/2.

    ; Use a series expansion the BesselI[0,x] appearing in the radial inegral.
    lastSum = replicate(0., n_elements(offsets))
    sum = 1-exp(-ssSqby2)
    prod = replicate(1., n_elements(offsets))
    k = 0
    while (++k lt 1000) do begin
        lastSum = sum
        prod *= tSqby2
        term = prod*igamma(k+1,ssSqby2)/gamma(k+1)
        sum += term 
        index = where(abs(term) lt 1e-2*sum, count)
        if (count eq n_elements(offsets)) then break 
    endwhile

    return, sum*exp(-tSqby2)
end
