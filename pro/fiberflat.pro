;+
; NAME:
;   fiberflat
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   fflat = fiberflat( flat_flux, wset )
;
; INPUTS:
;   flat_flux  - Array of extracted flux from a flat-field image [Nrow,Ntrace]
;   wset       - Trace set for wavelengths
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   fflat      - Array of flat-field flat-field vectors for each fiber
;                that remove relative flat-field variations as a function
;                of wavelength between fibers [Nrow, Ntrace]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The user should first "flat-field" the input array to take out
;   pixel-to-pixel variations.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   14-Oct-1999  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------

function fiberflat, flat_flux, wset

   ntrace = (size(flat_flux))[2]

   ; Compute the wavelengths for the extracted spectra
   traceset2xy, wset, xx, yy
   waves = 10^yy
   ny = (size(waves))[1]

   ; Find the largest possible wavelength range.
   ; Also find the minimum wavelength separation between pixels.
   wavemin = min(waves)
   wavemax = max(waves)
   if (waves[0,0] LT waves[ny-1,0]) then $
    deltaw = min( waves[1:ny-1,*] - waves[0:ny-2,*] ) $ ; ascending
   else $
    deltaw = min( waves[0:ny-2,*] - waves[1:ny-1,*] )   ; descending

   ; Linearly interpolate all of the flat-field vectors onto a common
   ; wavelength scale.  Set the mask array =1 for anywhere on a fiber
   ; that's missing wavelength coverage.
   wtemp = wavemin + findgen(fix(wavemax-wavemin)/deltaw)
   flattemp = fltarr(N_elements(wtemp), ntrace)
   masktemp = fltarr(N_elements(wtemp), ntrace)
   for i=0, ntrace-1 do begin
      flattemp[*,i] = interpol(flat_flux[*,i], waves[*,i], wtemp)
      masktemp[*,i] = wtemp LT min(waves[*,i]) OR wtemp GT max(waves[*,i])
flattemp[*,i] = smooth(flattemp[*,i],100)
      fullbkpt = efc(wtemp, flattemp[*,i], coeff, nbkpts=100)
      yfit = bvalu(wtemp, fullbkpt, coeff)
;plot, wtemp, flattemp[*,i]
;djs_oplot, wtemp, yfit, color='red'
;plot,wtemp,flattemp[*,i]/yfit,yr=[0.9,1.1],ps=3
   endfor

   ; Find wavelengths where all fibers have data
   mashmask = total(masktemp,2)
   wgood = where(mashmask EQ ntrace)

   mm = total(flattemp,2) / ntrace

   for i=0, ntrace-1 do $
    flattemp[*,i] = flattemp[*,i] / mm

   ; Linearly interpolate back to the original wavelength scale for each fiber
   fflux = fltarr(ny, ntrace)
   for i=0, ntrace-1 do $
      fflux[*,i] = interpol(flattemp[*,i], wtemp, waves[*,i])

; ???
stop
   return, fflat
end
;------------------------------------------------------------------------------
