;
;+
; NAME:
;	masklines.pro

;
; PURPOSE:
;	Identify and mask emission/absorption features in a spectrum
;
;
; CALLING SEQUENCE:
;
;	pixel_mask = MASKLINES(Flux, Wave, Linefile, Newflux, /SILENT)
;
;
; INPUTS:
;	Flux:       spectrum
;       Wave:       rest wavelength array, Angstroms (same length as Flux)
;       Linefile:   Path and name of file containing array of
;                   wavelengths to be masked.
;
; KEYWORD PARAMETERS: 
;       SILENT:    stop program writing warnings
;
; OUTPUTS:
;	This function returns an array the same length of flux where 1
;	identifies a pixel with a line, and 0 a pixel without a line.
;
; OPTIONAL OUTPUT:
;       Newflux:    The pixels identified as containing lines are
;                   replaced with the median of 50 nearby pixels. NOTE
;                   variable newflux must be defined as non-zero proir
;                   to calling program. 
;                  
;	
; MODIFICATION HISTORY: Vivienne Wild, vw@ast.cam.ac.uk, 05/01/05
;
;****************************************************************************************


function MASKLINES, flux, wave, linefile, newflux, silent=silent

if N_ELEMENTS(wave) ne N_ELEMENTS(flux) then message, 'flux and wave must be same size'

;read in line file

OPENR, lun, linefile, /get_lun
nlines = 0
READF, lun, nlines  
emlines = FLTARR(2,nlines)
READF, lun, emlines
FREE_LUN, lun

;create pixel mask
nbin = N_ELEMENTS(wave)
pixel_mask = FLTARR(nbin)
for i=0, nlines-1 do begin

    index = where(wave gt emlines[0,i] and wave lt emlines[1,i],count) 
    if count gt 0 then pixel_mask[index] = 1.

endfor

;warnings
junk = WHERE(pixel_mask eq 1, n)
if n eq 0 and NOT(KEYWORD_SET(SILENT)) then print,'masklines.pro: no em lines masked'

;mask flux array if required
if N_ELEMENTS(newflux) ne 0 and n ne 0 then begin
    
    goodpix = FINDGEN(nbin)
    badpix = FINDGEN(nbin)
    index = WHERE(pixel_mask eq 1,compl=compl)
    goodpix = goodpix[compl]
    badpix = badpix[index]
    newflux =flux

    for i = 0, N_ELEMENTS(badpix) -1 do begin
        up = WHERE(goodpix gt badpix[i]) ;non-line pixels above line
        low = WHERE(goodpix lt badpix[i],count) ;non-line pixels below line
        if up[0] ne -1 then begin
            flag1 = 1
            if n_elements(up) gt 25 then index1 = up[0:24] else index1 = up ;+/-25 pixels chosen to match typical filter scale 
        endif else flag1 = 0    ;line sits at far red of spectrum
        if low[0] ne -1 then begin
            flag2 = 1
            if n_elements(low) gt 25 then index2 = low[count-25:count-1] else index2 = low
        endif else flag2 = 0    ;line sits at far blue of spectrum

        if flag1 and flag2 then index = [index1,index2]
        if flag1 and NOT(flag2) then index = [index1]
        if flag2 and NOT(flag1) then index = [index2]
        if NOT(flag1) and NOT(flag2) then begin
            if NOT(KEYWORD_SET(SILENT)) then print,'something wrong with line mask'
            continue
        endif

        newflux[badpix[i]] = MEDIAN(flux[goodpix[index]])
    endfor
        
endif else if N_ELEMENTS(newflux) ne 0 and n eq 0 then newflux = flux

return, pixel_mask

end
